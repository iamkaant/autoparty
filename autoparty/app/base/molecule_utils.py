
import pandas as pd

import os
import glob
import zipfile
import datetime

import json
import numpy as np
import rdkit.Chem as Chem
from rdkit.Chem import AllChem, Lipinski, rdMolDescriptors, ChemicalFeatures
from rdkit import RDConfig

from sqlalchemy import desc, and_, or_, func
from sqlalchemy.sql import exists, false
from sqlalchemy.orm.exc import NoResultFound

from app import celery, db
from app.base.database import Molecule, Grade, Prediction, HPRunSetting, Model
from app.base.defaults import *

orderby_dict = {"score": [Molecule.score], "uncertainty": [desc(Prediction.uncertainty)], 
				"prediction": [Prediction.prediction, Prediction.uncertainty], 
				"disagreement": [desc(Prediction.error)]}

LUNA_IFP_LENGTH = 4096
LIGAND_MORGAN_BITS = 2048
LIGAND_DESC_LEN = 10
HBOND_INTERACTION_TYPES = {"Hydrogen bond", "Weak hydrogen bond", "Water-bridged hydrogen bond"}

_FEATURE_FACTORY = None


def _get_feature_factory():
	global _FEATURE_FACTORY
	if _FEATURE_FACTORY is None:
		fdef = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
		_FEATURE_FACTORY = ChemicalFeatures.BuildFeatureFactory(fdef)
	return _FEATURE_FACTORY


def _get_ligand_hbond_feature_counts(mol):
	"""Count ligand donor/acceptor feature groups using RDKit feature definitions."""
	factory = _get_feature_factory()
	feats = factory.GetFeaturesForMol(mol)

	donor_count = 0
	acceptor_count = 0
	for feat in feats:
		family = feat.GetFamily()
		if family == "Donor":
			donor_count += 1
		elif family == "Acceptor":
			acceptor_count += 1

	return donor_count, acceptor_count


def _normalize_inters(inters):
	if inters is None:
		return []
	if isinstance(inters, str):
		try:
			inters = json.loads(inters)
		except Exception:
			return []
	if isinstance(inters, str):
		try:
			inters = json.loads(inters)
		except Exception:
			return []
	return inters if isinstance(inters, list) else []


def _ligand_group_key(group):
	centroid = group.get('centroid', [])
	if len(centroid) == 3:
		return (round(float(centroid[0]), 3), round(float(centroid[1]), 3), round(float(centroid[2]), 3))

	atom_names = []
	for atom in group.get('atoms', []):
		name = atom.get('name')
		if isinstance(name, list) and len(name) > 0:
			atom_names.append(str(name[0]).strip())
		elif isinstance(name, str):
			atom_names.append(name.strip())
	if len(atom_names) == 0:
		return None
	return tuple(sorted(atom_names))


def _count_ligand_hbond_features_involved(inters):
	"""Count unique ligand donor/acceptor feature groups participating in H-bonds."""
	inters = _normalize_inters(inters)
	involved_donors = set()
	involved_acceptors = set()

	for inter in inters:
		if inter.get('type') not in HBOND_INTERACTION_TYPES:
			continue
		for grp in (inter.get('src_grp', {}), inter.get('trgt_grp', {})):
			compounds = grp.get('compounds', [])
			if not any(comp.get('chain') == 'z' for comp in compounds):
				continue
			group_key = _ligand_group_key(grp)
			if group_key is None:
				continue
			features = set(grp.get('features', []))
			if 'Donor' in features or 'WeakDonor' in features:
				involved_donors.add(group_key)
			if 'Acceptor' in features:
				involved_acceptors.add(group_key)

	return len(involved_donors), len(involved_acceptors)


def _ligand_feature_vector(smiles, inters=None):
	"""Ligand-only structural features for moieties and H-bond tendencies."""
	mol = Chem.MolFromSmiles(smiles) if smiles else None
	if mol is None:
		return np.zeros(LIGAND_MORGAN_BITS + LIGAND_DESC_LEN, dtype=np.float32)

	morgan = np.array(AllChem.GetMorganFingerprintAsBitVect(
		mol,
		radius=2,
		nBits=LIGAND_MORGAN_BITS,
		useFeatures=True
	), dtype=np.float32)

	num_hbd, num_hba = _get_ligand_hbond_feature_counts(mol)
	involved_d, involved_a = _count_ligand_hbond_features_involved(inters)
	stranded_d = max(0.0, float(num_hbd) - float(involved_d))
	stranded_a = max(0.0, float(num_hba) - float(involved_a))

	desc = np.array([
		float(len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 35])),  # bromine count
		float(len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 17])),  # chlorine count
		float(len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 53])),  # iodine count
		float(rdMolDescriptors.CalcNumSaturatedRings(mol)),
		float(rdMolDescriptors.CalcNumAromaticRings(mol)),
		float(num_hbd),
		float(num_hba),
		float(stranded_d),
		float(stranded_a),
		float(rdMolDescriptors.CalcFractionCSP3(mol)),
	], dtype=np.float32)

	return np.concatenate((morgan, desc), axis=0)


def _combine_ifp_and_ligand_features(ifp_json, smiles, inters=None, ifp_length=LUNA_IFP_LENGTH):
	ifp = _unpack_fp(ifp_json, fp_length=ifp_length)
	lig = _ligand_feature_vector(smiles, inters=inters)
	return np.concatenate((ifp, lig), axis=0)

# no reason for these to be tasks really, they're synchronous
def get_molecule_by_id(mol_id, return_grade = False, return_pred = False, party_id = -1):
	"""
	Returns single molecule by id and grade if requested and available.
	"""
	mol = Molecule.query.get(mol_id)

	if return_grade:
		try:
			grade = Grade.query.filter_by(mol_id = mol_id
				).filter_by(hp_settings_id = party_id).one()
		except:
			grade = None
		if not return_pred:
			return mol, grade
	if return_pred:
		try:
			pred = Prediction.query.filter_by(mol_id = mol_id
				).filter_by(hp_settings_id = party_id).one()
		except:
			pred = None
		if not return_grade:
			return mol, pred
		else:
			return mol, pred, grade
			
	return mol


def get_ordered_molecules(screen_id, party_id, name = None, orderby = "score", mode = "annotate", modetime = None,
		limit = None, offset = None):

	#all_mols = db.session.query(Molecule).all()

	"""
	Return molecules that do NOT already have grades, sorted by requested method
	"""
	if name is not None: # if molecule is requested by name, exclude nothing
		exclude = false()

	else: # either in annotation or review mode, exclude mols accordingly
		if mode == "annotate": # exclude all molecules with grades
			exclude = exists().where(Molecule.id == Grade.mol_id
				).where(Grade.hp_settings_id == party_id)
		elif mode == "review": # exclude molecules with NEW grades
			exclude = exists().where(Molecule.id == Grade.mol_id
				).where(Grade.timestamp > modetime)

	# base query - Molecules that belong to this run and don't meet the conditions above
	query = db.session.query(Molecule
		).filter(Molecule.run_id == screen_id # limit this sreen
		).filter(~exclude) 

	if name is not None: #check if we need a name filter
		# add name filter to query
		query = query.filter(Molecule.name.ilike(f"%{name}%"))

	if orderby != "score" and orderby != "name": # need to get predictions to order - uncertainty, prediction 
		query = query.outerjoin(Prediction
			).filter(Prediction.hp_settings_id == party_id)

	#unordered_mols = query.offset(offset).limit(limit).all()
	
	# decide ordering - either score, uncertainty, prediction, or disagreement(error)
	if orderby != 'name':
		query = query.order_by(*orderby_dict[orderby])
	#total = query.count()
	total = db.session.query(func.count(Molecule.id)).filter(
		Molecule.run_id == screen_id).filter(~exclude).scalar()

	mols = query.offset(offset).limit(limit).all()

	grades = get_relations(party_id, mols, "annotations", Grade)
	preds = get_relations(party_id, mols, "uncertains", Prediction)

	return mols, grades, preds, total

def get_relations(party_id, mols, relation, dbtable):
	relations = []
	for mol in mols:
		try:
			if relation == "annotations": dbrelat = mol.annotations 
			else: dbrelat = mol.uncertains

			relat = dbrelat.filter(dbtable.hp_settings_id == party_id).one_or_none()

		except:
			relat = None
		relations.extend([relat])
	return relations

def get_predictions(party_id, sort_by = "prediction"):
	if sort_by not in ['uncertainty', 'prediction']:
		sort_by = 'prediction'

	order = [Prediction.prediction, Prediction.uncertainty]
	if sort_by != 'prediction':
		order = order[::-1]

	preds = Prediction.query.filter_by(hp_settings_id = party_id).order_by(*order).all()
	return preds

def _unpack_fp(ifp_json, fp_length = LUNA_IFP_LENGTH):
	fp = np.zeros(fp_length)
	for key in ifp_json:
		fp[int(key)] = ifp_json[key]
	return fp

def get_num_grades(party_id):
	return Grade.query.filter_by(hp_settings_id = party_id).count()

def get_grades_for_training(party_id, format_dataframe = True, fp_col = 'fp', label_col = 'label', fp_length = LUNA_IFP_LENGTH):
	"""
	Gets all the grades for a run, formats into expected dataframe for FingerprintDataset if requested
	"""

	grades = Grade.query.filter_by(hp_settings_id = party_id).all()

	if format_dataframe:
		data = [
			(
				int(grade.mol_id),
				_combine_ifp_and_ligand_features(json.loads(grade.ifp), grade.molecule.smi, grade.molecule.inters, fp_length),
				grade.grade,
			)
			for grade in grades
		]
		df = pd.DataFrame(data, columns = ["id", fp_col, label_col])
		df.set_index("id", inplace=True)
		return df
	return grades

def get_molecules_for_predicting(mol_ids, format_dataframe = True, fp_col = "fp"):
	molecules = Molecule.query.filter(Molecule.id.in_(mol_ids)).all()

	if format_dataframe:
		data = [
			(
				int(mol.id),
				_combine_ifp_and_ligand_features(mol.ifp.counts, mol.smi, mol.inters),
			)
			for mol in molecules if mol.ifp is not None
		] # filter out mols with no ifp if present
		df = pd.DataFrame(data, columns = ["id", fp_col])
		df.set_index("id", inplace=True)
		return df.reindex(mol_ids).dropna()
	return molecules

# putting these here cause I dont know where else they fit
def prediction_csv(user_id, screen_id, party_id):
	"""
	Creates output csv file containing moelcule ids, names, grades, predictions, smiles
	"""
	prefix = "-".join([f"{val_id}{val}"for val_id, val in zip(['u', 's', 'p'], [user_id, screen_id, party_id])])

	# get existing grades - write these molecules first, and there are fewer of them so less db queries that doing it the other way around
	data = {}

	grades = get_grades_for_training(party_id, format_dataframe = False)
	for grade in grades:
		mol = grade.molecule
		data[mol.id] = [mol.name, mol.score, grade.grade, mol.smi]

	preds = get_predictions(party_id)

	for pred in preds:
		mol = pred.molecule
		if mol.id not in data.keys():
			data[mol.id] = [mol.name, mol.score, '',  mol.smi, pred.prediction, pred.uncertainty]
		else:
			data[mol.id] += [pred.prediction, pred.uncertainty]

	data_df = pd.DataFrame.from_dict(data, orient='index',
		columns = ['name', 'score', 'grade', 'smi', 'prediction', 'uncertainty'])

	data_df.to_csv(f"{UPLOAD_FOLDER}/{prefix}_predictions.csv", index_label='id')
	return True


def update_history(user_id, screen_id, hp_settings_id, update_info):

	new_model = Model(
		user_id = user_id,
		run_id = screen_id,
		hp_run_id = hp_settings_id,
		train_loss = update_info['train_loss'],
		val_loss = update_info['val_loss'],
		num_grades = update_info['train_size'])

	db.session.add(new_model)
	db.session.commit()

def recover_history(hp_settings_id):
	history_dict = {}
	history_dict["Training"] = []
	history_dict["Validation"] = []
	history_dict["Time"] = []
	try:
		# flip dict, seperate training, validation, num_sgrades
		models = Model.query.filter_by(hp_run_id = hp_settings_id
			).order_by(Model.timestamp).all()

		for model in models:
			history_dict['Training'].append(model.train_loss)
			history_dict['Validation'].append(model.val_loss)
			history_dict['Time'].append(str(model.timestamp).split('.')[0].replace(' ', '<br>')) #html syntax
		return history_dict

	except Exception as e:
		print(e)
		return {}
		