# AutoParty Run and Debug Guide (macOS)

This guide captures a working setup for running AutoParty from source on macOS (including Apple Silicon) and a checklist for debugging common startup/runtime issues.

## 1) Prerequisites

- Conda (Anaconda/Miniconda)
- Homebrew
- Redis server

Install Redis (if needed):

### D) Water-mediated H-bonds not showing
```bash
brew install redis
```

## 2) Environment Setup

From the repository root:

```bash
cd /Users/ak87/ownCloud/owncloud/Calculations/TMP/AutoParty
```

Create the project environment.

For Apple Silicon (M1/M2/M3), use x86_64 conda packages for compatibility:

```bash
CONDA_SUBDIR=osx-64 mamba env create -f autoparty/autoparty-env.yml
```

If you do not use `mamba`, use:

```bash
CONDA_SUBDIR=osx-64 conda env create -f autoparty/autoparty-env.yml
```

Install LUNA into the environment (required by AutoParty):

### E) Config file precedence (common confusion)
```bash
conda run -n autoparty-env pip install luna
```

## 3) Start Services
### F) Hydrogens visible in files but not in 3D viewer

AutoParty needs three pieces:

1. Redis
2. Celery worker
3. Flask app

Start Redis:

```bash
brew services start redis
```

Run Celery (terminal A):
### G) After changing defaults/config parsing code

```bash
cd /Users/ak87/ownCloud/owncloud/Calculations/TMP/AutoParty/autoparty
conda run -n autoparty-env bash start_celery.sh 4
```

Run Flask app (terminal B):

```bash
cd /Users/ak87/ownCloud/owncloud/Calculations/TMP/AutoParty/autoparty
conda run -n autoparty-env bash start_autoparty.sh
```

Open AutoParty:

- `http://127.0.0.1:5000`

## 4) Port Conflict on 5000 (Common on macOS)

If you see:

- `Address already in use`
- `Port 5000 is in use`

Run AutoParty on another port:

```bash
cd /Users/ak87/ownCloud/owncloud/Calculations/TMP/AutoParty/autoparty
conda run -n autoparty-env python manage.py 5001
```

Then open:

- `http://127.0.0.1:5001`

Find what is using 5000:

```bash
lsof -nP -iTCP:5000 -sTCP:LISTEN
```

On macOS, this is often Control Center / AirPlay Receiver.

## 5) Quick Health Checks

Check Redis:

```bash
redis-cli ping
```

Expected:

```text
PONG
```

Check web app response:

```bash
curl -sS -m 5 -o /dev/null -w '%{http_code}\n' http://127.0.0.1:5001
```

A `302` or `403` still confirms the app is up (depends on route/auth state).

## 6) Debugging Checklist

### A) `ModuleNotFoundError: No module named 'luna.MyBio'`

Cause:
- Newer LUNA package layout moved modules from `luna.MyBio.*` to `luna.pdb.*`.

Fix applied in this repo:
- Compatibility imports were added in [autoparty/app/base/luna_utils.py](autoparty/app/base/luna_utils.py) so both old and new LUNA layouts work.

### B) Flask starts but Celery tasks fail

Checks:

1. Confirm Redis is running (`redis-cli ping`)
2. Confirm Celery terminal is running without immediate exit
3. Inspect Celery log output:

```bash
cd /Users/ak87/ownCloud/owncloud/Calculations/TMP/AutoParty/autoparty
tail -f outputs/celery.logs
```

### C) Environment/package mismatch

Validate key imports in the conda env:

```bash
conda run -n autoparty-env python -c "import torch, rdkit, openbabel, luna, flask, celery, redis; print('ok')"
```

If this fails, reinstall into the same env:

```bash
conda run -n autoparty-env pip install --upgrade pip
conda run -n autoparty-env pip install luna
```

### D) Water-mediated H-bonds not showing

Symptoms:
- Structural waters are visible in the 3D viewer, but no water H-bond annotations/lines appear.

Important behavior:
- Interaction JSON is precomputed by Celery when a screen is processed.
- Changing LUNA configs does not retroactively update existing screens.

Required settings (in effective `*_luna.cfg`):

```ini
add_h2o_pairs_with_no_target = True
lazy_comps_list = NH3,NH4
```

Verify the latest generated config actually used by the run:

```bash
cd /Users/ak87/ownCloud/owncloud/Calculations/TMP/AutoParty/autoparty
ls -lt outputs/*_luna.cfg | head -n 5
sed -n '1,80p' outputs/u1-sXX_luna.cfg
```

If settings are wrong, then:
1. Update defaults in `autoparty/defaults/LUNA_default.cfg` (or upload custom LUNA configs at screen creation).
2. Restart Flask and Celery.
3. Re-upload/reprocess the screen so interactions are recalculated.

### E) Config file precedence (common confusion)

- Files placed in `autoparty/inputs/` are input data by convention and are not automatically used as uploaded LUNA config overrides.
- For custom LUNA behavior, provide config files through the UI upload flow when creating a new screen.
- Otherwise AutoParty uses defaults from `autoparty/defaults/`.

### F) Hydrogens visible in files but not in 3D viewer

Known required conditions in this repo:
- PDB must be parsed with `keepH: true` in viewer setup.
- Ligands must be loaded with explicit H preserved (`removeHs=False` in RDKit loaders).
- Interacting-residue H display is style-driven and separate from interaction calculation.

Quick sanity checks:

```bash
cd /Users/ak87/ownCloud/owncloud/Calculations/TMP/AutoParty/autoparty
grep -n "keepH: true" app/home/templates/hitpicking_in_progress.html
grep -n "removeHs=False" app/base/io_utils.py app/base/tasks.py
```

### G) After changing defaults/config parsing code

Always restart both services:

1. Celery worker
2. Flask app

Reason:
- Worker/app processes keep old imported config/code until restart.

Path reliability note:
- Ensure `DEFAULT_FOLDER` in `autoparty/app/base/defaults.py` resolves to an absolute path (project `defaults/` directory), so Flask/Celery do not load different config files based on current working directory.

## 7) Stop/Cleanup

Stop Redis service:

```bash
brew services stop redis
```

Stop app/worker:
- Press `Ctrl+C` in the terminals running Flask and Celery.

## 8) Optional: Run Tests

From repository root:

```bash
conda run -n autoparty-env pip install pytest
mv tests/ autoparty/
mv testing_inputs/* autoparty/inputs
cd autoparty
conda run -n autoparty-env python -m pytest tests
```

Note: tests can take several minutes.
