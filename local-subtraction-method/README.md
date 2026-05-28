# lsmethod

Research code for symbolic and numerical checks of local subtraction method kernels.

## Setup

Create the environment:

```bash
mamba env create -f environment.yml
mamba activate lsmethod
python -m pip install -e .
```

Run tests:

```bash
pytest
```

Run the first numerical checks:

```bash
python scripts/check_closed_vs_w.py
python scripts/check_closed_vs_u.py
```

## Structure

- `src/lsmethod/kinematics.py`: variables and derived kinematic quantities.
- `src/lsmethod/kernels.py`: closed form and parameter-integral representations.
- `src/lsmethod/equations.py`: map from code names to LaTeX labels.
- `tests/`: numerical and structural checks.
- `scripts/`: commands for manual checks and scans.
