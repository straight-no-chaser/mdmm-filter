# MD/MM batch evaluation of binder-receptor complexes

This repository contains `mdmm_filter_cuilab`, which is a minimal Python library for batch evaluation of binder-receptor complexes. It preprocesses Rosetta-relaxed designs, runs explicit-solvent OpenMM MD, launches AmberTools `MMPBSA.py` for MM/GBSA or MM/PBSA, and writes ranked summaries.

This repository is a part of ongoing research efforts in ADC design and molecular modeling in the [Cui Lab](https://cuilab.stanford.edu/research), Stanford University.

## Scope

- Input format: one `.cif` or `.mmcif` file per design.
- Chain convention: receptor is chain `A`, binder is chain `B` by default.
- Intended use case: protein-only binder design screening from RFD, MPNN, or AF3 outputs after relax.
- Intermediate format: preprocessing writes protein-only PDB files for OpenMM and AmberTools.

Explicit-solvent dynamics are delegated to OpenMM, and MM/PBSA or MM/GBSA is delegated to AmberTools.

## Dependencies

Python dependency:

- `openmm`

System dependencies expected on `PATH`:

- `tleap`
- `MMPBSA.py`

In practice those two executables come from AmberTools.

### Responsibility split

- OpenMM: add hydrogens, build the explicit-solvent system, minimize, equilibrate, run production MD, and write trajectories.
- AmberTools: build AMBER topology files with `tleap` and run MM/GBSA or MM/PBSA with `MMPBSA.py`.

## Installation

Install the library after OpenMM is available in your Python environment:

```bash
pip install -e .
```

AmberTools must already be installed separately so that `tleap` and `MMPBSA.py` are callable from the shell.

## Python API

```python
from mdmm_filter_cuilab import EvalConfig, run_pipeline

cfg = EvalConfig(
    input_dir="designs",
    output_dir="results",
    receptor_chain="A",
    binder_chain="B",
    do_md=True,
    do_mmpbsa=True,
    mmpbsa_mode="gbsa",
    md_time_ps=2000,
    temperature_k=300,
    ionic_strength_m=0.15,
    keep_intermediate=False,
    n_workers=4,
)

results = run_pipeline(cfg)
```

Single-design mode:

```python
from mdmm_filter_cuilab import EvalConfig, evaluate_one

cfg = EvalConfig(output_dir="results_single", do_md=False, do_mmpbsa=True, mmpbsa_mode="pbsa")
result = evaluate_one("designs/model_001.cif", cfg)
print(result)
```

If `do_md=False`, the package still hydrogenates the solute with OpenMM and runs a static single-frame MMPBSA estimate. That is much rougher than running MD and evaluating multiple extracted frames.

## CLI

Batch mode:

```bash
python -m mdmm_filter_cuilab.cli run \
  --input_dir designs \
  --output_dir results \
  --receptor_chain A \
  --binder_chain B \
  --do_md \
  --do_mmpbsa \
  --mmpbsa_mode gbsa \
  --md_time_ps 2000 \
  --n_workers 4
```

Single-design mode:

```bash
python -m mdmm_filter_cuilab.cli one \
  --cif designs/model_001.cif \
  --output_dir results_single \
  --mmpbsa_mode pbsa
```

## Output layout

For a design named `model_001`:

```text
results/
  model_001/
    design.log
    preprocess/
      complex.pdb
      receptor.pdb
      binder.pdb
    md/
      complex_h.pdb
      receptor_h.pdb
      binder_h.pdb
      production_solute.dcd
      final_complex.pdb
      final_solvated.pdb
      md.log
    mmpbsa/
      tleap.in
      complex.prmtop
      receptor.prmtop
      binder.prmtop
      mmpbsa_gbsa.in
      FINAL_RESULTS_MMPBSA.dat
      FINAL_ENERGY.csv
  summary.csv
  summary.json
```

`summary.csv` and `summary.json` contain at least:

- `design_name`
- `input_cif`
- `status`
- `md_completed`
- `mmpbsa_completed`
- `mmpbsa_mode`
- `deltaG_total`
- `n_frames`
- `md_final_potential_energy`
- `notes`

Ranking is ascending on `deltaG_total`, so more negative values rank better.

## Input requirements

- Each mmCIF file must contain both chains `A` and `B`, unless you override the chain IDs in `EvalConfig`.
- The workflow expects protein residues that Amber protein force fields can parameterize.
- Residues outside the standard protein set are rejected during preprocessing.
- Duplicate atom names within the same residue are rejected.

## Notes on MD defaults

The default workflow is intentionally short for screening:

- energy minimization: `5000` iterations
- equilibration: `100 ps`
- production: `2000 ps`
- frame target: every `10 ps`, capped at `100` frames

Those settings are useful for triage, not for publication-grade binding free energies.

## Debugging

- OpenMM is not installed in the active Python environment.
- `tleap` or `MMPBSA.py` is missing from `PATH`.
- The input CIF is missing chain `A` or `B`.
- A chain is present but empty after filtering.
- A residue is non-standard and cannot be parameterized by the selected protein force field.
- `MMPBSA.py` fails because AmberTools cannot interpret the generated topology or trajectory.
- Each design writes its own `design.log`, so a single failure does not abort the full batch.
