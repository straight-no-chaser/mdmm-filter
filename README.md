# MD/MM evaluator for binder-receptor complexes

A library for batch evaluation of binder-receptor complexes. Preprocesses Rosetta-relaxed designs, runs explicit-solvent OpenMM MD, and calls AmberTools for MM/GBSA or MM/PBSA.

## API

Batch mode:

```python
from mdmm_filter import EvalConfig, run_pipeline

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
from mdmm_filter import EvalConfig, evaluate_one

cfg = EvalConfig(output_dir="results_single", do_md=False, do_mmpbsa=True, mmpbsa_mode="pbsa")
result = evaluate_one("designs/model_001.cif", cfg)
print(result)
```

`do_md=False` hydrogenates the solute with OpenMM and runs a static single-frame MMPBSA estimate. 

## CLI

Batch mode:

```bash
python -m mdmm_filter.cli run \
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
python -m mdmm_filter.cli one \
  --cif designs/model_001.cif \
  --output_dir results_single \
  --mmpbsa_mode pbsa
```