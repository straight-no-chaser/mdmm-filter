from __future__ import annotations

import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

from .config import EvalConfig
from .io_utils import append_log, ensure_dir, list_cif_files
from .md import MDResult, SoluteArtifacts, assert_openmm_available, prepare_hydrogenated_structures, run_explicit_solvent_md
from .mmpbsa import MMPBSAResult, assert_ambertools_available, run_mmpbsa
from .preprocess import preprocess_structure
from .ranking import write_summary


@dataclass(slots=True)
class DesignResult:
    design_name: str
    input_cif: str
    status: str
    md_completed: bool
    mmpbsa_completed: bool
    mmpbsa_mode: str
    deltaG_total: float | None = None
    n_frames: int = 0
    md_final_potential_energy: float | None = None
    notes: str = ""
    output_dir: str = ""
    energy_terms: dict[str, float] = field(default_factory=dict)


def run_pipeline(cfg: EvalConfig) -> list[DesignResult]:
    cfg.validate_batch()
    _check_runtime_requirements(cfg)
    output_dir = ensure_dir(Path(cfg.output_dir))
    cif_paths = list_cif_files(Path(cfg.input_dir))
    if not cif_paths:
        raise ValueError(f"No .cif or .mmcif files found in {cfg.input_dir}")

    results: list[DesignResult] = []
    if cfg.n_workers == 1 or len(cif_paths) == 1:
        for cif_path in cif_paths:
            results.append(evaluate_one(cif_path, cfg))
    else:
        with ProcessPoolExecutor(max_workers=cfg.n_workers) as executor:
            futures = {
                executor.submit(_evaluate_one_worker, str(cif_path), cfg.as_dict()): cif_path
                for cif_path in cif_paths
            }
            for future in as_completed(futures):
                results.append(future.result())

    write_summary(results, output_dir)
    return sorted(
        results,
        key=lambda item: (
            item.deltaG_total is None,
            float("inf") if item.deltaG_total is None else item.deltaG_total,
            item.design_name,
        ),
    )


def evaluate_one(cif_path: str | Path, cfg: EvalConfig) -> DesignResult:
    cfg.validate()
    _check_runtime_requirements(cfg)
    cif_path = Path(cif_path).resolve()
    design_name = cif_path.stem
    design_dir = ensure_dir(Path(cfg.output_dir) / design_name)
    log_path = design_dir / "design.log"
    result = DesignResult(
        design_name=design_name,
        input_cif=str(cif_path),
        status="failed",
        md_completed=False,
        mmpbsa_completed=False,
        mmpbsa_mode=cfg.mmpbsa_mode,
        output_dir=str(design_dir.resolve()),
    )
    try:
        append_log(log_path, f"Starting evaluation for {design_name}", echo=cfg.verbose)
        prepared = preprocess_structure(
            cif_path=cif_path,
            output_root=design_dir,
            receptor_chain=cfg.receptor_chain,
            binder_chain=cfg.binder_chain,
            log_path=log_path,
        )

        solute_artifacts: SoluteArtifacts | None = None
        md_result: MDResult | None = None
        if cfg.do_md:
            md_result = run_explicit_solvent_md(prepared, design_dir / "md", cfg, log_path)
            solute_artifacts = md_result.solute
            result.md_completed = True
            result.n_frames = md_result.solute.n_frames
            result.md_final_potential_energy = md_result.final_potential_energy_kj_mol
        elif cfg.do_mmpbsa:
            solute_artifacts = prepare_hydrogenated_structures(prepared, design_dir / "static", cfg, log_path)
            result.n_frames = solute_artifacts.n_frames

        if cfg.do_mmpbsa and solute_artifacts is not None:
            mmpbsa_result = run_mmpbsa(
                complex_pdb=solute_artifacts.complex_pdb,
                receptor_pdb=solute_artifacts.receptor_pdb,
                binder_pdb=solute_artifacts.binder_pdb,
                trajectory_path=solute_artifacts.trajectory_path,
                n_frames=solute_artifacts.n_frames,
                output_dir=design_dir / "mmpbsa",
                cfg=cfg,
                log_path=log_path,
            )
            _apply_mmpbsa_result(result, mmpbsa_result)

        result.status = _final_status(result, cfg)
        append_log(log_path, f"Finished {design_name} with status={result.status}", echo=cfg.verbose)
    except Exception as exc:
        result.notes = str(exc)
        result.status = "failed"
        append_log(log_path, traceback.format_exc(), echo=cfg.verbose)
    return result


def _apply_mmpbsa_result(result: DesignResult, mmpbsa_result: MMPBSAResult) -> None:
    result.mmpbsa_completed = mmpbsa_result.completed
    result.deltaG_total = mmpbsa_result.delta_g_total
    result.n_frames = mmpbsa_result.n_frames
    result.energy_terms = dict(mmpbsa_result.energy_terms)


def _final_status(result: DesignResult, cfg: EvalConfig) -> str:
    requested = []
    if cfg.do_md:
        requested.append(result.md_completed)
    if cfg.do_mmpbsa:
        requested.append(result.mmpbsa_completed)
    if requested and all(requested):
        return "completed"
    if any(requested):
        return "partial"
    return "completed"


def _evaluate_one_worker(cif_path: str, cfg_dict: dict) -> DesignResult:
    cfg = EvalConfig.from_dict(cfg_dict)
    return evaluate_one(cif_path, cfg)


def _check_runtime_requirements(cfg: EvalConfig) -> None:
    if cfg.do_md or cfg.do_mmpbsa:
        assert_openmm_available()
    if cfg.do_mmpbsa:
        assert_ambertools_available(cfg.mmpbsa_executable)
