from __future__ import annotations

import re
import shutil
from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path

from .config import EvalConfig
from .io_utils import append_log, ensure_dir, run_command, write_text


@dataclass(slots=True)
class MMPBSAResult:
    mode: str
    completed: bool
    delta_g_total: float | None
    n_frames: int
    energy_terms: dict[str, float] = field(default_factory=dict)
    output_file: Path | None = None


def assert_ambertools_available(mmpbsa_executable: str) -> None:
    if shutil.which("tleap") is None:
        raise RuntimeError("AmberTools 'tleap' was not found in PATH")
    if shutil.which(mmpbsa_executable) is None:
        raise RuntimeError(f"AmberTools '{mmpbsa_executable}' was not found in PATH")


def run_mmpbsa(
    *,
    complex_pdb: Path,
    receptor_pdb: Path,
    binder_pdb: Path,
    trajectory_path: Path,
    n_frames: int,
    output_dir: Path,
    cfg: EvalConfig,
    log_path: Path,
) -> MMPBSAResult:
    output_dir = ensure_dir(output_dir)
    assert_ambertools_available(cfg.mmpbsa_executable)

    leap_input = output_dir / "tleap.in"
    mmpbsa_input = output_dir / f"mmpbsa_{cfg.mmpbsa_mode}.in"
    results_file = output_dir / "FINAL_RESULTS_MMPBSA.dat"
    energy_csv = output_dir / "FINAL_ENERGY.csv"

    _write_leap_input(leap_input, complex_pdb, receptor_pdb, binder_pdb)
    leap_result = run_command(["tleap", "-f", str(leap_input)], cwd=output_dir, log_path=log_path, description="tleap")
    if leap_result.returncode != 0:
        raise RuntimeError("tleap failed; see design log for details")

    _write_mmpbsa_input(mmpbsa_input, cfg)
    command = [
        cfg.mmpbsa_executable,
        "-O",
        "-i",
        str(mmpbsa_input),
        "-o",
        str(results_file),
        "-eo",
        str(energy_csv),
        "-cp",
        "complex.prmtop",
        "-rp",
        "receptor.prmtop",
        "-lp",
        "binder.prmtop",
        "-y",
        str(trajectory_path),
    ]
    mmpbsa_result = run_command(command, cwd=output_dir, log_path=log_path, description="MMPBSA.py")
    if mmpbsa_result.returncode != 0:
        raise RuntimeError("MMPBSA.py failed; see design log for details")
    if not results_file.exists():
        raise RuntimeError("MMPBSA.py did not produce FINAL_RESULTS_MMPBSA.dat")

    energy_terms = parse_mmpbsa_output(results_file)
    delta_total = energy_terms.get("delta_total")
    if delta_total is None:
        raise RuntimeError("Could not parse DELTA TOTAL from FINAL_RESULTS_MMPBSA.dat")
    append_log(log_path, f"MMPBSA {cfg.mmpbsa_mode} DELTA TOTAL = {delta_total:.4f} kcal/mol")
    return MMPBSAResult(
        mode=cfg.mmpbsa_mode,
        completed=True,
        delta_g_total=delta_total,
        n_frames=n_frames,
        energy_terms=energy_terms,
        output_file=results_file,
    )


def parse_mmpbsa_output(results_file: Path) -> dict[str, float]:
    text = results_file.read_text(encoding="utf-8", errors="replace")
    terms: dict[str, float] = {}
    in_diff_block = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("Differences (Complex - Receptor - Ligand):"):
            in_diff_block = True
            continue
        if not in_diff_block:
            continue
        if not stripped or stripped.startswith("-") or stripped.startswith("Energy Component"):
            continue
        match = re.match(r"^([A-Z][A-Z0-9_ ]+?)\s+(-?\d+(?:\.\d+)?)", stripped)
        if not match:
            continue
        key = match.group(1).strip().lower().replace(" ", "_")
        terms[key] = float(match.group(2))
    return terms


def _write_leap_input(path: Path, complex_pdb: Path, receptor_pdb: Path, binder_pdb: Path) -> None:
    text = "\n".join(
        [
            "source leaprc.protein.ff14SB",
            "set default PBRadii mbondi3",
            f'complex = loadpdb "{complex_pdb.as_posix()}"',
            'saveamberparm complex "complex.prmtop" "complex.inpcrd"',
            f'receptor = loadpdb "{receptor_pdb.as_posix()}"',
            'saveamberparm receptor "receptor.prmtop" "receptor.inpcrd"',
            f'binder = loadpdb "{binder_pdb.as_posix()}"',
            'saveamberparm binder "binder.prmtop" "binder.inpcrd"',
            "quit",
            "",
        ]
    )
    write_text(path, text)


def _write_mmpbsa_input(path: Path, cfg: EvalConfig) -> None:
    template_name = "mmpbsa_gbsa.in" if cfg.mmpbsa_mode == "gbsa" else "mmpbsa_pbsa.in"
    template = (files("mdmm_filter_cuilab") / "templates" / template_name).read_text(encoding="utf-8")
    rendered = template.format(
        startframe=1,
        interval=1,
        keep_files=1 if cfg.keep_intermediate else 0,
        saltcon=f"{cfg.ionic_strength_m:.3f}",
        istrng=f"{cfg.ionic_strength_m:.3f}",
    )
    write_text(path, rendered)
