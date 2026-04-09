from __future__ import annotations

from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Iterable

from .io_utils import ensure_dir, write_csv, write_json

SUMMARY_FIELDS = [
    "design_name",
    "input_cif",
    "status",
    "md_completed",
    "mmpbsa_completed",
    "mmpbsa_mode",
    "deltaG_total",
    "n_frames",
    "md_final_potential_energy",
    "notes",
]


def write_summary(results: Iterable[Any], output_dir: Path) -> tuple[Path, Path]:
    output_dir = ensure_dir(output_dir)
    rows = sorted((_result_to_row(result) for result in results), key=_sort_key)
    csv_path = output_dir / "summary.csv"
    json_path = output_dir / "summary.json"
    write_csv(csv_path, rows, SUMMARY_FIELDS)
    write_json(json_path, rows)
    return csv_path, json_path


def _result_to_row(result: Any) -> dict[str, Any]:
    if is_dataclass(result):
        data = asdict(result)
    elif isinstance(result, dict):
        data = dict(result)
    else:
        data = {name: getattr(result, name) for name in SUMMARY_FIELDS if hasattr(result, name)}
    return {
        "design_name": data.get("design_name", ""),
        "input_cif": data.get("input_cif", ""),
        "status": data.get("status", ""),
        "md_completed": data.get("md_completed", False),
        "mmpbsa_completed": data.get("mmpbsa_completed", False),
        "mmpbsa_mode": data.get("mmpbsa_mode", ""),
        "deltaG_total": data.get("deltaG_total"),
        "n_frames": data.get("n_frames", 0),
        "md_final_potential_energy": data.get("md_final_potential_energy"),
        "notes": data.get("notes", ""),
    }


def _sort_key(row: dict[str, Any]) -> tuple[int, float, str]:
    delta = row.get("deltaG_total")
    if delta is None:
        return (1, float("inf"), row.get("design_name", ""))
    return (0, float(delta), row.get("design_name", ""))
