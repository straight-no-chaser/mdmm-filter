from __future__ import annotations

import csv
import json
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def timestamp() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def append_log(log_path: Path, message: str, *, echo: bool = False) -> None:
    ensure_dir(log_path.parent)
    line = f"[{timestamp()}] {message.rstrip()}\n"
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(line)
    if echo:
        print(line, end="")


def list_cif_files(input_dir: Path) -> list[Path]:
    files = sorted(input_dir.glob("*.cif"))
    files.extend(sorted(input_dir.glob("*.mmcif")))
    deduped: dict[Path, None] = {}
    for path in files:
        deduped[path.resolve()] = None
    return list(deduped.keys())


def run_command(
    command: list[str],
    *,
    cwd: Path,
    log_path: Path,
    description: str,
) -> subprocess.CompletedProcess[str]:
    append_log(log_path, f"{description}: {shlex.join(command)}")
    result = subprocess.run(
        command,
        cwd=str(cwd),
        check=False,
        capture_output=True,
        text=True,
    )
    if result.stdout:
        append_log(log_path, f"{description} stdout:\n{result.stdout}")
    if result.stderr:
        append_log(log_path, f"{description} stderr:\n{result.stderr}")
    append_log(log_path, f"{description} return code: {result.returncode}")
    return result


def write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


def write_json(path: Path, payload: Any) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
