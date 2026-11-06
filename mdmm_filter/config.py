from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class EvalConfig:
    input_dir: str = ""
    output_dir: str = "results"
    receptor_chain: str = "A"
    binder_chain: str = "B"
    do_md: bool = True
    do_mmpbsa: bool = True
    mmpbsa_mode: str = "gbsa"
    md_time_ps: int = 2000
    equil_time_ps: int = 100
    minimization_steps: int = 5000
    temperature_k: float = 300.0
    friction_per_ps: float = 1.0
    timestep_fs: float = 2.0
    ionic_strength_m: float = 0.15
    water_padding_nm: float = 1.0
    frame_interval_ps: float = 10.0
    max_frames: int = 100
    keep_intermediate: bool = False
    n_workers: int = 1
    openmm_platform: str | None = None
    mmpbsa_executable: str = "MMPBSA.py"
    verbose: bool = True

    def validate(self) -> None:
        if not self.output_dir:
            raise ValueError("output_dir must not be empty")
        if not self.receptor_chain:
            raise ValueError("receptor_chain must not be empty")
        if not self.binder_chain:
            raise ValueError("binder_chain must not be empty")
        if self.receptor_chain == self.binder_chain:
            raise ValueError("receptor_chain and binder_chain must differ")
        if not (self.do_md or self.do_mmpbsa):
            raise ValueError("At least one of do_md or do_mmpbsa must be enabled")
        if self.mmpbsa_mode.lower() not in {"gbsa", "pbsa"}:
            raise ValueError("mmpbsa_mode must be 'gbsa' or 'pbsa'")
        if self.md_time_ps <= 0:
            raise ValueError("md_time_ps must be positive")
        if self.equil_time_ps < 0:
            raise ValueError("equil_time_ps must be non-negative")
        if self.minimization_steps <= 0:
            raise ValueError("minimization_steps must be positive")
        if self.temperature_k <= 0:
            raise ValueError("temperature_k must be positive")
        if self.friction_per_ps <= 0:
            raise ValueError("friction_per_ps must be positive")
        if self.timestep_fs <= 0:
            raise ValueError("timestep_fs must be positive")
        if self.ionic_strength_m < 0:
            raise ValueError("ionic_strength_m must be non-negative")
        if self.water_padding_nm <= 0:
            raise ValueError("water_padding_nm must be positive")
        if self.frame_interval_ps <= 0:
            raise ValueError("frame_interval_ps must be positive")
        if self.max_frames <= 0:
            raise ValueError("max_frames must be positive")
        if self.n_workers <= 0:
            raise ValueError("n_workers must be positive")
        if not self.mmpbsa_executable:
            raise ValueError("mmpbsa_executable must not be empty")
        self.mmpbsa_mode = self.mmpbsa_mode.lower()

    def validate_batch(self) -> None:
        self.validate()
        if not self.input_dir:
            raise ValueError("input_dir must not be empty for batch runs")
        if not Path(self.input_dir).exists():
            raise ValueError(f"input_dir does not exist: {self.input_dir}")

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "EvalConfig":
        return cls(**data)
