from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from .config import EvalConfig
from .io_utils import append_log, ensure_dir
from .preprocess import PreparedStructure, split_pdb_file_by_chains


@dataclass(slots=True)
class SoluteArtifacts:
    complex_pdb: Path
    receptor_pdb: Path
    binder_pdb: Path
    trajectory_path: Path
    n_frames: int


@dataclass(slots=True)
class MDResult:
    md_dir: Path
    solute: SoluteArtifacts
    final_complex_pdb: Path
    final_solvated_pdb: Path
    trajectory_dcd: Path
    final_potential_energy_kj_mol: float | None


def assert_openmm_available() -> None:
    _load_openmm()


def prepare_hydrogenated_structures(
    prepared: PreparedStructure,
    output_dir: Path,
    cfg: EvalConfig,
    log_path: Path,
) -> SoluteArtifacts:
    app, _, _, _ = _load_openmm()
    output_dir = ensure_dir(output_dir)
    append_log(log_path, f"Preparing hydrogenated solute from {prepared.complex_pdb}")
    forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
    pdb = app.PDBFile(str(prepared.complex_pdb))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    complex_pdb = output_dir / "complex_h.pdb"
    receptor_pdb = output_dir / "receptor_h.pdb"
    binder_pdb = output_dir / "binder_h.pdb"
    with complex_pdb.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    split_pdb_file_by_chains(complex_pdb, cfg.receptor_chain, cfg.binder_chain, receptor_pdb, binder_pdb)
    append_log(log_path, f"Wrote hydrogenated solute PDBs to {output_dir}")
    return SoluteArtifacts(
        complex_pdb=complex_pdb,
        receptor_pdb=receptor_pdb,
        binder_pdb=binder_pdb,
        trajectory_path=complex_pdb,
        n_frames=1,
    )


def run_explicit_solvent_md(
    prepared: PreparedStructure,
    output_dir: Path,
    cfg: EvalConfig,
    log_path: Path,
) -> MDResult:
    app, unit, LangevinMiddleIntegrator, Platform = _load_openmm()
    output_dir = ensure_dir(output_dir)
    append_log(log_path, f"Running explicit-solvent OpenMM MD in {output_dir}")

    forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
    pdb = app.PDBFile(str(prepared.complex_pdb))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)

    hydrogenated_complex = output_dir / "complex_h.pdb"
    hydrogenated_receptor = output_dir / "receptor_h.pdb"
    hydrogenated_binder = output_dir / "binder_h.pdb"
    with hydrogenated_complex.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    split_pdb_file_by_chains(
        hydrogenated_complex,
        cfg.receptor_chain,
        cfg.binder_chain,
        hydrogenated_receptor,
        hydrogenated_binder,
    )

    solute_topology = modeller.topology
    solute_atom_count = modeller.topology.getNumAtoms()

    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=cfg.water_padding_nm * unit.nanometer,
        ionicStrength=cfg.ionic_strength_m * unit.molar,
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
    )
    integrator = LangevinMiddleIntegrator(
        cfg.temperature_k * unit.kelvin,
        cfg.friction_per_ps / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )

    platform = Platform.getPlatformByName(cfg.openmm_platform) if cfg.openmm_platform else None
    if platform is None:
        simulation = app.Simulation(modeller.topology, system, integrator)
    else:
        simulation = app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=cfg.minimization_steps)
    simulation.context.setVelocitiesToTemperature(cfg.temperature_k * unit.kelvin)

    report_interval = min(_report_interval_steps(cfg), _ps_to_steps(cfg.md_time_ps, cfg.timestep_fs))
    report_interval = max(1, report_interval)
    equil_steps = _ps_to_steps(cfg.equil_time_ps, cfg.timestep_fs)
    production_steps = _ps_to_steps(cfg.md_time_ps, cfg.timestep_fs)
    if production_steps <= 0:
        raise RuntimeError("Production MD must run for at least one step")

    state_log = output_dir / "md.log"
    simulation.reporters.append(
        app.StateDataReporter(
            str(state_log),
            max(1, min(report_interval, 1000)),
            step=True,
            potentialEnergy=True,
            temperature=True,
            speed=True,
            volume=True,
            density=True,
        )
    )
    dcd_path = output_dir / "production_solute.dcd"
    solute_reporter = _SoluteDCDReporter(
        dcd_path=dcd_path,
        topology=solute_topology,
        atom_count=solute_atom_count,
        report_interval=report_interval,
        timestep_fs=cfg.timestep_fs,
        app_module=app,
    )

    if equil_steps:
        append_log(log_path, f"Equilibrating for {cfg.equil_time_ps} ps ({equil_steps} steps)")
        simulation.step(equil_steps)
    simulation.reporters.append(solute_reporter)
    append_log(log_path, f"Producing for {cfg.md_time_ps} ps ({production_steps} steps)")
    simulation.step(production_steps)

    final_state = simulation.context.getState(getPositions=True, getEnergy=True)
    final_solvated_pdb = output_dir / "final_solvated.pdb"
    final_complex_pdb = output_dir / "final_complex.pdb"
    with final_solvated_pdb.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(modeller.topology, final_state.getPositions(), handle)
    with final_complex_pdb.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(solute_topology, final_state.getPositions()[:solute_atom_count], handle)
    solute_reporter.close()

    final_potential = final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    append_log(log_path, f"Final potential energy: {final_potential:.3f} kJ/mol")
    return MDResult(
        md_dir=output_dir,
        solute=SoluteArtifacts(
            complex_pdb=hydrogenated_complex,
            receptor_pdb=hydrogenated_receptor,
            binder_pdb=hydrogenated_binder,
            trajectory_path=dcd_path,
            n_frames=solute_reporter.frames_written,
        ),
        final_complex_pdb=final_complex_pdb,
        final_solvated_pdb=final_solvated_pdb,
        trajectory_dcd=dcd_path,
        final_potential_energy_kj_mol=final_potential,
    )


class _SoluteDCDReporter:
    def __init__(
        self,
        *,
        dcd_path: Path,
        topology,
        atom_count: int,
        report_interval: int,
        timestep_fs: float,
        app_module,
    ) -> None:
        self._dcd_path = dcd_path
        self._topology = topology
        self._atom_count = atom_count
        self._report_interval = report_interval
        self._timestep_fs = timestep_fs
        self._app = app_module
        self._handle = None
        self._writer = None
        self.frames_written = 0

    def describeNextReport(self, simulation):
        steps = self._report_interval - simulation.currentStep % self._report_interval
        return steps, True, False, False, False, None

    def report(self, simulation, state) -> None:
        if self._writer is None:
            self._handle = self._dcd_path.open("wb")
            dt_ps = (self._report_interval * self._timestep_fs) / 1000.0
            self._writer = self._app.DCDFile(self._handle, self._topology, dt_ps)
        positions = state.getPositions(asNumpy=True)[: self._atom_count]
        self._writer.writeModel(positions, periodicBoxVectors=state.getPeriodicBoxVectors())
        self.frames_written += 1

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None


def _ps_to_steps(time_ps: float, timestep_fs: float) -> int:
    return int(round(time_ps * 1000.0 / timestep_fs))


def _report_interval_steps(cfg: EvalConfig) -> int:
    target = max(1, int(round(cfg.frame_interval_ps * 1000.0 / cfg.timestep_fs)))
    total = _ps_to_steps(cfg.md_time_ps, cfg.timestep_fs)
    cap = max(1, math.ceil(total / cfg.max_frames))
    return max(target, cap)


def _load_openmm():
    try:
        from openmm import LangevinMiddleIntegrator, Platform, unit
        from openmm import app
    except ImportError as exc:
        raise RuntimeError(
            "OpenMM is required for hydrogenation and MD. Install the Python 'openmm' package before running this workflow."
        ) from exc
    return app, unit, LangevinMiddleIntegrator, Platform
