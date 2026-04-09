from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from .config import EvalConfig
from .pipeline import evaluate_one, run_pipeline


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    cfg = config_from_args(args)
    if args.command == "run":
        results = run_pipeline(cfg)
        print(json.dumps([asdict(result) for result in results], indent=2))
    elif args.command == "one":
        result = evaluate_one(args.cif, cfg)
        print(json.dumps(asdict(result), indent=2))
    else:
        parser.error("A subcommand is required")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mdmm_filter_cuilab", description="Batch binder evaluation workflow")
    subparsers = parser.add_subparsers(dest="command")

    run_parser = subparsers.add_parser("run", help="Evaluate all CIF files in a directory")
    run_parser.add_argument("--input_dir", required=True)
    add_shared_arguments(run_parser)

    one_parser = subparsers.add_parser("one", help="Evaluate a single CIF file")
    one_parser.add_argument("--cif", required=True)
    add_shared_arguments(one_parser)

    return parser


def add_shared_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--receptor_chain", default="A")
    parser.add_argument("--binder_chain", default="B")
    parser.add_argument("--do_md", action="store_true")
    parser.add_argument("--no_md", action="store_true")
    parser.add_argument("--do_mmpbsa", action="store_true")
    parser.add_argument("--no_mmpbsa", action="store_true")
    parser.add_argument("--mmpbsa_mode", choices=["gbsa", "pbsa"], default="gbsa")
    parser.add_argument("--md_time_ps", type=int, default=2000)
    parser.add_argument("--equil_time_ps", type=int, default=100)
    parser.add_argument("--minimization_steps", type=int, default=5000)
    parser.add_argument("--temperature_k", type=float, default=300.0)
    parser.add_argument("--friction_per_ps", type=float, default=1.0)
    parser.add_argument("--timestep_fs", type=float, default=2.0)
    parser.add_argument("--ionic_strength_m", type=float, default=0.15)
    parser.add_argument("--water_padding_nm", type=float, default=1.0)
    parser.add_argument("--frame_interval_ps", type=float, default=10.0)
    parser.add_argument("--max_frames", type=int, default=100)
    parser.add_argument("--keep_intermediate", action="store_true")
    parser.add_argument("--n_workers", type=int, default=1)
    parser.add_argument("--openmm_platform")
    parser.add_argument("--mmpbsa_executable", default="MMPBSA.py")
    parser.add_argument("--quiet", action="store_true")


def config_from_args(args: argparse.Namespace) -> EvalConfig:
    do_md = not args.no_md if not args.do_md else True
    do_mmpbsa = not args.no_mmpbsa if not args.do_mmpbsa else True
    input_dir = args.input_dir if hasattr(args, "input_dir") else str(Path(args.cif).resolve().parent)
    return EvalConfig(
        input_dir=input_dir,
        output_dir=args.output_dir,
        receptor_chain=args.receptor_chain,
        binder_chain=args.binder_chain,
        do_md=do_md,
        do_mmpbsa=do_mmpbsa,
        mmpbsa_mode=args.mmpbsa_mode,
        md_time_ps=args.md_time_ps,
        equil_time_ps=args.equil_time_ps,
        minimization_steps=args.minimization_steps,
        temperature_k=args.temperature_k,
        friction_per_ps=args.friction_per_ps,
        timestep_fs=args.timestep_fs,
        ionic_strength_m=args.ionic_strength_m,
        water_padding_nm=args.water_padding_nm,
        frame_interval_ps=args.frame_interval_ps,
        max_frames=args.max_frames,
        keep_intermediate=args.keep_intermediate,
        n_workers=args.n_workers,
        openmm_platform=args.openmm_platform,
        mmpbsa_executable=args.mmpbsa_executable,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
