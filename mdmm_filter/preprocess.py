from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .io_utils import append_log, ensure_dir

STANDARD_PROTEIN_RESIDUES = {
    "ACE",
    "ALA",
    "ARG",
    "ASH",
    "ASN",
    "ASP",
    "CYS",
    "CYX",
    "GLH",
    "GLN",
    "GLU",
    "GLY",
    "HID",
    "HIE",
    "HIP",
    "HIS",
    "ILE",
    "LEU",
    "LYN",
    "LYS",
    "MET",
    "NME",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


class PreprocessError(RuntimeError):
    pass


@dataclass(slots=True)
class CifAtom:
    record_type: str
    atom_name: str
    residue_name: str
    chain_id: str
    residue_id: str
    insertion_code: str
    x: float
    y: float
    z: float
    element: str
    alt_id: str = ""
    model_number: str = "1"


@dataclass(slots=True)
class PreparedStructure:
    design_name: str
    source_cif: Path
    design_dir: Path
    preprocess_dir: Path
    complex_pdb: Path
    receptor_pdb: Path
    binder_pdb: Path
    atom_count: int
    receptor_atom_count: int
    binder_atom_count: int


def parse_mmcif_atoms(cif_path: Path) -> list[CifAtom]:
    lines = cif_path.read_text(encoding="utf-8").splitlines()
    headers: list[str] = []
    rows: list[list[str]] = []
    index = 0
    while index < len(lines):
        stripped = lines[index].strip()
        if stripped != "loop_":
            index += 1
            continue
        probe = index + 1
        loop_headers: list[str] = []
        while probe < len(lines) and lines[probe].lstrip().startswith("_"):
            loop_headers.append(lines[probe].strip())
            probe += 1
        if loop_headers and loop_headers[0].startswith("_atom_site."):
            headers = loop_headers
            buffer: list[str] = []
            while probe < len(lines):
                candidate = lines[probe].strip()
                if not candidate or candidate.startswith("#"):
                    probe += 1
                    if rows:
                        break
                    continue
                if candidate == "loop_" or candidate.startswith("data_") or candidate.startswith("_"):
                    break
                buffer.extend(shlex.split(candidate, posix=True))
                while len(buffer) >= len(headers):
                    rows.append(buffer[: len(headers)])
                    buffer = buffer[len(headers) :]
                probe += 1
            break
        index = probe
    if not headers or not rows:
        raise PreprocessError(f"No _atom_site loop found in {cif_path}")

    header_index = {name.removeprefix("_atom_site."): idx for idx, name in enumerate(headers)}

    def value(row: list[str], *keys: str) -> str:
        for key in keys:
            idx = header_index.get(key)
            if idx is not None and idx < len(row):
                token = row[idx]
                if token not in {".", "?"}:
                    return token
        return ""

    atoms: list[CifAtom] = []
    seen_atom_keys: set[tuple[str, str, str, str, str]] = set()
    for row in rows:
        record_type = value(row, "group_PDB").upper() or "ATOM"
        if record_type not in {"ATOM", "HETATM"}:
            continue
        model_number = value(row, "pdbx_PDB_model_num") or "1"
        if model_number != "1":
            continue
        alt_id = value(row, "label_alt_id")
        if alt_id and alt_id not in {"A", "1"}:
            continue
        atom_name = value(row, "auth_atom_id", "label_atom_id")
        atom = CifAtom(
            record_type=record_type,
            atom_name=atom_name,
            residue_name=value(row, "auth_comp_id", "label_comp_id"),
            chain_id=value(row, "auth_asym_id", "label_asym_id"),
            residue_id=value(row, "auth_seq_id", "label_seq_id"),
            insertion_code=value(row, "pdbx_PDB_ins_code"),
            x=float(value(row, "Cartn_x")),
            y=float(value(row, "Cartn_y")),
            z=float(value(row, "Cartn_z")),
            element=(value(row, "type_symbol") or atom_name_to_element(atom_name)).upper(),
            alt_id=alt_id,
            model_number=model_number,
        )
        if not atom.atom_name or not atom.residue_name or not atom.chain_id or not atom.residue_id:
            raise PreprocessError(f"Incomplete atom_site record in {cif_path}")
        atom_key = (
            atom.chain_id,
            atom.residue_id,
            atom.insertion_code,
            atom.residue_name,
            atom.atom_name,
        )
        if atom_key in seen_atom_keys:
            raise PreprocessError(
                f"Duplicate atom name conflict for {atom.chain_id}:{atom.residue_name}{atom.residue_id} {atom.atom_name}"
            )
        seen_atom_keys.add(atom_key)
        atoms.append(atom)
    if not atoms:
        raise PreprocessError(f"No usable atoms found in {cif_path}")
    return atoms


def atom_name_to_element(atom_name: str) -> str:
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if not letters:
        return "X"
    if len(letters) >= 2 and letters[1].islower():
        return letters[:2]
    return letters[0]


def split_chain_atoms(
    atoms: Iterable[CifAtom],
    receptor_chain: str,
    binder_chain: str,
) -> tuple[list[CifAtom], list[CifAtom], list[CifAtom]]:
    selected = [atom for atom in atoms if atom.chain_id in {receptor_chain, binder_chain}]
    receptor_atoms = [atom for atom in selected if atom.chain_id == receptor_chain]
    binder_atoms = [atom for atom in selected if atom.chain_id == binder_chain]
    if not receptor_atoms:
        raise PreprocessError(f"Missing or empty receptor chain: {receptor_chain}")
    if not binder_atoms:
        raise PreprocessError(f"Missing or empty binder chain: {binder_chain}")
    _check_standard_residues(receptor_atoms + binder_atoms)
    return selected, receptor_atoms, binder_atoms


def _check_standard_residues(atoms: Iterable[CifAtom]) -> None:
    seen_residues: set[tuple[str, str, str, str]] = set()
    for atom in atoms:
        residue_key = (atom.chain_id, atom.residue_name, atom.residue_id, atom.insertion_code)
        if residue_key in seen_residues:
            continue
        seen_residues.add(residue_key)
        if atom.residue_name.upper() not in STANDARD_PROTEIN_RESIDUES:
            raise PreprocessError(
                f"Non-standard residue cannot be parameterized: {atom.chain_id}:{atom.residue_name}{atom.residue_id}"
            )


def preprocess_structure(
    cif_path: Path,
    output_root: Path,
    receptor_chain: str,
    binder_chain: str,
    log_path: Path,
) -> PreparedStructure:
    design_name = cif_path.stem
    preprocess_dir = ensure_dir(output_root / "preprocess")
    append_log(log_path, f"Preprocessing {cif_path}")
    atoms = parse_mmcif_atoms(cif_path)
    complex_atoms, receptor_atoms, binder_atoms = split_chain_atoms(atoms, receptor_chain, binder_chain)
    complex_pdb = preprocess_dir / "complex.pdb"
    receptor_pdb = preprocess_dir / "receptor.pdb"
    binder_pdb = preprocess_dir / "binder.pdb"
    write_pdb(complex_atoms, complex_pdb)
    write_pdb(receptor_atoms, receptor_pdb)
    write_pdb(binder_atoms, binder_pdb)
    append_log(log_path, f"Wrote preprocessed PDBs to {preprocess_dir}")
    return PreparedStructure(
        design_name=design_name,
        source_cif=cif_path.resolve(),
        design_dir=output_root,
        preprocess_dir=preprocess_dir,
        complex_pdb=complex_pdb,
        receptor_pdb=receptor_pdb,
        binder_pdb=binder_pdb,
        atom_count=len(complex_atoms),
        receptor_atom_count=len(receptor_atoms),
        binder_atom_count=len(binder_atoms),
    )


def write_pdb(atoms: Iterable[CifAtom], output_path: Path) -> None:
    ensure_dir(output_path.parent)
    atoms_list = list(atoms)
    lines: list[str] = []
    previous_chain: str | None = None
    serial = 1
    for atom in atoms_list:
        if previous_chain is not None and atom.chain_id != previous_chain:
            lines.append("TER")
        residue_id = atom.residue_id[-4:].rjust(4)
        insertion_code = (atom.insertion_code or " ")[:1]
        chain_id = atom.chain_id[:1]
        atom_name = format_atom_name(atom.atom_name)
        element = atom.element[:2].rjust(2)
        lines.append(
            f"{atom.record_type:<6}{serial:>5} {atom_name}{atom.residue_name[:3]:>4} {chain_id}"
            f"{residue_id}{insertion_code}   {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}"
            f"{1.00:>6.2f}{0.00:>6.2f}          {element}"
        )
        serial += 1
        previous_chain = atom.chain_id
    if atoms_list:
        lines.append("TER")
    lines.append("END")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def format_atom_name(atom_name: str) -> str:
    trimmed = atom_name[:4]
    if len(trimmed) == 4:
        return trimmed
    if trimmed and trimmed[0].isdigit():
        return trimmed.ljust(4)
    return trimmed.rjust(4)


def split_pdb_file_by_chains(
    pdb_path: Path,
    receptor_chain: str,
    binder_chain: str,
    receptor_out: Path,
    binder_out: Path,
) -> None:
    receptor_lines: list[str] = []
    binder_lines: list[str] = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain_id = line[21:22].strip()
        if chain_id == receptor_chain:
            receptor_lines.append(line)
        elif chain_id == binder_chain:
            binder_lines.append(line)
    if not receptor_lines:
        raise PreprocessError(f"Hydrogenated PDB is missing receptor chain {receptor_chain}")
    if not binder_lines:
        raise PreprocessError(f"Hydrogenated PDB is missing binder chain {binder_chain}")
    receptor_out.write_text("\n".join(receptor_lines + ["TER", "END"]) + "\n", encoding="utf-8")
    binder_out.write_text("\n".join(binder_lines + ["TER", "END"]) + "\n", encoding="utf-8")
