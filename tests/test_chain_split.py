from __future__ import annotations

import unittest
from pathlib import Path

from mdmm_filter_cuilab.preprocess import PreprocessError, parse_mmcif_atoms, split_chain_atoms


MMCIF_TEXT = """data_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 1 ? 11.104 13.207 2.100 1.00 10.00 1 GLY A N 1
ATOM 2 C CA . GLY A 1 1 ? 12.100 13.500 2.500 1.00 10.00 1 GLY A CA 1
ATOM 3 C C . GLY A 1 1 ? 12.700 12.300 3.200 1.00 10.00 1 GLY A C 1
ATOM 4 O O . GLY A 1 1 ? 12.200 11.200 2.900 1.00 10.00 1 GLY A O 1
ATOM 5 N N . ALA B 1 1 ? 15.104 16.207 4.100 1.00 10.00 1 ALA B N 1
ATOM 6 C CA . ALA B 1 1 ? 16.100 16.500 4.500 1.00 10.00 1 ALA B CA 1
ATOM 7 C C . ALA B 1 1 ? 16.700 15.300 5.200 1.00 10.00 1 ALA B C 1
ATOM 8 O O . ALA B 1 1 ? 16.200 14.200 4.900 1.00 10.00 1 ALA B O 1
#
"""


class ChainSplitTests(unittest.TestCase):
    def _write_case(self, name: str) -> Path:
        case_dir = Path("tests_artifacts")
        case_dir.mkdir(exist_ok=True)
        cif_path = case_dir / f"{name}.cif"
        cif_path.write_text(MMCIF_TEXT, encoding="utf-8")
        self.addCleanup(lambda: cif_path.unlink(missing_ok=True))
        self.addCleanup(self._cleanup_case_dir, case_dir)
        return cif_path

    @staticmethod
    def _cleanup_case_dir(case_dir: Path) -> None:
        try:
            case_dir.rmdir()
        except OSError:
            pass

    def test_parse_and_split(self) -> None:
        cif_path = self._write_case("model_001")
        atoms = parse_mmcif_atoms(cif_path)
        complex_atoms, receptor_atoms, binder_atoms = split_chain_atoms(atoms, "A", "B")
        self.assertEqual(len(complex_atoms), 8)
        self.assertEqual(len(receptor_atoms), 4)
        self.assertEqual(len(binder_atoms), 4)

    def test_missing_chain_raises(self) -> None:
        cif_path = self._write_case("model_missing_chain")
        atoms = parse_mmcif_atoms(cif_path)
        with self.assertRaises(PreprocessError):
            split_chain_atoms(atoms, "A", "C")


if __name__ == "__main__":
    unittest.main()
