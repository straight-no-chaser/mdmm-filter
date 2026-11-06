from __future__ import annotations

import unittest

from mdmm_filter.config import EvalConfig


class ConfigTests(unittest.TestCase):
    def test_default_config_is_valid(self) -> None:
        cfg = EvalConfig()
        cfg.validate()
        self.assertEqual(cfg.mmpbsa_mode, "gbsa")

    def test_invalid_mode_raises(self) -> None:
        cfg = EvalConfig(mmpbsa_mode="bad")
        with self.assertRaises(ValueError):
            cfg.validate()

    def test_batch_validation_requires_input_dir(self) -> None:
        cfg = EvalConfig(input_dir="")
        with self.assertRaises(ValueError):
            cfg.validate_batch()


if __name__ == "__main__":
    unittest.main()
