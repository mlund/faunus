import pathlib
import unittest

import yaml
from jsonschema import Draft7Validator


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCHEMA = yaml.safe_load((ROOT / "docs" / "schema.yml").read_text())
VALIDATOR = Draft7Validator(SCHEMA)


def assert_valid(testcase: unittest.TestCase, obj: dict):
    errors = list(VALIDATOR.iter_errors(obj))
    testcase.assertFalse(errors, "\n".join(str(e) for e in errors))


def assert_invalid(testcase: unittest.TestCase, obj: dict):
    errors = list(VALIDATOR.iter_errors(obj))
    testcase.assertTrue(errors, "Expected validation to fail")


class SchemaExamplesTest(unittest.TestCase):
    def base_cfg(self):
        """Minimal configuration mirroring minimal.yml."""
        return {
            "insertmolecules": [{"salt": {"N": 20}}],
            "temperature": 300,
        }

    def base_no_scatter(self):
        cfg = self.base_cfg()
        cfg["analysis"] = [{"systemenergy": {"file": "energy.dat", "nstep": 100}}]
        return cfg

    def test_minimal_valid(self):
        assert_valid(self, self.base_cfg())

    # analysis -> scatter
    def test_scatter_debye_with_qgrid_valid(self):
        cfg = self.base_cfg()
        cfg["analysis"] = [
            {
                "scatter": {
                    "molecules": ["salt"],
                    "nstep": 10,
                    "qmin": 0,
                    "qmax": 10,
                    "dq": 0.1,
                    "scheme": "debye",
                }
            }
        ]
        assert_valid(self, cfg)

    def test_scatter_qgrid_requires_debye(self):
        cfg = self.base_cfg()
        cfg["analysis"] = [{"scatter": {"molecules": ["salt"], "nstep": 10, "qmin": 0}}]
        assert_invalid(self, cfg)

    def test_scatter_explicit_with_pmax_valid(self):
        cfg = self.base_cfg()
        cfg["analysis"] = [{"scatter": {"molecules": ["salt"], "nstep": 10, "pmax": 5}}]
        assert_valid(self, cfg)

    def test_scatter_pmax_forbidden_in_debye(self):
        cfg = self.base_cfg()
        cfg["analysis"] = [
            {
                "scatter": {
                    "molecules": ["salt"],
                    "nstep": 10,
                    "scheme": "debye",
                    "pmax": 1,
                }
            }
        ]
        assert_invalid(self, cfg)

    # moves -> transrot
    def test_transrot_missing_molecule_invalid(self):
        cfg = self.base_cfg()
        cfg["moves"] = [{"transrot": {}}]
        assert_invalid(self, cfg)

    # mcloop
    def test_mcloop_missing_micro_invalid(self):
        cfg = self.base_no_scatter()
        cfg["mcloop"] = {"macro": 10}
        assert_invalid(self, cfg)

    # geometry
    def test_geometry_missing_length_invalid(self):
        cfg = self.base_no_scatter()
        cfg["geometry"] = {"type": "cuboid"}
        assert_invalid(self, cfg)

    # insertmolecules
    def test_insertmolecules_missing_count_invalid(self):
        cfg = self.base_no_scatter()
        cfg["insertmolecules"] = [{"salt": {}}]
        assert_invalid(self, cfg)


if __name__ == "__main__":
    unittest.main()
