import os
import sys
import unittest
import copy
from pathlib import Path
from unittest.mock import patch

import yaml
from jsonschema import Draft202012Validator, ValidationError

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from meew_metabolic_paths import (
    expand_path,
    resolve_cplex_path,
    resolve_resource_path,
)


class PathResolutionTests(unittest.TestCase):
    def test_resources_root_schema_is_optional_and_typed(self):
        schema = yaml.safe_load(
            (ROOT / "workflow" / "schemas" / "config.schema.yaml").read_text(
                encoding="utf-8"
            )
        )
        config = yaml.safe_load(
            (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
        )
        validator = Draft202012Validator(schema)
        validator.validate(config)
        self.assertEqual(config["resources_root"], "")
        legacy = copy.deepcopy(config)
        legacy.pop("resources_root")
        validator.validate(legacy)
        invalid = copy.deepcopy(config)
        invalid["resources_root"] = 42
        with self.assertRaises(ValidationError):
            validator.validate(invalid)

    def test_checkm2_precedence_and_spaces(self):
        with patch.dict(
            os.environ,
            {
                "HOME": "/tmp/metabolic home",
                "XDG_DATA_HOME": "/tmp/metabolic xdg",
                "MEEW_RESOURCES": "/tmp/metabolic resources",
                "CHECKM2DB": "/tmp/tool CheckM2",
            },
            clear=True,
        ):
            self.assertEqual(
                resolve_resource_path(
                    "/tmp/explicit CheckM2",
                    "/tmp/configured root",
                    "CHECKM2DB",
                    "checkm2_db",
                ),
                "/tmp/explicit CheckM2",
            )
            self.assertEqual(
                resolve_resource_path(
                    "", "/tmp/configured root", "CHECKM2DB", "checkm2_db"
                ),
                "/tmp/configured root/checkm2_db",
            )
            self.assertEqual(
                resolve_resource_path("", "", "CHECKM2DB", "checkm2_db"),
                "/tmp/tool CheckM2",
            )
            del os.environ["CHECKM2DB"]
            self.assertEqual(
                resolve_resource_path("", "", "CHECKM2DB", "checkm2_db"),
                "/tmp/metabolic resources/checkm2_db",
            )
            del os.environ["MEEW_RESOURCES"]
            self.assertEqual(
                resolve_resource_path("", "", "CHECKM2DB", "checkm2_db"),
                "/tmp/metabolic xdg/meew/checkm2_db",
            )
            del os.environ["XDG_DATA_HOME"]
            self.assertEqual(
                resolve_resource_path("", "", "CHECKM2DB", "checkm2_db"),
                "/tmp/metabolic home/.local/share/meew/checkm2_db",
            )
            self.assertEqual(
                expand_path("relative FAA/input.faa"), "relative FAA/input.faa"
            )

    def test_cplex_is_gated_by_solver_and_explicit_precedes_environment(self):
        with patch.dict(
            os.environ, {"CPLEX_STUDIO_DIR": "/tmp/environment CPLEX"}, clear=True
        ):
            self.assertEqual(resolve_cplex_path("scip", ""), "")
            self.assertEqual(resolve_cplex_path("cplex", ""), "/tmp/environment CPLEX")
            self.assertEqual(
                resolve_cplex_path("cplex", "/tmp/explicit CPLEX"),
                "/tmp/explicit CPLEX",
            )


if __name__ == "__main__":
    unittest.main()
