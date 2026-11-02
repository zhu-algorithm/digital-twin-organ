import json
from pathlib import Path

from scripts.export_tdln_dashboard import build_dashboard_data


ROOT = Path(__file__).parents[1]


def test_dashboard_export_schema_and_generated_file_are_current():
    expected = build_dashboard_data()
    committed = json.loads((ROOT / "docs/data/tdln_dashboard.json").read_text())
    assert committed == expected
    assert expected["schema_version"] == 1
    assert expected["research_only"] is True
    assert len(expected["zones"]) == 6
    assert set(expected["scenarios"]) == {"control", "anti_pd1"}
    assert expected["species"] == ["antigen", "mature_dc", "tpex", "effector"]


def test_dashboard_has_accessible_controls_and_local_model_data():
    page = (ROOT / "docs/tdln.html").read_text()
    for marker in ("id=\"scenario-control\"", "id=\"scenario-treated\"", "id=\"day-slider\"",
                   "aria-live=\"polite\"", "<title>", "<desc ", "data/tdln_dashboard.json"):
        assert marker in page
    assert "研究用途" in page
