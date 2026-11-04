import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).parents[1]
SPEC = importlib.util.spec_from_file_location("export_liver_tumor_dashboard", ROOT / "scripts/export_liver_tumor_dashboard.py")
assert SPEC and SPEC.loader
EXPORTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORTER)


def test_liver_tumour_dashboard_export_is_current():
    expected = EXPORTER.build_dashboard_data()
    committed = json.loads((ROOT / "docs/data/liver_tumor_dashboard.json").read_text())
    assert committed == expected
    assert expected["schema_version"] == 1
    assert expected["research_only"] is True
    assert set(expected["scenarios"]) == {"untreated", "antiangiogenic", "low_perfusion",
                                          "portal_restriction", "arterial_restriction"}
    assert all(len(records) == 61 for records in expected["scenarios"].values())
    assert set(expected["immune_coupling"]) == {"control", "anti_pd1"}
    assert all(len(records) == 43 for records in expected["immune_coupling"].values())
    treated = expected["immune_coupling"]["anti_pd1"]
    control = expected["immune_coupling"]["control"]
    assert max(row["pd1_occupancy"] for row in treated) > max(
        row["pd1_occupancy"] for row in control)
    assert len(expected["tdln_zones"]) == 6
    assert set(expected["liver_driven_tdln_spatial"]) == {"control", "anti_pd1"}
    spatial_control = expected["liver_driven_tdln_spatial"]["control"]
    assert len(spatial_control) == 15
    assert spatial_control[-1]["values"]["antigen"]["subcapsular_sinus"] > 0
    spatial_treated = expected["liver_driven_tdln_spatial"]["anti_pd1"]
    assert spatial_treated[-1]["values"]["effector"]["efferent_lymph"] > (
        spatial_control[-1]["values"]["effector"]["efferent_lymph"])


def test_liver_tumour_page_is_animated_accessible_and_local():
    page = (ROOT / "docs/liver-tumor.html").read_text()
    for marker in ('id="play"', 'id="day-slider"', 'aria-live="polite"', '<title>', '<desc ',
                   'data/liver_tumor_dashboard.json', '@keyframes breathe', '@keyframes flow',
                   '@media(prefers-reduced-motion:reduce)', '正常肝细胞', '肿瘤细胞',
                   'id="scenario-anti_pd1"', 'id="tdln-node"', 'immune_coupling'):
        assert marker in page
    assert 'id="scenario-portal_restriction"' in page
    assert 'id="scenario-arterial_restriction"' in page
    assert "研究用途" in page
