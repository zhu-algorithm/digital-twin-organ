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
    assert set(expected["scenarios"]) == {"untreated", "antiangiogenic", "low_perfusion"}
    assert all(len(records) == 61 for records in expected["scenarios"].values())


def test_liver_tumour_page_is_animated_accessible_and_local():
    page = (ROOT / "docs/liver-tumor.html").read_text()
    for marker in ('id="play"', 'id="day-slider"', 'aria-live="polite"', '<title>', '<desc ',
                   'data/liver_tumor_dashboard.json', '@keyframes breathe', '@keyframes flow',
                   '@media(prefers-reduced-motion:reduce)', '正常肝细胞', '肿瘤细胞'):
        assert marker in page
    assert "研究用途" in page
