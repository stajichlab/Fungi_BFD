"""MERGE_ASM_STATS and MERGE_BUSCO_GENOME pass --allow-shrink to publish_table.py when
params.allow_shrink is true (same wiring as MERGE_SAMPLES). Static check of the module text."""
import re
from pathlib import Path

MODULES = Path(__file__).resolve().parents[1] / "modules" / "BFD"
EXPR = (
    "def shrinkFlag = !params.merge_all.toBoolean() ? '' :\n"
    "        (params.allow_shrink.toBoolean() ? '--enforce-no-shrink --allow-shrink' : '--enforce-no-shrink')"
)


def _text(name):
    return (MODULES / name / "main.nf").read_text()


def test_modules_wire_allow_shrink():
    for name in ("MERGE_ASM_STATS", "MERGE_BUSCO_GENOME"):
        t = _text(name)
        assert EXPR in t, name
        assert "${shrinkFlag}" in t, name
        assert "merge_all.toBoolean() ? '--enforce-no-shrink' : ''" not in t, name


def test_publish_table_accepts_the_flags():
    src = (Path(__file__).resolve().parents[1] / "bin" / "publish_table.py").read_text()
    assert re.search(r'add_argument\("--allow-shrink"', src)
    assert re.search(r'add_argument\("--enforce-no-shrink"', src)
