"""The SQL files applied by build_BFD_duckDB.sh declare their own default
data path on line 1, as '-- TABLES_DIR=<path>'.

The builder rewrites that path to the one the build actually uses with a global
substitution of the line-1 string:

    sed "s#${DEFAULT_P}#${NEWPATH}#g" "$SQL_DIR/$f"

So line 1 is not documentation -- it is the substitution key. If the body
references a path that line 1 does not spell exactly, that reference survives
the rewrite and the built database silently reads whatever is at the stale
path. That happened when these files moved from BFD_dashboard into Fungi_BFD
(2026-11-19): line 1 of bfd_traits.sql was repointed at lib/refs while its four
read_csv() calls still named BFD_dashboard/config/refs, and the traits tests
read the real FungalTraits table instead of their fixture.
"""
import re
from pathlib import Path

import pytest

SQL_DIR = Path(__file__).resolve().parents[2] / "sql"
FILES = ["bfd_compat_views.sql", "bfd_extra_views.sql"]
KEYS = {}  # all use TABLES_DIR

# read_parquet('/abs/path/...'), read_csv('/abs/path/...'), read_csv_auto(...)
_READ = re.compile(r"read_(?:parquet|csv|csv_auto)\(\s*'([^']+)'")


def _header(path: Path) -> tuple[str, str]:
    first = path.read_text().splitlines()[0]
    m = re.fullmatch(r"-- (TABLES_DIR|REF_DIR)=(\S+)", first)
    assert m, f"line 1 of {path.name} must be '-- TABLES_DIR=<path>' or '-- REF_DIR=<path>', got {first!r}"
    return m.group(1), m.group(2)


@pytest.mark.parametrize("name", FILES)
def test_header_key_is_the_expected_one(name):
    key, _ = _header(SQL_DIR / name)
    assert key == KEYS.get(name, "TABLES_DIR")


@pytest.mark.parametrize("name", FILES)
def test_every_absolute_read_path_is_under_the_declared_prefix(name):
    path = SQL_DIR / name
    _, declared = _header(path)
    bad = [p for p in _READ.findall(path.read_text())
           if p.startswith("/") and not p.startswith(declared + "/")]
    assert not bad, (
        f"{name}: these absolute read paths are not under the line-1 prefix "
        f"{declared!r}, so the builder's substitution will not rewrite them: {bad}"
    )


@pytest.mark.parametrize("name", FILES)
def test_declared_prefix_is_used_or_the_file_reads_nothing(name):
    """Either the body actually uses the declared prefix, or it contains no
    absolute read at all.

    bfd_extra_views.sql is the second case: every one of its views is built on
    other views, so it never names a parquet path and its '-- TABLES_DIR='
    header is vestigial. The builder still requires the header (it errors when
    line 1 lacks the key) and the substitution is simply a no-op. What must
    never happen is the third case -- a body that reads from an absolute path
    the header does not cover -- which is what the test above catches.
    """
    path = SQL_DIR / name
    _, declared = _header(path)
    text = path.read_text()
    body = "\n".join(text.splitlines()[1:])
    absolute_reads = [p for p in _READ.findall(text) if p.startswith("/")]
    if not absolute_reads:
        return  # nothing to rewrite; header is inert by design
    assert declared in body, (
        f"{name}: line 1 declares {declared!r} and the body has absolute reads "
        f"{absolute_reads}, but none of them sit under that prefix"
    )
