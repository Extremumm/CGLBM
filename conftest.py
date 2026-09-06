"""Pytest configuration shared by every test under ``programs``."""

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent

# pycglbm is normally installed (`pip install -e PyCGLBM`). Fall back to the
# in-tree copy so that the suite also runs straight from a fresh clone.
try:  # pragma: no cover - trivial
    import pycglbm  # noqa: F401
except ImportError:  # pragma: no cover - trivial
    sys.path.insert(0, str(PROJECT_ROOT / "PyCGLBM"))

#: Markers a test function must carry, mirroring the directory it lives in.
MARKERS = ("unit_test", "validation", "verification")


def pytest_addoption(parser):
    parser.addoption("--runlong", action="store_true", default=False, help="run long tests")
    parser.addoption(
        "--allow_plot",
        action="store_true",
        default=False,
        help="keep the interactive matplotlib backend, allowing to plot",
    )


def pytest_configure(config):
    for marker in MARKERS:
        config.addinivalue_line("markers", f"{marker}: only run {marker} tests")
    config.addinivalue_line("markers", "long: test that runs a full simulation")

    if not config.getoption("--allow_plot"):
        import matplotlib

        matplotlib.use("Agg")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--runlong"):
        # --runlong not given in cli: skip long tests
        skip_long = pytest.mark.skip(reason="need --runlong option to run")
        for item in items:
            if "long" in item.keywords:
                item.add_marker(skip_long)
