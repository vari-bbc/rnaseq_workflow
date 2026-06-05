from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_CONFIG = Path("tests/test_config/config.yaml")


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def test_config(repo_root):
    yaml = pytest.importorskip("yaml")
    with (repo_root / TEST_CONFIG).open() as handle:
        return yaml.safe_load(handle)
