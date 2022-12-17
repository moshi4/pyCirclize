from pathlib import Path

import pytest


@pytest.fixture
def testdata_dir() -> Path:
    """Testdata directory fixture"""
    return Path(__file__).parent / "testdata"
