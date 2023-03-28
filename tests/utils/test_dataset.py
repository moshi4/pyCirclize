from pathlib import Path
from urllib.request import urlopen

import pytest

from pycirclize.utils import (
    fetch_genbank_by_accid,
    load_eukaryote_example_dataset,
    load_example_image_file,
    load_prokaryote_example_file,
)


def check_network_conn(url: str = "https://github.com/moshi4/pyCirclize") -> bool:
    """Check network connection

    Parameters
    ----------
    url : str, optional
        Check target URL

    Returns
    -------
    check_result : bool
        Network connection check result
    """
    try:
        with urlopen(url) as _:
            return True
    except Exception:
        return False


@pytest.mark.skipif(
    condition=not check_network_conn(),
    reason="No network connection.",
)
def test_load_prokaryote_example_file():
    """Test `load_prokaryote_example_file()`"""
    gbk_file = load_prokaryote_example_file("enterobacteria_phage.gbk")
    assert gbk_file.exists()


@pytest.mark.skipif(
    condition=not check_network_conn(),
    reason="No network connection.",
)
def test_load_eukaryote_example_dataset():
    """Test `load_eukaryote_example_dataset()`"""
    bed_file, cytoband_file, _ = load_eukaryote_example_dataset("hg38")
    assert bed_file.exists()
    assert cytoband_file.exists()


@pytest.mark.skipif(
    condition=not check_network_conn(),
    reason="No network connection.",
)
def test_fetch_genbank_by_accid(tmp_path: Path):
    accid = "JX128258.1"
    # Case1. Download as textio
    _ = fetch_genbank_by_accid(accid)
    # Case2. Download as file
    gbk_outfile = tmp_path / "out.gbk"
    fetch_genbank_by_accid(accid, gbk_outfile=gbk_outfile)
    assert gbk_outfile.exists()


def test_load_example_image_file():
    """Test `load_example_image_file()`"""
    # 1. Normal scenario
    image_file = load_example_image_file("python_logo.png")
    assert image_file.exists()

    # 2. Exception scenario
    with pytest.raises(FileNotFoundError):
        load_example_image_file("noexists.png")
