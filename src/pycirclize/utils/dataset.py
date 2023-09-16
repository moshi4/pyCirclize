from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from io import StringIO, TextIOWrapper
from pathlib import Path
from urllib.request import urlretrieve

from Bio import Entrez

from pycirclize import config


def load_prokaryote_example_file(
    filename: str,
    cache_dir: str | Path | None = None,
    overwrite_cache: bool = False,
) -> Path:
    """Load pycirclize example Genbank or GFF file

    Load example file from <https://github.com/moshi4/pycirclize-data/>
    and cache file in local directory (Default: `~/.cache/pycirclize/`).

    List of example Genbank or GFF filename

    - `enterobacteria_phage.[gbk|gff]`
    - `mycoplasma_alvi.[gbk|gff]`
    - `escherichia_coli.[gbk|gff].gz`

    Parameters
    ----------
    filename : str
        Genbank or GFF filename (e.g. `enterobacteria_phage.gff`)
    cache_dir : str | Path | None, optional
        Output cache directory (Default: `~/.cache/pycirclize/`)
    overwrite_cache : bool, optional
        If True, overwrite cache file.
        Assumed to be used when cache file is corrupt.

    Returns
    -------
    file_path : Path
        Genbank or GFF file
    """
    # Check specified filename exists or not
    if filename not in config.PROKARYOTE_FILES:
        err_msg = f"{filename=} not found."
        raise ValueError(err_msg)

    # Cache local directory
    if cache_dir is None:
        package_name = __name__.split(".")[0]
        cache_base_dir = Path.home() / ".cache" / package_name
        cache_dir = cache_base_dir / "prokaryote"
        os.makedirs(cache_dir, exist_ok=True)
    else:
        cache_dir = Path(cache_dir)
        if not cache_dir.exists():
            raise ValueError(f"{cache_dir=} not exists.")

    # Download file
    file_url = config.GITHUB_DATA_URL + f"prokaryote/{filename}"
    file_path = cache_dir / filename
    if overwrite_cache or not file_path.exists():
        urlretrieve(file_url, file_path)

    return file_path


def load_eukaryote_example_dataset(
    name: str = "hg38",
    cache_dir: str | Path | None = None,
    overwrite_cache: bool = False,
) -> tuple[Path, Path, list[ChrLink]]:
    """Load pycirclize eukaryote example dataset

    Load example file from <https://github.com/moshi4/pycirclize-data/>
    and cache file in local directory (Default: `~/.cache/pycirclize/`).

    List of dataset contents (download from UCSC)

    1. Chromosome BED file (e.g. `chr1 0 248956422`)
    2. Cytoband file (e.g. `chr1 0 2300000 p36.33 gneg`)
    3. Chromosome links (e.g. `chr1 1000 4321 chr3 8000 5600`)

    Parameters
    ----------
    name : str, optional
        Dataset name (`hg38` or `mm10`)
    cache_dir : str | Path | None, optional
        Output cache directory (Default: `~/.cache/pycirclize/`)
    overwrite_cache : bool
        If True, overwrite cache dataset.
        Assumed to be used when cache dataset is corrupt.

    Returns
    -------
    chr_bed_file, cytoband_file, chr_links : tuple[Path, Path, list[ChrLink]]
        BED file, Cytoband file, Chromosome links
    """
    # Check specified name dataset exists or not
    if name not in config.EUKARYOTE_DATASET:
        raise ValueError(f"{name=} dataset not found.")

    # Dataset cache local directory
    if cache_dir is None:
        package_name = __name__.split(".")[0]
        cache_base_dir = Path.home() / ".cache" / package_name
        cache_dir = cache_base_dir / "eukaryote" / name
        os.makedirs(cache_dir, exist_ok=True)
    else:
        cache_dir = Path(cache_dir)
        if not cache_dir.exists():
            raise ValueError(f"{cache_dir=} not exists.")

    # Download & cache dataset
    eukaryote_files: list[Path] = []
    chr_links: list[ChrLink] = []
    for filename in config.EUKARYOTE_DATASET[name]:
        file_url = config.GITHUB_DATA_URL + f"eukaryote/{name}/{filename}"
        file_path = cache_dir / filename
        if overwrite_cache or not file_path.exists():
            urlretrieve(file_url, file_path)
        if str(file_path).endswith("link.tsv"):
            chr_links = ChrLink.load(file_path)
        else:
            eukaryote_files.append(file_path)

    return eukaryote_files[0], eukaryote_files[1], chr_links


def load_example_image_file(filename: str) -> Path:
    """Load example image file from local package data

    e.g. `python_logo.png`

    Parameters
    ----------
    filename : str
        Image file name

    Returns
    -------
    image_file_path : Path
        Image file path
    """
    image_dir = Path(__file__).parent / "images"
    image_filenames = [f.name for f in image_dir.glob("*.png")]

    if filename.lower() in image_filenames:
        return image_dir / filename.lower()
    else:
        err_msg = f"{filename=} is not found.\n"
        err_msg += f"Available filenames = {image_filenames}"
        raise FileNotFoundError(err_msg)


def fetch_genbank_by_accid(
    accid: str,
    gbk_outfile: str | Path | None = None,
    email: str | None = None,
) -> TextIOWrapper:
    """Fetch genbank text by `Accession ID`

    Parameters
    ----------
    accid : str
        Accession ID
    gbk_outfile : str | Path | None, optional
        If file path is set, write fetch data to file
    email : str | None, optional
        Email address to notify download limitation (Required for bulk download)

    Returns
    -------
    TextIOWrapper
        Genbank data

    Examples
    --------
    >>> gbk_fetch_data = fetch_genbank_by_accid("NC_002483")
    >>> gbk = Genbank(gbk_fetch_data)
    """
    Entrez.email = "" if email is None else email
    gbk_fetch_data: TextIOWrapper = Entrez.efetch(
        db="nucleotide",
        id=accid,
        rettype="gbwithparts",
        retmode="text",
    )
    if gbk_outfile is not None:
        gbk_text = gbk_fetch_data.read()
        with open(gbk_outfile, "w") as f:
            f.write(gbk_text)
        gbk_fetch_data = StringIO(gbk_text)

    return gbk_fetch_data


@dataclass
class ChrLink:
    """Chromosome Link DataClass"""

    query_chr: str
    query_start: int
    query_end: int
    ref_chr: str
    ref_start: int
    ref_end: int

    @staticmethod
    def load(chr_link_file: str | Path) -> list[ChrLink]:
        """Load chromosome link file

        Parameters
        ----------
        chr_link_file : str | Path
            Chromosome link file

        Returns
        -------
        chr_link_list : list[ChrLink]
            Chromosome link list
        """
        chr_link_list = []
        with open(chr_link_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                qchr, qstart, qend = row[0], int(row[1]), int(row[2])
                rchr, rstart, rend = row[3], int(row[4]), int(row[5])
                chr_link_list.append(ChrLink(qchr, qstart, qend, rchr, rstart, rend))
        return chr_link_list
