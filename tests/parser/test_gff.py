from pathlib import Path

import pytest

from pycirclizely.parser import Gff


def test_parse_complete_genome(prokaryote_testdata_dir: Path):
    """Test parse complete genome"""
    gff_file = prokaryote_testdata_dir / "enterobacteria_phage.gff"
    gff = Gff(gff_file)
    seqid = "NC_000902.1"
    assert gff.target_seqid == seqid
    max_genome_size = 60942
    assert gff.range_size == max_genome_size
    assert gff.seq_region == (0, max_genome_size)
    assert gff.seqid_list == [seqid]
    assert gff.get_seqid2size() == {seqid: max_genome_size}


@pytest.mark.parametrize(
    "gff_filename",
    [
        ("mycoplasma_alvi.gff.gz"),
        ("mycoplasma_alvi_nocomment.gff.gz"),
    ],
)
def test_parse_contig_genomes(prokaryote_testdata_dir: Path, gff_filename: str):
    """Test parse contig genomes"""
    gff_file = prokaryote_testdata_dir / gff_filename
    gff = Gff(gff_file)
    seqid2size = {
        "NZ_JNJU01000001.1": 264665,
        "NZ_JNJU01000002.1": 190782,
        "NZ_KL370824.1": 158240,
        "NZ_KL370825.1": 155515,
        "NZ_JNJU01000007.1": 67647,
        "NZ_JNJU01000008.1": 2683,
        "NZ_JNJU01000009.1": 1108,
    }
    seqid_list = list(seqid2size.keys())
    size_list = list(seqid2size.values())

    assert gff.target_seqid == list(seqid2size.keys())[0]
    assert gff.range_size == size_list[0]
    assert gff.seq_region == (0, size_list[0])
    assert gff.seqid_list == seqid_list
    assert gff.get_seqid2size() == seqid2size

    seqid2cds_features = gff.get_seqid2features()
    first_contig_cds_features = list(seqid2cds_features.values())[0]
    assert len(first_contig_cds_features) == 204

    seqid2trna_features = gff.get_seqid2features("tRNA")
    first_contig_trna_features = list(seqid2trna_features.values())[0]
    assert len(first_contig_trna_features) == 12
