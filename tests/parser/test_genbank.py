from pathlib import Path

from pycirclizely.parser import Genbank


def test_parse_complete_genome(prokaryote_testdata_dir: Path):
    """Test parse complete genome"""
    gbk_file = prokaryote_testdata_dir / "enterobacteria_phage.gbk"
    gbk = Genbank(gbk_file)
    seqid = "NC_000902.1"
    assert gbk.name == "enterobacteria_phage"
    max_genome_size = 60942
    assert gbk.range_size == max_genome_size
    assert gbk.genome_length == max_genome_size
    assert gbk.full_genome_length == max_genome_size
    assert gbk.get_seqid2size() == {seqid: max_genome_size}


def test_parse_contig_genomes(prokaryote_testdata_dir: Path):
    """Test parse contig genomes"""
    gbk_file = prokaryote_testdata_dir / "mycoplasma_alvi.gbk.gz"
    gbk = Genbank(gbk_file)
    seqid2size = {
        "NZ_JNJU01000001.1": 264665,
        "NZ_JNJU01000002.1": 190782,
        "NZ_KL370824.1": 158240,
        "NZ_KL370825.1": 155515,
        "NZ_JNJU01000007.1": 67647,
        "NZ_JNJU01000008.1": 2683,
        "NZ_JNJU01000009.1": 1108,
    }
    size_list = list(seqid2size.values())

    assert gbk.name == "mycoplasma_alvi"
    assert gbk.range_size == size_list[0]
    assert gbk.get_seqid2size() == seqid2size

    seqid2cds_features = gbk.get_seqid2features()
    first_contig_cds_features = list(seqid2cds_features.values())[0]
    assert len(first_contig_cds_features) == 204

    seqid2trna_features = gbk.get_seqid2features("tRNA")
    first_contig_trna_features = list(seqid2trna_features.values())[0]
    assert len(first_contig_trna_features) == 12


def test_calc_gc(prokaryote_testdata_dir: Path):
    """Test GCskew, GCcontent calculation"""
    gbk_file = prokaryote_testdata_dir / "enterobacteria_phage.gbk"
    gbk = Genbank(gbk_file)
    gbk.calc_gc_skew()
    gbk.calc_gc_content()

    # Set user specified short sequence (< 1000)
    # Check if zero length error of window_size & step_size is not occured
    seq = "ATGC" * 100
    gbk.calc_gc_skew(seq=seq)
    gbk.calc_gc_content(seq=seq)


def test_write_cds_fasta(prokaryote_testdata_dir: Path, fig_outfile: Path):
    """Test `gbk.write_cds_fasta()` (only check properly output file)"""
    gbk_file = prokaryote_testdata_dir / "enterobacteria_phage.gbk"
    gbk = Genbank(gbk_file)
    gbk.write_cds_fasta(fig_outfile)
    assert fig_outfile.exists()


def test_genome_fasta(prokaryote_testdata_dir: Path, fig_outfile: Path):
    """Test `gbk.write_genome_fasta()` (only check properly output file)"""
    gbk_file = prokaryote_testdata_dir / "enterobacteria_phage.gbk"
    gbk = Genbank(gbk_file)
    gbk.write_genome_fasta(fig_outfile)
    assert fig_outfile.exists()
