from __future__ import annotations

import bz2
import gzip
import warnings
import zipfile
from collections import defaultdict
from io import StringIO, TextIOWrapper
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import Seq, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from numpy.typing import NDArray


class Genbank:
    """Genbank Parser Class"""

    def __init__(
        self,
        gbk_source: str | Path | TextIOWrapper | list[SeqRecord],
        *,
        name: str | None = None,
        min_range: None = None,
        max_range: None = None,
    ):
        """
        Parameters
        ----------
        gbk_source : str | Path | TextIOWrapper | list[SeqRecord]
            Genbank file or source
            (`*.gz`, `*.bz2`, `*.zip` compressed file can be readable)
        name : str | None, optional
            name (If None, `file name` or `record name` is set)
        min_range : None, optional
            No longer used. Left for backward compatibility.
        max_range : None, optional
            No longer used. Left for backward compatibility.
        """
        self._gbk_source = gbk_source
        if isinstance(gbk_source, (str, Path, StringIO, TextIOWrapper)):
            self._records = self._parse_gbk_source(gbk_source)
        else:
            self._records = gbk_source

        # Set genbank name
        if name is not None:
            self._name = name
        elif isinstance(self._gbk_source, (str, Path)):
            gbk_file = Path(self._gbk_source)
            if gbk_file.suffix in (".gz", ".bz2", ".zip"):
                self._name = gbk_file.with_suffix("").with_suffix("").name
            else:
                self._name = gbk_file.with_suffix("").name
        elif isinstance(self._gbk_source, (StringIO, TextIOWrapper)):
            self._name = self._records[0].name
        else:
            raise ValueError("Failed to get genbank name.")

        if min_range or max_range:
            warnings.warn("min_range & max_range is no longer used in Genbank parser.")

        if len(self.records) == 0:
            raise ValueError(f"Failed to parse {gbk_source} as Genbank file.")

    ############################################################
    # Property
    ############################################################

    @property
    def name(self) -> str:
        """Name"""
        return self._name

    @property
    def records(self) -> list[SeqRecord]:
        """Genbank records"""
        return self._records

    @property
    def genome_seq(self) -> str:
        """Genome sequence (only first record)"""
        return str(self.records[0].seq)

    @property
    def genome_length(self) -> int:
        """Genome length (only first record)"""
        return len(self.genome_seq)

    @property
    def range_size(self) -> int:
        """Same as `self.genome_length` (Left for backward compatibility)"""
        return self.genome_length

    @property
    def full_genome_seq(self) -> str:
        """Full genome sequence (concatenate all records)"""
        return "".join(str(r.seq) for r in self.records)

    @property
    def full_genome_length(self) -> int:
        """Full genome length (concatenate all records)"""
        return len(self.full_genome_seq)

    ############################################################
    # Public Method
    ############################################################

    def calc_genome_gc_content(self, seq: str | None = None) -> float:
        """Calculate genome GC content

        Parameters
        ----------
        seq : str | None, optional
            Sequence for GC content calculation (Default: `self.genome_seq`)

        Returns
        -------
        gc_content : float
            GC content
        """
        seq = self.genome_seq if seq is None else seq
        gc_content = SeqUtils.gc_fraction(seq) * 100
        return gc_content

    def calc_gc_skew(
        self,
        window_size: int | None = None,
        step_size: int | None = None,
        *,
        seq: str | None = None,
    ) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
        """Calculate GC skew in sliding window

        Parameters
        ----------
        window_size : int | None, optional
            Window size (Default: `genome_size / 500`)
        step_size : int | None, optional
            Step size (Default: `genome_size / 1000`)
        seq : str | None, optional
            Sequence for GCskew calculation (Default: `self.genome_seq`)

        Returns
        -------
        pos_list : NDArray[np.int64]
            Position list
        gc_skew_list : NDArray[np.float64]
            GC skew list
        """
        pos_list, gc_skew_list = [], []
        seq = self.genome_seq if seq is None else seq
        if window_size is None:
            window_size = int(len(seq) / 500)
        if step_size is None:
            step_size = int(len(seq) / 1000)
        if window_size == 0 or step_size == 0:
            window_size, step_size = len(seq), int(len(seq) / 2)
        pos_list = list(range(0, len(seq), step_size)) + [len(seq)]
        for pos in pos_list:
            window_start_pos = pos - int(window_size / 2)
            window_end_pos = pos + int(window_size / 2)
            window_start_pos = 0 if window_start_pos < 0 else window_start_pos
            window_end_pos = len(seq) if window_end_pos > len(seq) else window_end_pos

            subseq = seq[window_start_pos:window_end_pos]
            g = subseq.count("G") + subseq.count("g")
            c = subseq.count("C") + subseq.count("c")
            try:
                skew = (g - c) / float(g + c)
            except ZeroDivisionError:
                skew = 0.0
            gc_skew_list.append(skew)

        pos_list = np.array(pos_list).astype(np.int64)
        gc_skew_list = np.array(gc_skew_list).astype(np.float64)

        return pos_list, gc_skew_list

    def calc_gc_content(
        self,
        window_size: int | None = None,
        step_size: int | None = None,
        *,
        seq: str | None = None,
    ) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
        """Calculate GC content in sliding window

        Parameters
        ----------
        window_size : int | None, optional
            Window size (Default: `genome_size / 500`)
        step_size : int | None, optional
            Step size (Default: `genome_size / 1000`)
        seq : str | None, optional
            Sequence for GC content calculation (Default: `self.genome_seq`)

        Returns
        -------
        pos_list : NDArray[np.int64]
            Position list
        gc_content_list : NDArray[np.float64]
            GC content list
        """
        pos_list, gc_content_list = [], []
        seq = self.genome_seq if seq is None else seq
        if window_size is None:
            window_size = int(len(seq) / 500)
        if step_size is None:
            step_size = int(len(seq) / 1000)
        if window_size == 0 or step_size == 0:
            window_size, step_size = len(seq), int(len(seq) / 2)
        pos_list = list(range(0, len(seq), step_size)) + [len(seq)]
        for pos in pos_list:
            window_start_pos = pos - int(window_size / 2)
            window_end_pos = pos + int(window_size / 2)
            window_start_pos = 0 if window_start_pos < 0 else window_start_pos
            window_end_pos = len(seq) if window_end_pos > len(seq) else window_end_pos

            subseq = seq[window_start_pos:window_end_pos]
            gc_content = SeqUtils.gc_fraction(subseq) * 100
            gc_content_list.append(gc_content)

        pos_list = np.array(pos_list).astype(np.int64)
        gc_content_list = np.array(gc_content_list).astype(np.float64)

        return pos_list, gc_content_list

    def get_seqid2seq(self) -> dict[str, str]:
        """Get seqid & complete/contig/scaffold genome sequence dict

        Returns
        -------
        seqid2seq : dict[str, int]
            seqid & genome sequence dict
        """
        return {str(rec.id): rec.seq for rec in self.records}

    def get_seqid2size(self) -> dict[str, int]:
        """Get seqid & complete/contig/scaffold genome size dict

        Returns
        -------
        seqid2size : dict[str, int]
            seqid & genome size dict
        """
        return {seqid: len(seq) for seqid, seq in self.get_seqid2seq().items()}

    def get_seqid2features(
        self,
        feature_type: str | list[str] | None = "CDS",
        target_strand: int | None = None,
    ) -> dict[str, list[SeqFeature]]:
        """Get seqid & features in target seqid genome dict

        Parameters
        ----------
        feature_type : str | list[str] | None, optional
            Feature type (`CDS`, `gene`, `mRNA`, etc...)
            If None, extract regardless of feature type.
        target_strand : int | None, optional
            Extract target strand. If None, extract regardless of strand.

        Returns
        -------
        seqid2features : dict[str, list[SeqFeature]]
            seqid & features dict
        """
        if isinstance(feature_type, str):
            feature_type = [feature_type]

        seqid2features = defaultdict(list)
        for rec in self.records:
            feature: SeqFeature
            for feature in rec.features:
                strand = feature.location.strand
                if feature_type is not None and feature.type not in feature_type:
                    continue
                if target_strand is not None and strand != target_strand:
                    continue
                # Exclude feature which straddle genome start position
                if self._is_straddle_feature(feature):
                    continue
                start = int(feature.location.start)  # type: ignore
                end = int(feature.location.end)  # type: ignore
                seqid2features[rec.id].append(
                    SeqFeature(
                        location=SimpleLocation(start, end, strand),
                        type=feature.type,
                        qualifiers=feature.qualifiers,
                    ),
                )
        return seqid2features

    def extract_features(
        self,
        feature_type: str | list[str] | None = "CDS",
        *,
        target_strand: int | None = None,
        target_range: tuple[int, int] | None = None,
    ) -> list[SeqFeature]:
        """Extract features (only first record)

        Parameters
        ----------
        feature_type : str | list[str] | None, optional
            Feature type (`CDS`, `gene`, `mRNA`, etc...)
            If None, extract regardless of feature type.
        target_strand : int | None, optional
            Extract target strand. If None, extract regardless of strand.
        target_range : tuple[int, int] | None, optional
            Extract target range. If None, extract regardless of range.

        Returns
        -------
        features : list[SeqFeature]
            Extracted features
        """
        seqid2features = self.get_seqid2features(feature_type, target_strand)
        first_record_features = list(seqid2features.values())[0]
        if target_range:
            target_features = []
            for feature in first_record_features:
                start = int(feature.location.start)  # type: ignore
                end = int(feature.location.end)  # type: ignore
                if min(target_range) <= start <= end <= max(target_range):
                    target_features.append(feature)
            return target_features
        else:
            return first_record_features

    def write_cds_fasta(self, outfile: str | Path) -> None:
        """Write CDS fasta file

        Parameters
        ----------
        outfile : str | Path
            Output CDS fasta file
        """
        cds_records: list[SeqRecord] = []
        counter = 0
        seqid2cds_features = self.get_seqid2features(feature_type="CDS")
        for seqid, cds_features in seqid2cds_features.items():
            for cds_feature in cds_features:
                # Ignore no translation feature
                translation = cds_feature.qualifiers.get("translation", [None])[0]
                if translation is None:
                    continue
                # Get feature location
                start = int(cds_feature.location.start)  # type: ignore
                end = int(cds_feature.location.end)  # type: ignore
                strand = -1 if cds_feature.location.strand == -1 else 1
                # Set feature id
                location_id = f"|{seqid}|{start}_{end}_{strand}|"
                protein_id = cds_feature.qualifiers.get("protein_id", [None])[0]
                if protein_id is None:
                    feature_id = f"GENE{counter:06d}{location_id}"
                else:
                    feature_id = f"GENE{counter:06d}_{protein_id}{location_id}"
                counter += 1
                # Add SeqRecord of CDS feature
                seq = Seq(translation)
                product = cds_feature.qualifiers.get("product", [""])[0]
                seq_record = SeqRecord(seq, feature_id, description=product)
                cds_records.append(seq_record)
        # Write CDS file
        SeqIO.write(cds_records, outfile, "fasta-2line")

    def write_genome_fasta(self, outfile: str | Path) -> None:
        """Write genome fasta file

        Parameters
        ----------
        outfile : str | Path
            Output genome fasta file
        """
        with open(outfile, "w") as f:
            for seqid, seq in self.get_seqid2seq().items():
                f.write(f">{seqid}\n{seq}\n")

    ############################################################
    # Private Method
    ############################################################

    def _parse_gbk_source(
        self, gbk_source: str | Path | TextIOWrapper
    ) -> list[SeqRecord]:
        """Parse genbank source

        Parameters
        ----------
        gbk_source : str | Path | TextIOWrapper
            Genbank file or source

        Returns
        -------
        list[SeqRecord]
            Genbank SeqRecords
        """
        # Parse compressed file
        if isinstance(gbk_source, (str, Path)):
            if Path(gbk_source).suffix == ".gz":
                with gzip.open(gbk_source, mode="rt") as f:
                    return list(SeqIO.parse(f, "genbank"))
            elif Path(gbk_source).suffix == ".bz2":
                with bz2.open(gbk_source, mode="rt") as f:
                    return list(SeqIO.parse(f, "genbank"))
            elif Path(gbk_source).suffix == ".zip":
                with zipfile.ZipFile(gbk_source) as zip:
                    with zip.open(zip.namelist()[0]) as f:
                        return list(SeqIO.parse(TextIOWrapper(f), "genbank"))
        # Parse no compressed file or TextIOWrapper
        return list(SeqIO.parse(gbk_source, "genbank"))

    def _is_straddle_feature(self, feature: SeqFeature) -> bool:
        """Check target feature straddle genome start position or not

        Parameters
        ----------
        feature : SeqFeature
            Target feature

        Returns
        -------
        result : bool
            Check result
        """
        strand = feature.location.strand
        if strand == -1:
            start = int(feature.location.parts[-1].start)  # type: ignore
            end = int(feature.location.parts[0].end)  # type: ignore
        else:
            start = int(feature.location.parts[0].start)  # type: ignore
            end = int(feature.location.parts[-1].end)  # type: ignore
        return True if start > end else False

    def __str__(self):
        text = f"{self.name}: {len(self.records)} records\n"
        for num, (seqid, size) in enumerate(self.get_seqid2size().items(), 1):
            text += f"{num:02d}. {seqid} ({size:,} bp)\n"
        return text
