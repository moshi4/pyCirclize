from __future__ import annotations

import bz2
import gzip
import zipfile
from functools import lru_cache
from io import TextIOWrapper
from pathlib import Path
from typing import Any

import numpy as np
from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation, Seq, SeqFeature
from Bio.SeqRecord import SeqRecord


class Genbank:
    """Genbank Parser Class"""

    def __init__(
        self,
        gbk_source: str | Path | TextIOWrapper,
        name: str | None = None,
        reverse: bool = False,
        min_range: int | None = None,
        max_range: int | None = None,
    ):
        """
        Parameters
        ----------
        gbk_source : str | Path | TextIOWrapper
            Genbank file or source
            (`*.gz`, `*.bz2`, `*.zip` compressed file can be readable)
        name : str | None
            name (If None, `file name` or `record name` is set)
        reverse : bool, optional
            If True, reverse complement genome is used
        min_range : int | None, optional
            Min range to be extracted (Default: `0`)
        max_range : int | None, optional
            Max range to be extracted (Default: `genome length`)
        """
        self._gbk_source = gbk_source
        self._name = name
        self._records = self._parse_gbk_source(gbk_source)
        self.reverse = reverse
        self.min_range = 0 if min_range is None else min_range
        self.max_range = self.full_genome_length if max_range is None else max_range

        if not 0 <= self.min_range <= self.max_range <= self.full_genome_length:
            err_msg = f"min_range={min_range}, max_range={max_range} is invalid. \n"
            err_msg += "Range must be "
            err_msg += f"'0 <= min_range <= max_range <= {self.full_genome_length}'"
            raise ValueError(err_msg)

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

    @property
    def name(self) -> str:
        """Name"""
        if self._name is not None:
            return self._name
        if isinstance(self._gbk_source, (str, Path)):
            gbk_file = Path(self._gbk_source)
            if gbk_file.suffix in (".gz", ".bz2", ".zip"):
                return gbk_file.with_suffix("").with_suffix("").name
            else:
                return gbk_file.with_suffix("").name
        elif isinstance(self._gbk_source, TextIOWrapper):
            return self._records[0].name
        else:
            raise NotImplementedError()

    @property
    def records(self) -> list[SeqRecord]:
        """Genbank records"""
        if self.reverse:
            return list(reversed([r.reverse_complement() for r in self._records]))
        else:
            return self._records

    @property
    def full_genome_length(self) -> int:
        """Full genome sequence length"""
        return len(self.full_genome_seq)

    @property
    def genome_length(self) -> int:
        """Range genome sequence length (Same as `range_size`)"""
        return len(self.genome_seq)

    @property
    def range_size(self) -> int:
        """Range size (`max_range - min_range`)"""
        return self.max_range - self.min_range

    @property
    def full_genome_seq(self) -> str:
        """Full genome sequence"""
        return "".join(str(r.seq) for r in self.records)

    @property
    def genome_seq(self) -> str:
        """Range genome sequence"""
        seq = "".join(str(r.seq) for r in self.records)
        return seq[self.min_range : self.max_range]

    @lru_cache(maxsize=None)
    def calc_genome_gc_content(self) -> float:
        """Calculate genome GC content"""
        return SeqUtils.GC(self.genome_seq)

    def calc_gc_skew(
        self,
        window_size: int | None = None,
        step_size: int | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculate GC skew in sliding window

        Parameters
        ----------
        window_size : int | None, optional
            Window size (Default: `genome_size / 1000`)
        step_size : int | None, optional
            Step size (Default: `genome_size / 2000`)

        Returns
        -------
        gc_skew_result_tuple : tuple[np.ndarray, np.ndarray]
            Position list & GC skew list
        """
        pos_list, gc_skew_list = [], []
        seq = self.genome_seq
        if window_size is None:
            window_size = int(len(seq) / 1000)
        if step_size is None:
            step_size = int(len(seq) / 2000)
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

        return (np.array(pos_list), np.array(gc_skew_list))

    def calc_gc_content(
        self,
        window_size: int | None = None,
        step_size: int | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculate GC content in sliding window

        Parameters
        ----------
        window_size : int | None, optional
            Window size (Default: `genome_size / 1000`)
        step_size : int | None, optional
            Step size (Default: `genome_size / 2000`)

        Returns
        -------
        gc_content_result_tuple : tuple[np.ndarray, np.ndarray]
            Position list & GC content list
        """
        pos_list, gc_content_list = [], []
        seq = self.genome_seq
        if window_size is None:
            window_size = int(len(seq) / 1000)
        if step_size is None:
            step_size = int(len(seq) / 2000)
        pos_list = list(range(0, len(seq), step_size)) + [len(seq)]
        for pos in pos_list:
            window_start_pos = pos - int(window_size / 2)
            window_end_pos = pos + int(window_size / 2)
            window_start_pos = 0 if window_start_pos < 0 else window_start_pos
            window_end_pos = len(seq) if window_end_pos > len(seq) else window_end_pos

            subseq = seq[window_start_pos:window_end_pos]
            gc_content_list.append(SeqUtils.GC(subseq))

        return (np.array(pos_list), np.array(gc_content_list))

    def extract_features(
        self,
        feature_type: str = "CDS",
        target_strand: int | None = None,
        fix_position: bool = False,
        allow_partial: bool = False,
        pseudogene: bool = False,
    ) -> list[SeqFeature]:
        """Extract features within min-max range

        Parameters
        ----------
        feature_type : str, optional
            Extract feature type
        target_strand : int | None, optional
            Extract target strand
        fix_position : bool, optional
            If True, fix feature start & end position by specified min_range parameter
            (fixed_start = start - min_range, fixed_end = end - min_range)
        allow_partial : bool, optional
            If True, allow extraction of features that are partially included in range
        pseudogene : bool, optional
            If True and `feature_type='CDS'`, only extract CDS features with
            `/pseudo` or `/pseudogene` qualifiers.

        Returns
        -------
        features : list[SeqFeature]
            Extracted features
        """
        extract_features = []
        min_range, max_range = self.min_range, self.max_range
        base_len = 0
        for record in self.records:
            features = [f for f in record.features if f.type == feature_type]
            for f in features:
                if feature_type == "CDS":
                    if pseudogene:
                        # Only extract pseudogene
                        qual = f.qualifiers
                        if "pseudo" not in qual and "pseudogene" not in qual:
                            continue
                    else:
                        # Exclude pseudogene (no translated gene)
                        translation = f.qualifiers.get("translation", [None])[0]
                        if translation is None:
                            continue
                if f.strand == -1:
                    # Handle rare case (complement & join)
                    # Found in NC_00913 protein_id=NP_417367.1
                    start = self._to_int(f.location.parts[-1].start) + base_len
                    end = self._to_int(f.location.parts[0].end) + base_len
                else:
                    start = self._to_int(f.location.parts[0].start) + base_len
                    end = self._to_int(f.location.parts[-1].end) + base_len
                # Restrict features in range
                if allow_partial:
                    if (
                        not min_range <= start <= max_range
                        and not min_range <= end <= max_range
                    ):
                        # Ignore completely out of range features
                        continue
                    # If partially within range, fix position to within range
                    if start <= min_range <= end <= max_range:
                        start = min_range
                    if min_range <= start <= max_range <= end:
                        end = max_range
                else:
                    if not min_range <= start <= end <= max_range:
                        # Ignore out of range features
                        continue
                # Extract only target strand feature
                if target_strand is not None and f.strand != target_strand:
                    continue
                # Fix start & end position by min_range
                if fix_position:
                    start -= min_range
                    end -= min_range

                extract_features.append(
                    SeqFeature(
                        location=FeatureLocation(start, end, f.strand),
                        type=f.type,
                        qualifiers=f.qualifiers,
                    ),
                )
            base_len += len(record.seq)

        return extract_features

    def write_cds_fasta(
        self,
        fasta_outfile: str | Path,
        seqtype: str = "protein",
        fix_position: bool = False,
        allow_partial: bool = False,
    ):
        """Write CDS protein features fasta file

        Parameters
        ----------
        fasta_outfile : str | Path
            CDS fasta file
        seqtype : str, optional
            Sequence type (`protein`|`nucleotide`)
        fix_position : bool, optional
            If True, fix feature start & end position by specified min_range parameter
            (fixed_start = start - min_range, fixed_end = end - min_range)
        allow_partial : bool, optional
            If True, features that are partially included in range are also extracted
        """
        features = self.extract_features("CDS", None, fix_position, allow_partial)
        cds_seq_records: list[SeqRecord] = []
        for idx, feature in enumerate(features, 1):
            qualifiers = feature.qualifiers
            protein_id = qualifiers.get("protein_id", [None])[0]
            product = qualifiers.get("product", [""])[0]
            translation = qualifiers.get("translation", [None])[0]

            start = self._to_int(feature.location.start)
            end = self._to_int(feature.location.end)
            strand = -1 if feature.strand == -1 else 1

            location_id = f"|{start}_{end}_{strand}|"
            if protein_id is None:
                seq_id = f"GENE{idx:06d}{location_id}"
            else:
                seq_id = f"GENE{idx:06d}_{protein_id}{location_id}"

            if seqtype == "protein":
                seq = Seq(translation)
            elif seqtype == "nucleotide":
                seq = Seq(feature.location.extract(self.genome_seq))
            else:
                raise ValueError(f"seqtype='{seqtype}' is invalid.")

            cds_seq_record = SeqRecord(seq=seq, id=seq_id, description=product)
            cds_seq_records.append(cds_seq_record)

        SeqIO.write(cds_seq_records, fasta_outfile, "fasta-2line")

    def write_genome_fasta(
        self,
        outfile: str | Path,
    ) -> None:
        """Write genome fasta file

        Parameters
        ----------
        outfile : str | Path
            Output genome fasta file
        """
        write_seq = self.genome_seq
        with open(outfile, "w") as f:
            f.write(f">{self.name}\n{write_seq}\n")

    def _to_int(self, value: Any) -> int:
        """Convert to int (Required for AbstractPosition|ExactPosition)"""
        return int(str(value).replace("<", "").replace(">", ""))

    def __str__(self):
        name = str(self._gbk_source)
        min_max_range = f"{self.min_range:,} - {self.max_range:,} bp"
        return f"{name} ({min_max_range})"
