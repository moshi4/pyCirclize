from __future__ import annotations

import bz2
import gzip
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from io import TextIOWrapper
from pathlib import Path
from typing import Any, TextIO

from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature


class Gff:
    """GFF Parser Class"""

    def __init__(
        self,
        gff_file: str | Path,
        name: str | None = None,
        target_seqid: str | None = None,
        min_range: int | None = None,
        max_range: int | None = None,
    ):
        """
        Parameters
        ----------
        gff_file : str | Path
            GFF file (`*.gz`, `*.bz2`, `*.zip` compressed file can be readable)
        name : str | None, optional
            name (If None, `file name` is set)
        target_seqid : str | None, optional
            Target seqid to be extracted. If None, only first seqid record is extracted.
        min_range : int | None, optional
            Min range to be extracted.
            If None, appropriate value is taken from GFF records.
        max_range : int | None, optional
            Max range to be extracted.
            If None, appropriate value is taken from GFF records.
        """
        self._gff_file = Path(gff_file)
        self._name = name
        self._records, start, end = self._parse_gff(gff_file, target_seqid)
        self.min_range = start if min_range is None else min_range
        self.max_range = end if max_range is None else max_range
        self._seq_region = (start, end)

        if not 0 <= self.min_range <= self.max_range:
            err_msg = "Range must be '0 <= min_range <= max_range' "
            err_msg += f"(min_range={self.min_range}, max_range={self.max_range})"
            raise ValueError(err_msg)

    def _parse_gff(
        self,
        gff_file: str | Path,
        target_seqid: str | None,
    ) -> tuple[list[GffRecord], int, int]:
        """Parse GFF file

        Only parse target seqid record.
        If target_record is None, only parse first seqid record.

        Parameters
        ----------
        gff_file : str | Path
            GFF file
        target_seqid : str | None
            Target seqid to be extracted

        Returns
        -------
        gff_records, start, end : tuple[list[GffRecord], int, int]
            GFF record list, start, end
        """
        gff_file = Path(gff_file)
        if gff_file.suffix == ".gz":
            with gzip.open(gff_file, mode="rt") as f:
                gff_records, start, end = self._parse_gff_textio(f, target_seqid)
        elif gff_file.suffix == ".bz2":
            with bz2.open(gff_file, mode="rt") as f:
                gff_records, start, end = self._parse_gff_textio(f, target_seqid)
        elif gff_file.suffix == ".zip":
            with zipfile.ZipFile(gff_file) as zip:
                with zip.open(zip.namelist()[0]) as f:
                    io = TextIOWrapper(f)
                    gff_records, start, end = self._parse_gff_textio(io, target_seqid)
        else:
            with open(gff_file) as f:
                gff_records, start, end = self._parse_gff_textio(f, target_seqid)

        return gff_records, start, end

    def _parse_gff_textio(
        self,
        handle: TextIO,
        target_seqid: str | None = None,
    ) -> tuple[list[GffRecord], int, int]:
        """Parse GFF file TextIO

        Parameters
        ----------
        handle : TextIO
            GFF TextIO handle
        target_seqid : str | None, optional
            GFF target seqid

        Returns
        -------
        gff_records, start, end : tuple[list[GffRecord], int, int]
            GFF record list, start, end
        """
        gff_all_lines = handle.read().splitlines()
        gff_record_lines = filter(GffRecord.is_gff_line, gff_all_lines)
        gff_records = list(map(GffRecord.parse_gff_line, gff_record_lines))
        if len(gff_records) == 0:
            err_msg = f"Failed to parse '{self._gff_file}' as GFF file "
            raise ValueError(err_msg)

        seqid_list = list(dict.fromkeys([rec.seqid for rec in gff_records]))
        self._seqid_list = seqid_list
        if target_seqid is None:
            target_seqid = seqid_list[0]
        self._target_seqid = target_seqid
        if target_seqid not in seqid_list:
            err_msg = f"Not found target_seqid='{target_seqid}' in '{self._gff_file}'"
            raise ValueError(err_msg)
        gff_records = [rec for rec in gff_records if rec.seqid == target_seqid]

        # Try to get start-end region from '##sequence-region' annotation line
        start, end = None, None
        for line in gff_all_lines:
            if line.startswith("##sequence-region") and target_seqid in line:
                if len(line.split()) == 4:
                    start, end = line.split()[2:4]
                    start, end = int(start) - 1, int(end)
                    break
        if start is None or end is None:
            start, end = 0, max([r.end for r in gff_records])

        return gff_records, start, end

    @property
    def name(self) -> str:
        """Name"""
        if self._name is not None:
            return self._name
        if self._gff_file.suffix in (".gz", ".bz2", ".zip"):
            return self._gff_file.with_suffix("").with_suffix("").name
        else:
            return self._gff_file.with_suffix("").name

    @property
    def seq_region(self) -> tuple[int, int]:
        """GFF sequence-region start & end tuple

        If `##sequence-region` pragma is not found, seq_region=`(0, max_coords_value)`
        """
        return self._seq_region

    @property
    def records(self) -> list[GffRecord]:
        """GFF records"""
        return self._records

    @property
    def records_within_range(self) -> list[GffRecord]:
        """GFF records within min-max range"""
        target_gff_records: list[GffRecord] = []
        for rec in self.records:
            if rec.is_within_range(self.min_range, self.max_range):
                target_gff_records.append(rec)
        return target_gff_records

    @property
    def range_size(self) -> int:
        """Range size (`max_range - min_range`)"""
        return self.max_range - self.min_range

    @property
    def target_seqid(self) -> str:
        """Target seqid"""
        return self._target_seqid

    @property
    def seqid_list(self) -> list[str]:
        """seqid list"""
        return self._seqid_list

    def extract_features(
        self,
        feature_type: str = "CDS",
        target_strand: int | None = None,
        pseudogene: bool = False,
    ) -> list[SeqFeature]:
        """Extract features within min-max range

        Parameters
        ----------
        feature_type : str, optional
            Feature type (`CDS`, `gene`, `mRNA`, etc...)
        target_strand : int | None, optional
            Extract target strand
        pseudogene : bool, optional
            If True, `pseudo=`, `pseudogene=` tagged record only extract.
            If False, `pseudo=`, `pseudogene=` not tagged record only extract.

        Returns
        -------
        features : list[SeqFeature]
            Feature list
        """
        features: list[SeqFeature] = []
        for rec in self.records_within_range:
            if rec.type != feature_type:
                continue
            if target_strand is not None and rec.strand != target_strand:
                continue
            if pseudogene and ("pseudo" in rec.attrs or "pseudogene" in rec.attrs):
                features.append(rec.to_seq_feature())
            else:
                features.append(rec.to_seq_feature())

        return features

    def extract_exon_features(self, feature_type: str = "mRNA") -> list[SeqFeature]:
        """Extract exon structure features within min-max range

        Extract exons based on `parent feature` and `exon` ID-Parent relation

        Parameters
        ----------
        feature_type : str, optional
            Feature type (e.g. `mRNA`, `ncRNA` , etc...)

        Returns
        -------
        features : list[SeqFeature]
            Feature list
        """
        # Extract exon features by mRNA-exon relation
        parent_id = None
        parent_id2record: dict[str, GffRecord] = {}
        parent_id2exons: dict[str, list[GffRecord]] = defaultdict(list)
        for rec in self.records_within_range:
            if rec.type == feature_type:
                parent_id = rec.attrs.get("ID", [None])[0]
                if parent_id is None:
                    continue
                parent_id2record[parent_id] = rec
            if rec.type == "exon":
                if parent_id is not None and parent_id == rec.attrs["Parent"][0]:
                    parent_id2exons[parent_id].append(rec)

        # Set exon features
        exon_features: list[SeqFeature] = []
        for parent_id in parent_id2record.keys():
            parent_record = parent_id2record[parent_id]
            exons = parent_id2exons[parent_id]

            parent_feature = parent_record.to_seq_feature()
            if len(exons) == 1:
                exon_feature = parent_feature
            elif len(exons) >= 2:
                exons = sorted(exons, key=lambda e: e.start)
                locs = [e.to_feature_location() for e in exons]
                exon_feature = SeqFeature(
                    location=CompoundLocation(locs),
                    type=parent_feature.type,
                    strand=parent_feature.strand,
                    id=parent_feature.id,
                    qualifiers=parent_feature.qualifiers,
                )
            else:
                # If no exon exists, skip feature extraction
                continue

            exon_features.append(exon_feature)

        return exon_features

    def __str__(self):
        return f"{self.name} ({self.min_range:,} - {self.max_range:,} bp)"


@dataclass
class GffRecord:
    """GFF Record DataClass"""

    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: float | None
    strand: int
    phase: int | None
    attrs: dict[str, list[str]]

    def is_within_range(self, min_range: int, max_range: int) -> bool:
        """Check within target range or not

        Parameters
        ----------
        min_range : int
            Min range
        max_range : int
            Max range

        Returns
        -------
        check_result : bool
            Check result
        """
        if min_range <= self.start <= self.end <= max_range:
            return True
        else:
            return False

    def to_seq_feature(self) -> SeqFeature:
        """Convert GffRecord to SeqFeature"""
        return SeqFeature(
            location=self.to_feature_location(),
            type=self.type,
            strand=self.strand,
            id=self.attrs.get("ID", [""])[0],
            qualifiers=self.attrs,
        )

    def to_feature_location(self) -> FeatureLocation:
        """Convert GffRecord to FeatureLocation

        Returns
        -------
        feature_location : FeatureLocation
            Feature location
        """
        return FeatureLocation(self.start, self.end, self.strand)

    def to_gff_line(self) -> str:
        """Convert GffRecord to GFF record line

        Returns
        -------
        gff_line : str
            GFF record line (1-based coordinates)
        """
        return "\t".join(
            (
                self.seqid,
                self.source,
                self.type,
                str(self.start + 1),
                str(self.end),
                "." if self.score is None else str(self.score),
                "-" if self.strand == -1 else "+",
                "." if self.phase is None else str(self.phase),
                ";".join([f"{k}={','.join(v)}" for k, v in self.attrs.items()]),
            )
        )

    @staticmethod
    def is_gff_line(line: str) -> bool:
        """Check GFF record line or not

        Parameters
        ----------
        line : str
            GFF line

        Returns
        -------
        check_result : bool
            Check result
        """
        if line.startswith("#") or len(line.split("\t")) < 9:
            return False
        else:
            return True

    @staticmethod
    def parse_gff_line(gff_line: str) -> GffRecord:
        """Parse GFF record line

        Parameters
        ----------
        gff_line : str
            GFF record line (1-based coordinates)

        Returns
        -------
        gff_record : GffRecord
            GFF record (0-based coordinates)
        """
        gff_elms: list[Any] = gff_line.split("\t")[0:9]
        # start, end
        gff_elms[3], gff_elms[4] = int(gff_elms[3]) - 1, int(gff_elms[4])
        # score
        gff_elms[5] = None if gff_elms[5] in (".", "") else float(gff_elms[5])
        # strand
        if gff_elms[6] == "+":
            gff_elms[6] = 1
        elif gff_elms[6] == "-":
            gff_elms[6] = -1
        else:
            gff_elms[6] = 0
        # frame
        gff_elms[7] = None if gff_elms[7] in (".", "") else int(gff_elms[7])
        # attrs
        attr_dict: dict[str, list[str]] = {}
        attrs = str(gff_elms[8]).split(";")
        for attr in attrs:
            if attr != "" and "=" in attr:
                key, value = attr.split("=")
                attr_dict[key] = value.split(",")
        gff_elms[8] = attr_dict

        return GffRecord(*gff_elms)
