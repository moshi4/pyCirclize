from __future__ import annotations

import bz2
import gzip
import warnings
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from io import TextIOWrapper
from pathlib import Path
from typing import Any, TextIO

from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation


class Gff:
    """GFF Parser Class"""

    def __init__(
        self,
        gff_file: str | Path,
        *,
        name: str | None = None,
        target_seqid: str | None = None,
        min_range: None = None,
        max_range: None = None,
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
        min_range : None, optional
            No longer used. Left for backward compatibility.
        max_range : None, optional
            No longer used. Left for backward compatibility.
        """
        self._gff_file = Path(gff_file)
        self._name = name
        self._records, start, end = self._parse_gff(gff_file, target_seqid)
        self._seq_region = (start, end)

        if min_range or max_range:
            warnings.warn("min_range & max_range is no longer used in Gff parser.")

    ############################################################
    # Property
    ############################################################

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
        """GFF records (only target seqid)"""
        return self._records

    @property
    def all_records(self) -> list[GffRecord]:
        """All GFF records"""
        return self._all_records

    @property
    def target_seqid(self) -> str:
        """Target seqid"""
        return self._target_seqid

    @property
    def seqid_list(self) -> list[str]:
        """seqid list"""
        return self._seqid_list

    @property
    def genome_length(self) -> int:
        """Genome length (target seqid record)"""
        return max(self.seq_region)

    @property
    def range_size(self) -> int:
        """Same as `self.genome_length` (Left for backward compatibility)"""
        return self.genome_length

    @property
    def full_genome_length(self) -> int:
        """Full genome length (concatenate all records)"""
        return sum(list(self.get_seqid2size().values()))

    ############################################################
    # Public Method
    ############################################################

    def get_seqid2size(self) -> dict[str, int]:
        """Get seqid & complete/contig/scaffold genome size dict

        By default, size is defined by `##sequence-region` pragma of target seqid.
        If `##sequence-region` is not found, size is defined by max coordinate size in
        target seqid features. This may differ from actual genome size.

        Returns
        -------
        seqid2size : dict[str, int]
            seqid & genome size dict
        """
        return self._seqid2size

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

        gff_records = GffRecord.filter_records(
            self.all_records,
            feature_type=feature_type,
            target_strand=target_strand,
        )
        seqid2features: dict[str, list[SeqFeature]] = {}
        for seqid in self.seqid_list:
            seqid2features[seqid] = []
        for rec in gff_records:
            seqid2features[rec.seqid].append(rec.to_seq_feature())
        return seqid2features

    def extract_features(
        self,
        feature_type: str | list[str] | None = "CDS",
        *,
        target_strand: int | None = None,
        target_range: tuple[int, int] | None = None,
    ) -> list[SeqFeature]:
        """Extract features

        If `target_seqid` is specified when the Gff instance initialized,
        then the features of the target seqid are extracted.
        Otherwise, extract the features of the seqid in the first row.

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
            Feature list
        """
        gff_records = GffRecord.filter_records(
            self.records,
            feature_type=feature_type,
            target_strand=target_strand,
            target_range=target_range,
        )
        return [rec.to_seq_feature() for rec in gff_records]

    def extract_exon_features(self, feature_type: str = "mRNA") -> list[SeqFeature]:
        """Extract exon structure features

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
        for rec in self.records:
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
                    id=parent_feature.id,
                    qualifiers=parent_feature.qualifiers,
                )
            else:
                # If no exon exists, skip feature extraction
                continue

            exon_features.append(exon_feature)

        return exon_features

    ############################################################
    # Private Method
    ############################################################

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
        # Parse GFF lines
        gff_all_lines = handle.read().splitlines()
        gff_record_lines = filter(GffRecord.is_gff_line, gff_all_lines)
        gff_records = list(map(GffRecord.parse_gff_line, gff_record_lines))
        if len(gff_records) == 0:
            err_msg = f"Failed to parse '{self._gff_file}' as GFF file "
            raise ValueError(err_msg)

        # Get target seqid & GFF records
        seqid_list = list(dict.fromkeys([rec.seqid for rec in gff_records]))
        if target_seqid is None:
            target_seqid = seqid_list[0]
        if target_seqid not in seqid_list:
            err_msg = f"Not found {target_seqid=} in '{self._gff_file}'"
            raise ValueError(err_msg)
        target_gff_records = [rec for rec in gff_records if rec.seqid == target_seqid]

        # Try to get start-end region from '##sequence-region' annotation line
        # If not found, (0, max_coordinate) is set as start-end
        seqid2start_end: dict[str, tuple[int, int]] = {}
        for seqid in seqid_list:
            start, end = None, None
            for line in gff_all_lines:
                if line.startswith("##sequence-region"):
                    # e.g. `##sequence-region NC_XXXXXX 1 10000` (seqid, start, end)
                    if len(line.split()) == 4 and line.split()[1] == seqid:
                        start, end = line.split()[2:4]
                        start, end = int(start) - 1, int(end)
                        break
            if start is None or end is None:
                seqid_gff_records = [rec for rec in gff_records if rec.seqid == seqid]
                start, end = 0, max([r.end for r in seqid_gff_records])
            seqid2start_end[seqid] = (start, end)
        seqid2size = {seqid: e - s for seqid, (s, e) in seqid2start_end.items()}

        # Set properties
        self._target_seqid = target_seqid
        self._seqid_list = seqid_list
        self._all_records = gff_records
        self._seqid2size = seqid2size

        return target_gff_records, *seqid2start_end[target_seqid]


@dataclass
class GffRecord:
    """GFF Record DataClass"""

    seqid: str
    source: str
    type: str
    start: int  # 1-based coordinate
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
        """Convert GffRecord to SeqFeature (1-based to 0-based coordinate)"""
        return SeqFeature(
            location=self.to_feature_location(),
            type=self.type,
            id=self.attrs.get("ID", [""])[0],
            qualifiers=self.attrs,
        )

    def to_feature_location(self) -> SimpleLocation:
        """Convert GffRecord to SimpleLocation (1-based to 0-based coordinate)

        Returns
        -------
        feature_location : SimpleLocation
            Simple location
        """
        return SimpleLocation(self.start - 1, self.end, self.strand)

    def to_gff_line(self) -> str:
        """Convert GffRecord to GFF record line

        Returns
        -------
        gff_line : str
            GFF record line
        """
        return "\t".join(
            (
                self.seqid,
                self.source,
                self.type,
                str(self.start),
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
            GFF record line

        Returns
        -------
        gff_record : GffRecord
            GFF record
        """
        gff_elms: list[Any] = gff_line.split("\t")[0:9]
        # start, end
        gff_elms[3], gff_elms[4] = int(gff_elms[3]), int(gff_elms[4])
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

    @staticmethod
    def filter_records(
        gff_records: list[GffRecord],
        feature_type: str | list[str] | None = "CDS",
        target_strand: int | None = None,
        target_range: tuple[int, int] | None = None,
    ) -> list[GffRecord]:
        """Filter GFF records by feature_type, strand, range

        Parameters
        ----------
        gff_records : list[GffRecord]
            GFF records to be filterd
        feature_type : str | list[str] | None, optional
            Feature type (`CDS`, `gene`, `mRNA`, etc...). If None, no filter.
        target_strand : int | None, optional
            Target strand. If None, no filter.
        target_range : tuple[int, int] | None, optional
            Target range. If None, no filter.

        Returns
        -------
        filter_gff_records : list[SeqFeature]
            Filtered GFF records
        """
        if isinstance(feature_type, str):
            feature_type = [feature_type]

        filter_gff_records = []
        for rec in gff_records:
            if feature_type is not None and rec.type not in feature_type:
                continue
            if target_strand is not None and rec.strand != target_strand:
                continue
            if target_range is not None:
                min_range, max_range = min(target_range), max(target_range)
                if not min_range <= rec.start <= rec.end <= max_range:
                    continue
            filter_gff_records.append(rec)
        return filter_gff_records
