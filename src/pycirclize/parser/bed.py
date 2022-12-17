from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


class Bed:
    """BED Parser Class"""

    def __init__(self, bed_file: str | Path):
        """
        Parameters
        ----------
        bed_file : str | Path
            BED format file
        """
        self._bed_file = bed_file
        self._records = BedRecord.parse(bed_file)

    @property
    def records(self) -> list[BedRecord]:
        """BED records"""
        return self._records


@dataclass
class BedRecord:
    chr: str
    start: int
    end: int
    name: str | None = None
    score: str | None = None

    @property
    def size(self) -> int:
        """Size"""
        return self.end - self.start

    @staticmethod
    def parse(bed_file: str | Path) -> list[BedRecord]:
        """Parse BED format file

        Parameters
        ----------
        bed_file : str | Path
            BED format file

        Returns
        -------
        bed_records : list[BedRecord]
            BED records
        """
        bed_records = []
        with open(bed_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if row[0].startswith("#") or len(row) < 3:
                    continue
                try:
                    chr, start, end = row[0], int(row[1]), int(row[2])
                except Exception:
                    continue
                name, score = None, None
                if len(row) >= 5:
                    name, score = row[3], row[4]
                rec = BedRecord(chr, start, end, name, score)
                bed_records.append(rec)
        return bed_records
