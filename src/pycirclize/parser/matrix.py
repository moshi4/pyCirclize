from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd


class Matrix:
    """Matrix Parser Class"""

    def __init__(self, matrix: str | Path | pd.DataFrame, delimiter: str = "\t"):
        """
        Parameters
        ----------
        matrix : str | Path | pd.DataFrame
            Matrix file or Matrix DataFrame
        delimiter : str, optional
            Matrix file delimiter. By default, `tab` delimiter.
        """
        # If input matrix is file path, convert to pandas dataframe
        if isinstance(matrix, (str, Path)):
            matrix = pd.read_csv(matrix, delimiter=delimiter, index_col=0)

        # Calculate data size & link positions
        rev_matrix = matrix.iloc[::-1, ::-1]
        name2size, links = defaultdict(float), []
        for row_name, row in zip(rev_matrix.index, rev_matrix.values):
            for col_name, value in zip(rev_matrix.columns, row):
                row_size, col_size = name2size[row_name], name2size[col_name]
                if row_name == col_name:
                    link_row = (row_name, row_size, row_size + value)
                    link_col = (col_name, col_size + (value * 2), col_size + value)
                else:
                    link_row = (row_name, row_size, row_size + value)
                    link_col = (col_name, col_size + value, col_size)
                links.append((link_row, link_col))
                name2size[row_name] += value
                name2size[col_name] += value

        self._matrix = matrix
        self._col_names = list(matrix.columns)
        self._row_names = list(matrix.index)
        self._links = links
        self._name2size = name2size

    @property
    def all_names(self) -> list[str]:
        """Row + Column all names"""
        return list(self.to_sectors().keys())

    @property
    def col_names(self) -> list[str]:
        """Column names"""
        return self._col_names

    @property
    def row_names(self) -> list[str]:
        """Row names"""
        return self._row_names

    def to_sectors(self) -> dict[str, float]:
        """Convert matrix to sectors for Circos initialization

        >>> # Example usage
        >>> matrix = Matrix(matrix_file)
        >>> circos = Circos(matrix.to_sectors())

        Returns
        -------
        sectors : dict[str, float]
            Sector dict (e.g. `{"A": 12, "B": 15, "C":20, ...}`)
        """
        sectors = {}
        for row_name in self.row_names:
            sectors[row_name] = self._name2size[row_name]
        for col_name in self.col_names:
            sectors[col_name] = self._name2size[col_name]
        return sectors

    def to_links(
        self,
    ) -> list[tuple[tuple[str, float, float], tuple[str, float, float]]]:
        """Convert matrix to links data for `circos.link()` method

        >>> # Example usage
        >>> matrix = Matrix(matrix_file)
        >>> circos = Circos(matrix.to_sectors())
        >>> for link in matrix.to_links():
        >>>    circos.link(*link)

        Returns
        -------
        links : list[tuple[tuple[str, float, float], tuple[str, float, float]]]
            List of link `((name1, start1, end1), (name2, end2, start2))`
        """
        return self._links

    def __str__(self):
        return str(self._matrix)
