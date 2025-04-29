from __future__ import annotations

from pathlib import Path

import pandas as pd

from pycirclizely.utils import ColorCycler


class Table:
    """Table Parser Base Class"""

    def __init__(
        self,
        table_data: str | Path | pd.DataFrame,
        *,
        delimiter: str = "\t",
    ):
        """
        Parameters
        ----------
        table_data : str | Path | pd.DataFrame
            Table file or Table DataFrame
        delimiter : str, optional
            Table file delimiter. By default, `tab` delimiter.
        """
        if isinstance(table_data, (str, Path)):
            table_data = pd.read_csv(table_data, sep=delimiter, index_col=0)
        self._dataframe = table_data

    @property
    def dataframe(self) -> pd.DataFrame:
        """Table dataframe"""
        return self._dataframe

    @property
    def row_names(self) -> list[str]:
        """Table row names"""
        return list(map(str, self.dataframe.index))

    @property
    def col_names(self) -> list[str]:
        """Table column names"""
        return list(map(str, self.dataframe.columns))

    @property
    def row_num(self) -> int:
        """Table row count number"""
        return len(self.dataframe.index)

    @property
    def col_num(self) -> int:
        """Table column count number"""
        return len(self.dataframe.columns)

    def get_col_name2color(self, cmap: str = "tab10") -> dict[str, str]:
        """Get column name & color dict

        Parameters
        ----------
        cmap : str, optional
            Colormap (e.g. `tab10`, `Set3`)

        Returns
        -------
        col_name2color : dict[str, str]
            Column name & color dict
        """
        ColorCycler.set_cmap(cmap)
        return {n: ColorCycler.get_color() for n in self.col_names}

    def get_row_name2color(self, cmap: str = "tab10") -> dict[str, str]:
        """Get row name & color dict

        Parameters
        ----------
        cmap : str, optional
            Colormap (e.g. `tab10`, `Set3`)

        Returns
        -------
        col_name2color : dict[str, str]
            Column name & color dict
        """
        ColorCycler.set_cmap(cmap)
        return {n: ColorCycler.get_color() for n in self.row_names}

    def __str__(self):
        return str(self.dataframe)


class StackedBarTable(Table):
    """Table Parser Class

    Basically used for plotting stacked bar chart
    """

    @property
    def row_sum_vmax(self) -> float:
        """Max value in each row values sum"""
        return max(map(sum, self.dataframe.itertuples(index=False)))

    @property
    def row_name2sum(self) -> dict[str, float]:
        """Row name & sum dict"""
        row_sum_list = list(map(sum, self.dataframe.itertuples(index=False)))
        return dict(zip(self.row_names, row_sum_list))

    @property
    def stacked_bar_heights(self) -> list[list[float]]:
        """Stacked bar heights"""
        return [list(self.dataframe[col_name]) for col_name in self.col_names]

    @property
    def stacked_bar_bottoms(self) -> list[list[float]]:
        """Stacked bar bottoms"""
        bottoms: list[list[float]] = []
        row_name2stack_value = {name: 0.0 for name in self.row_names}
        for col_name in self.col_names:
            bottom = [row_name2stack_value[name] for name in self.row_names]
            for row_name in self.row_names:
                value = float(self.dataframe.at[row_name, col_name])
                row_name2stack_value[row_name] += value
            bottoms.append(bottom)
        return bottoms

    def calc_bar_label_x_list(
        self,
        track_size: float,
    ) -> list[float]:
        """Calculate list of x position for bar label plot

        Parameters
        ----------
        track_size : float
            Track size

        Returns
        -------
        x_list : list[float]
            List of x position for bar label plot
        """
        interval = track_size / len(self.row_names)
        return [cnt * interval + (interval / 2) for cnt in range(len(self.row_names))]

    def calc_barh_label_r_list(
        self,
        track_r_lim: tuple[float, float],
    ) -> list[float]:
        """Calculate list of radius position for horizontal bar label plot

        Parameters
        ----------
        track_r_lim : tuple[float, float]
            Track radius limit region

        Returns
        -------
        bar_label_r_list : list[float]
            List of radius position for horizontal bar label plot
        """
        rmin, rmax = track_r_lim
        interval = (rmax - rmin) / len(self.row_names)
        bar_label_r_list: list[float] = []
        for cnt in range(len(self.row_names)):
            r_center = rmax - (interval * cnt) - (interval / 2)
            bar_label_r_list.append(r_center)
        return bar_label_r_list

    def calc_barh_r_lim_list(
        self,
        track_r_lim: tuple[float, float],
        width: float = 0.8,
    ) -> list[tuple[float, float]]:
        """Calculate list of radius limit for horizontal bar plot

        Parameters
        ----------
        track_r_lim : tuple[float, float]
            Track radius limit region
        width : float, optional
            Bar width ratio (0.0 - 1.0)

        Returns
        -------
        list[tuple[float, float]]
            List of radius limit for horizontal bar plot
        """
        rmin, rmax = track_r_lim
        interval = (rmax - rmin) / len(self.row_names)
        bar_r_lim_list: list[tuple[float, float]] = []
        for cnt in range(len(self.row_names)):
            r_center = rmax - (interval * cnt) - (interval / 2)
            r_top = r_center + (interval / 2) * width
            r_bottom = r_center - (interval / 2) * width
            bar_r_lim_list.append((r_bottom, r_top))
        return bar_r_lim_list


class RadarTable(Table):
    """Radar Table Parser Class"""

    @property
    def row_name2values(self) -> dict[str, list[float]]:
        """Row name & values"""
        row_name2values = {}
        for row_name in self.row_names:
            row_name2values[row_name] = list(self.dataframe.loc[row_name])
        return row_name2values
