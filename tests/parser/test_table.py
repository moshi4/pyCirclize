from __future__ import annotations

import pandas as pd
import pytest

from pycirclizely.parser import RadarTable, StackedBarTable


class TestStackedBarTable:
    """Test StackedBarTable Class"""

    @pytest.fixture
    def sb_table(self):
        """Initialize stacked bar table fixture"""
        table_df = pd.DataFrame(
            [
                [1, 4, 7, 10],
                [2, 5, 8, 11],
                [3, 6, 9, 12],
            ],
            index=["row1", "row2", "row3"],
            columns=["col1", "col2", "col3", "col4"],
        )
        return StackedBarTable(table_df)

    def test_properties(self, sb_table: StackedBarTable):
        """Test properties"""
        assert sb_table.row_names == ["row1", "row2", "row3"]
        assert sb_table.col_names == ["col1", "col2", "col3", "col4"]
        assert sb_table.row_num == 3
        assert sb_table.col_num == 4
        assert sb_table.row_sum_vmax == 3 + 6 + 9 + 12
        assert sb_table.row_name2sum == dict(row1=22, row2=26, row3=30)
        assert sb_table.stacked_bar_heights == [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
            [10, 11, 12],
        ]
        assert sb_table.stacked_bar_bottoms == [
            [0, 0, 0],
            [1, 2, 3],
            [5, 7, 9],
            [12, 15, 18],
        ]

    @pytest.mark.parametrize(
        "track_size, expected_x_list",
        [
            (30, [5, 15, 25]),
            (150, [25, 75, 125]),
        ],
    )
    def test_calc_bar_label_x_list(
        self,
        sb_table: StackedBarTable,
        track_size: float,
        expected_x_list: list[float],
    ):
        """Test `calc_bar_label_x_list()`"""
        x_list = sb_table.calc_bar_label_x_list(track_size)
        assert x_list == expected_x_list

    @pytest.mark.parametrize(
        "track_r_lim, expected_r_list",
        [
            ((70, 100), [95, 85, 75]),
            ((10, 70), [60, 40, 20]),
        ],
    )
    def test_calc_barh_label_r_list(
        self,
        sb_table: StackedBarTable,
        track_r_lim: tuple[float, float],
        expected_r_list: list[float],
    ):
        """Test `calc_barh_label_r_list()`"""
        r_list = sb_table.calc_barh_label_r_list(track_r_lim)
        assert r_list == expected_r_list

    @pytest.mark.parametrize(
        "track_r_lim, width, expected_r_lim_list",
        [
            ((70, 100), 1.0, [(90, 100), (80, 90), (70, 80)]),
            ((10, 70), 0.8, [(52, 68), (32, 48), (12, 28)]),
        ],
    )
    def test_calc_barh_r_lim_list(
        self,
        sb_table: StackedBarTable,
        track_r_lim: tuple[float, float],
        width: float,
        expected_r_lim_list: list[tuple[float, float]],
    ):
        """Test `calc_barh_r_lim_list()`"""
        r_lim_list = sb_table.calc_barh_r_lim_list(track_r_lim, width)
        assert r_lim_list == expected_r_lim_list


class TestRaderTable:
    """Test RadarTable Class"""

    @pytest.fixture
    def radar_table(self):
        """Initialize radar table fixture"""
        table_df = pd.DataFrame(
            data=[
                [80, 80, 80, 80, 80, 80],
                [90, 95, 95, 30, 30, 80],
                [60, 20, 20, 100, 90, 50],
            ],
            index=["Hero", "Warrior", "Wizard"],
            columns=["HP", "ATK", "DEF", "SP.ATK", "SP.DEF", "SPD"],
        )
        return RadarTable(table_df)

    def test_row_name2values(self, radar_table: RadarTable):
        """Test `row_name2values()`"""
        assert radar_table.row_name2values == dict(
            Hero=[80, 80, 80, 80, 80, 80],
            Warrior=[90, 95, 95, 30, 30, 80],
            Wizard=[60, 20, 20, 100, 90, 50],
        )
