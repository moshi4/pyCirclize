from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def testdata_dir() -> Path:
    """Testdata directory fixture"""
    return Path(__file__).parent / "testdata"


@pytest.fixture
def matrix_df() -> pd.DataFrame:
    """Pandas matrix dataframe"""
    matrix_data = [
        [4, 14, 13, 17, 5, 2],
        [7, 1, 6, 8, 12, 15],
        [9, 10, 3, 16, 11, 18],
    ]
    row_names = ["S1", "S2", "S3"]
    col_names = ["E1", "E2", "E3", "E4", "E5", "E6"]
    matrix_df = pd.DataFrame(matrix_data, index=row_names, columns=col_names)
    return matrix_df


@pytest.fixture
def fromto_table_df() -> pd.DataFrame:
    """Pandas from-to table dataframe"""
    return pd.DataFrame(
        data=[
            ["A", "B", 10],
            ["A", "C", 5],
            ["A", "D", 15],
            ["A", "E", 20],
            ["A", "F", 3],
            ["B", "A", 3],
            ["B", "G", 15],
            ["F", "D", 13],
            ["F", "E", 2],
            ["E", "A", 20],
            ["E", "D", 6],
        ],
    )


@pytest.fixture
def csv_matrix_file(matrix_df: pd.DataFrame, tmp_path: Path) -> Path:
    """CSV matrix file fixture"""
    csv_matrix_file = tmp_path / "matrix.csv"
    matrix_df.to_csv(csv_matrix_file)
    return csv_matrix_file


@pytest.fixture
def tsv_matrix_file(matrix_df: pd.DataFrame, tmp_path: Path) -> Path:
    """TSV matrix file fixture"""
    tsv_matrix_file = tmp_path / "matrix.tsv"
    matrix_df.to_csv(tsv_matrix_file, sep="\t")
    return tsv_matrix_file
