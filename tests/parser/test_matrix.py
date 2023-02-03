from pathlib import Path

import pandas as pd

from pycirclize.parser import Matrix


def test_load_dataframe_matrix(matrix_df: pd.DataFrame):
    """Test load panda dataframe matrix"""
    # Load pandas matrix dataframe
    matrix = Matrix(matrix_df)

    # Test row & column names
    row_names = ["S1", "S2", "S3"]
    col_names = ["E1", "E2", "E3", "E4", "E5", "E6"]
    assert matrix.all_names == row_names + col_names
    assert matrix.row_names == row_names
    assert matrix.col_names == col_names

    # Only test successfully call function
    matrix.to_sectors()
    matrix.to_links()


def test_parse_fromto_table(fromto_table_df: pd.DataFrame):
    """Test parse from-to table"""
    # Parse from-to table dataframe
    matrix = Matrix.parse_fromto_table(fromto_table_df)

    # Test row & column names
    expected_names = list("ABCDEFG")
    assert matrix.all_names == expected_names
    assert matrix.row_names == expected_names
    assert matrix.col_names == expected_names

    # Only test successfully call function
    matrix.to_sectors()
    matrix.to_links()


def test_load_tsv_matrix(tsv_matrix_file: Path):
    """Test load tsv matrix"""
    # Load tsv format matrix file
    matrix = Matrix(tsv_matrix_file)

    # Test row & column names
    row_names = ["S1", "S2", "S3"]
    col_names = ["E1", "E2", "E3", "E4", "E5", "E6"]
    assert matrix.all_names == row_names + col_names
    assert matrix.row_names == row_names
    assert matrix.col_names == col_names

    # Only test successfully call function
    matrix.to_sectors()
    matrix.to_links()


def test_load_csv_matrix(csv_matrix_file: Path):
    """Test load csv matrix"""
    # Load csv format matrix file
    matrix = Matrix(csv_matrix_file, delimiter=",")

    # Test row & column names
    row_names = ["S1", "S2", "S3"]
    col_names = ["E1", "E2", "E3", "E4", "E5", "E6"]
    assert matrix.all_names == row_names + col_names
    assert matrix.row_names == row_names
    assert matrix.col_names == col_names

    # Only test successfully call function
    matrix.to_sectors()
    matrix.to_links()
