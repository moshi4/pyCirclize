from pathlib import Path

import pandas as pd

from pycirclizely.parser import Matrix


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


def test_matrix_sort():
    """Test `matrix.sort()`"""
    matrix_df = pd.DataFrame(
        [
            [1, 2],
            [3, 4],
        ],
        index=["R1", "R2"],
        columns=["C1", "C2"],
    )
    matrix = Matrix(matrix_df)

    # Ascending Sort
    expected_asc_matrix_df = pd.DataFrame(
        [
            [0, 1, 2, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 3, 4, 0],
        ],
        index=["R1", "C1", "C2", "R2"],
        columns=["R1", "C1", "C2", "R2"],
    )
    asc_matrix_df = matrix.sort("asc").dataframe
    assert asc_matrix_df.equals(expected_asc_matrix_df)

    # Descending Sort
    expected_desc_matrix_df = pd.DataFrame(
        [
            [0, 4, 3, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 2, 1, 0],
        ],
        index=["R2", "C2", "C1", "R1"],
        columns=["R2", "C2", "C1", "R1"],
    )
    desc_matrix_df = matrix.sort("desc").dataframe
    assert desc_matrix_df.equals(expected_desc_matrix_df)

    # User-specified Order Sort
    expected_sort_matrix_df = pd.DataFrame(
        [
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [1, 2, 0, 0],
            [3, 4, 0, 0],
        ],
        index=["C1", "C2", "R1", "R2"],
        columns=["C1", "C2", "R1", "R2"],
    )
    sort_matrix_df = matrix.sort(["C1", "C2", "R1", "R2"]).dataframe
    assert sort_matrix_df.equals(expected_sort_matrix_df)


def test_to_fromto_table(tsv_matrix_file: Path):
    """Test `matrix.to_fromto_table()`"""
    matrix = Matrix(tsv_matrix_file)
    expected_table_df = pd.DataFrame(
        [
            ["S1", "E1", 4],
            ["S1", "E2", 14],
            ["S1", "E3", 13],
            ["S1", "E4", 17],
            ["S1", "E5", 5],
            ["S1", "E6", 2],
            ["S2", "E1", 7],
            ["S2", "E2", 1],
            ["S2", "E3", 6],
            ["S2", "E4", 8],
            ["S2", "E5", 12],
            ["S2", "E6", 15],
            ["S3", "E1", 9],
            ["S3", "E2", 10],
            ["S3", "E3", 3],
            ["S3", "E4", 16],
            ["S3", "E5", 11],
            ["S3", "E6", 18],
        ],
        columns=["from", "to", "value"],
    )
    assert matrix.to_fromto_table().equals(expected_table_df)
