import pytest

from voidfindertk import box, read_table


def test_read_table_type_output(buffer):
    assert isinstance(read_table(buffer), box.Box)


def test_nullvalues_on_read():
    with pytest.raises(TypeError):
        read_table("./tests/datasets/data_random_2.txt")


def test_number_of_columns(buffer):
    with pytest.raises(ValueError):
        read_table(buffer, usecols=[1, 2, 3])
