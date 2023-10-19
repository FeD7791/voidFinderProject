import pandas as pd
from box import Box


def box_builder(input_path, data_cols):
    """
    Builds a Box objects based on the specified columns of an input

    Parameters
    ----------
    input_path : str
        path of the file
    data_cols : list(int)
        List of integers that references to columns with the needed input
        data: [x,y,z,vx,vy,vz,m] (in that specific order) where:
        (x,y,z) : position coordinates in Mpc
        (vx,vy,vz) : velocity coordinates in Mpc/h
        m : mass value in solar masses
    """
    path = rf"{input_path}"
    data = pd.read_table(path, sep="\s+")
    box = Box(
        data.iloc[:, data_cols[0]],
        data.iloc[:, data_cols[1]],
        data.iloc[:, data_cols[2]],
        data.iloc[:, data_cols[3]],
        data.iloc[:, data_cols[4]],
        data.iloc[:, data_cols[5]],
        data.iloc[:, data_cols[6]],
    )
    return box
