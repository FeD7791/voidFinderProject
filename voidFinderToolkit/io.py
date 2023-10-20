import pandas as pd
from .box import Box


def read_table(input_path, *args, **kwargs):
    """
    Builds a Box objects based on the specified columns of an input
    If *args or **kwargs specified takes the first 7 columns of the
    input

    Parameters
    ----------
    input_path : str
        path of the file
    *args : int
        List of integers that references to columns with the needed input
    **kwargs : int
        List of [x,y,z,vx,vy,vz,m] that references to columns with the 
        needed input
    """
    path = rf"{input_path}"
    data = pd.read_table(path, sep="\s+")
    data_cols = [0,1,2,3,4,5,6]

    if args and not kwargs:
        data_cols = list(args)[0:7]
    elif kwargs and not args:
        if  'x' in kwargs: data_cols[0] = kwargs['x'] 
        if  'y' in kwargs: data_cols[1] = kwargs['y']
        if  'z' in kwargs: data_cols[2] = kwargs['z']
        if  'vx' in kwargs: data_cols[3] = kwargs['vx']
        if  'vy' in kwargs: data_cols[4] = kwargs['vy']
        if  'vz' in kwargs: data_cols[5] = kwargs['vz']
        if  'm' in kwargs: data_cols[6] = kwargs['m']    
        
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

