import pandas as pd
from .box import Box
from Tools.classes import MissingValuesError


def read_table(input_path,**kwargs):
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
    path = f"{input_path}"
    keys_set = {'x','y','z','vx','vy','vz','m'}

   
    if kwargs:
        try:
            if ('x' and 'y' and 'z'and 'vx' 
                and 'vy' and 'vz' and 'm' not in kwargs):
                missing_variables = keys_set.difference(set(kwargs.keys()))
                raise MissingValuesError('Missing the following values:',missing_variables)
            data_cols = [
                kwargs['x'],kwargs['y'],kwargs['z'],
                kwargs['vx'],kwargs['vy'],kwargs['vz'],
                kwargs['m'] ]
            data = pd.read_table(path, sep="\s+",usecols=data_cols)
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
        except MissingValuesError as error:
            print('MissingValuesError:')
            print(error)
    else:
        data_cols = [0,1,2,3,4,5,6]
        data = pd.read_table(path, sep="\s+",usecols=data_cols)
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
        

    

