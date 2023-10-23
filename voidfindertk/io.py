import pandas as pd
from .box import Box
from tools.classes import MissingValuesError


def read_table(input_path,**kwargs):
    """
    Builds a Box objects based on the specified columns of an input
    If *args or **kwargs specified takes the first 7 columns of the
    input

    Parameters
    ----------
    input_path : str
        path of the file
    **kwargs : int
        List of [x,y,z,vx,vy,vz,m] that references to columns with the 
        needed input
    """
    path = input_path
    keys_set = {'x','y','z','vx','vy','vz','m'}

   
    if kwargs:
        if (not (all(key in kwargs for key in keys_set))): 
            
            missing_variables = keys_set.difference(set(kwargs.keys()))
            raise MissingValuesError('Missing the following values:',missing_variables) 
        
        data_cols = [
        kwargs['x'],kwargs['y'],kwargs['z'],
        kwargs['vx'],kwargs['vy'],kwargs['vz'],
        kwargs['m'] ]
        
    else:
        data_cols = [0,1,2,3,4,5,6]
        
    data = pd.read_table(path, sep="\s+",usecols=data_cols) 
    box = Box(
        data.iloc[:,0],
        data.iloc[:,1],
        data.iloc[:,2],
        data.iloc[:,3],
        data.iloc[:,4],
        data.iloc[:,5],
        data.iloc[:,6]
        )
    return box
        

    

