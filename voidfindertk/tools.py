import numpy as np
import ctypes
import pandas as pd
from . import box, data_box

def ctypes_input_params_builder(finder_type,**kwargs):
    if finder_type == 'spherical':
        # Input params
        params = {
        #Double Values
        "RadIncrement" : 0. ,
        "BoxSize": 1000 ,
        "MaxRadiusSearch": 40.0,
        "ProxyGridSize": 5.0,
        "FracRadius" : 0.5 ,
        "DeltaThreshold": -0.9,  
        "DeltaSeed": -0.7,
        "OverlapTol":0,
        "Redshift" : 0.99 ,
        "OmegaMatter" : 0.25,
        "OmegaLambda" : 0.75,
        "Hubble" : 0.73,
        "FidOmegaMatter": 0.2  ,
        "FidOmegaLambda" : 0.8 ,
        "FidHubble" : 0.7, 
        "MinProfileDist": 0.5,
        "MaxProfileDist" : 3.0,
        "ScalePos": 1,
        "ScaleVel" : 1,
        "InnerShellVel": 0.8,
        "OuterShellVel": 1.2,
        #Int Values
        "FormatTracers": 0,
        "NumFiles": 32,
        "NumRanWalk" : 75 , 
        "OMPcores": 8,
        "RSDist" : 0,
        "GDist": 0,
        "WriteProfiles" : 0 ,
        "NumProfileBins" : 100,   
            
        }
        #Replace params as needed
        for key,value in kwargs.items():
            params[key] = value
        
        # Load Input Params
        params1 = dict(list(params.items())[0:21]) # Params that are double type
        params2 = dict(list(params.items())[21:29]) # Params that are int type

        # Build ctypes structure for params
        p1 = [(key,ctypes.c_double) for key,value in params1.items()] 
        p2 = [(key,ctypes.c_int) for key,value in params2.items()]

        class InputParams(ctypes.Structure):
            _fields_ = p1 + p2
    
    return {'InputParams_class':InputParams,'params_dict':params}

def ctypes_output_builder(finder_type):
    if finder_type == 'spherical':
        # Create stucts --> class for output
        class voids(ctypes.Structure):
            _fields_ = [
                ("n_voids", ctypes.c_int),
                ("Rad", ctypes.c_float),
                # ("Rini", ctypes.c_float),
                # ("Ini", ctypes.c_float * 1),
                ("Pos", ctypes.c_float * 3),
                ("Vel", ctypes.c_float * 3),
                ("Dtype", ctypes.c_float),
                ("Delta", ctypes.c_float),
                ("Poisson", ctypes.c_float),
                # ("Dist4", ctypes.c_float),
                # ("ToF", ctypes.c_bool),
                ("Nran", ctypes.c_int),
            ]

        class voidArray(ctypes.Structure):
            _fields_ = [("voids", voids * 1)]
    return {'voidArray_class':voidArray}

def process_output_from_finder(finder_type, array_of_voids):
    va = array_of_voids
    if finder_type == 'spherical':
        va_arr = np.ctypeslib.as_array(va, shape=(10,))

        n_voids = va_arr[0][0][0][0] - 1

        full_arr = np.ctypeslib.as_array(va, shape=(n_voids,))

        radius = []
        x_coord = []
        y_coord = []
        z_coord = []
        vx_coord = []
        vy_coord = []
        vz_coord = []
        dtype = []
        delta = []
        poisson = []
        nran = []
        for i in range(0, n_voids, 1):
            radius.append(full_arr[i][0][0][1])
            x_coord.append(full_arr[i][0][0][2][0])
            y_coord.append(full_arr[i][0][0][2][1])
            z_coord.append(full_arr[i][0][0][2][2])
            vx_coord.append(full_arr[i][0][0][3][0])
            vy_coord.append(full_arr[i][0][0][3][1])
            vz_coord.append(full_arr[i][0][0][3][2])
            dtype.append(full_arr[i][0][0][4])
            delta.append(full_arr[i][0][0][5])
            poisson.append(full_arr[i][0][0][6])
            nran.append(full_arr[i][0][0][7])
    
    return {
        'rad':radius, 
        'x_void':x_coord, 'y_void':y_coord, 'z_void':z_coord,
        'vel_x_void':vx_coord, 'vel_y_void':vy_coord, 'vel_z_void':vz_coord,
        'delta':delta, 'dtype':dtype, 'poisson':poisson, 'nran':nran
        }

def preprocess_data_box(databox,**kwargs):
    kwargs.setdefault('m_min', 0)
    kwargs.setdefault('m_max', len(databox.box))
    b = databox.box
    ##
    df = pd.DataFrame(b.__dict__)
    df = df[df['m'] >= kwargs['m_min']]
    df = df[df['m'] <= kwargs['m_max']]
    df.reset_index(drop=True,inplace=True)
    df.drop(columns=['_len'],axis=1,inplace=True)
    df.drop_duplicates(inplace=True, ignore_index=True) #Drop any remaining duplicates
    box2 = box.Box(**df.to_dict(orient='list'))
    return data_box.DataBox(box2)
