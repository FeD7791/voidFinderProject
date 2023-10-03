from sys import path
import os 
funciones_path = os.path.realpath('../voidfinder_tk')
path.append(funciones_path)
from finders import spheric as spv

def test_voidfinder():
    test_void = spv.SphericalVoid(radius=np.nan,x=np.nan,y=np.nan,z=np.nan,vx=np.nan,vy=np.nan,vz=np.nan,delta=np.nan,d_type=np.nan)
    assert test_void.type_of_v() == 'No type'

def test_voidfinder():
    test_void = spv.SphericalVoid(17.328250,445.838043,168.413986,84.131096,-148.989853,-135.213730,111.886261,-0.901690,-0.025397)
    assert test_void.type_of_v() == 'R'