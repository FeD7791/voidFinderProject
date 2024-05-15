import numpy as np
from voidfindertk import box as b
from voidfindertk.zobovvf.zobov_tools import ZobovVoids, calculate_tracers_inside_void
from unittest.mock import MagicMock

def test_calculate_tracers_inside_void():
    a = np.array([0,1,2,3,4] + list(np.arange(5,1000)))
    box_params = ['x','y','z','vx','vy','vz','m']
    voids_params = [
        'Void_number', 
        'File_void_number', 
        'CoreParticle', 
        'CoreDens', 
        'ZoneVol', 
        'Zone_number_part', 
        'Void_number_Zones', 
        'VoidVol', 
        'Void_number_Part', 
        'VoidDensContrast', 
        'VoidProb']
    box = b.Box(**dict([(i,a) for i in box_params]))