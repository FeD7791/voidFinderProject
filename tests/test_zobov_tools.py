import numpy as np
from voidfindertk import box as b
from voidfindertk.zobovvf.zobov_tools import ZobovVoids, calculate_tracers_inside_void

def test_calculate_tracers_inside_void(mkbox, make_zobov_voids_params):
    