from voidfindertk import spherical_tools
from astropy import units as u
import pytest
from voidfindertk.spherical_tools import SphericalVoids
import numpy as np
import pandas as pd

def test_spherical_voids_output_units(make_spherical_voids_params):
    params = make_spherical_voids_params(n_voids=500)
    spherical_voids = SphericalVoids(**params)
    assert spherical_voids.rad.unit == u.Mpc
    assert spherical_voids.x_void.unit == u.Mpc
    assert spherical_voids.y_void.unit == u.Mpc
    assert spherical_voids.z_void.unit == u.Mpc
    assert spherical_voids.vel_x_void.unit == u.Mpc/u.h
    assert spherical_voids.vel_y_void.unit == u.Mpc/u.h
    assert spherical_voids.vel_z_void.unit == u.Mpc/u.h
    assert len(spherical_voids) == 500
    assert str(spherical_voids) == "<SphericalVoids size=500>"
    #wont test slice method because i dont know if i will keep it

