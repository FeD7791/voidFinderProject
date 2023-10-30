from astropy import units as u

import numpy as np

import pytest


from voidfindertk.box import Box


@pytest.fixture
def box_of_randoms(xyz_size=500, vxyz_size=220, m_size=12, n=1000):
    rng = np.random.default_rng(seed=42)
    x = xyz_size * rng.random(n)
    y = xyz_size * rng.random(n)
    z = xyz_size * rng.random(n)
    vx = vxyz_size * rng.random(n)
    vy = vxyz_size * rng.random(n)
    vz = vxyz_size * rng.random(n)
    m = m_size * rng.random(n)
    return Box(x, y, z, vx, vy, vz, m)


def test_box_properties(box_of_randoms):
    box = box_of_randoms

    units = [
        box.x.unit == u.Mpc,
        box.y.unit == u.Mpc,
        box.z.unit == u.Mpc,
        box.vx.unit == u.Mpc / u.h,
        box.vy.unit == u.Mpc / u.h,
        box.vz.unit == u.Mpc / u.h,
        box.m.unit == u.M_sun,
    ]
    assert all(units)


def test_box_inputlength():
    cols = []
    for e in np.arange(0, stop=7, step=1):
        cols.append(
            np.random.normal(0, 500, size=np.random.randint(250, size=1))
        )
    np.random.shuffle(cols)

    with pytest.raises(ValueError, match="Arrays should be of the same size"):
        Box(cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6])


def test_box_length(box_of_randoms):
    box = box_of_randoms
    assert (
        len(box) == list(
            {
                len(e)
                for e in [
                    box.x,
                    box.y,
                    box.z,
                    box.vx,
                    box.vy,
                    box.vz,
                    box.m,
                ]
            }
        )[0]
    )


def test_all_float(box_of_randoms):
    box = box_of_randoms
    arr_box = [box.x, box.y, box.z, box.vx, box.vy, box.vz, box.m]
    bool_list = []
    for arr in arr_box:
        bool_list.append(
            all([isinstance(e, u.quantity.Quantity) for e in arr])
        )
    assert all(bool_list)
