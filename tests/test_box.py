import pytest
import numpy as np
from astropy import units as u
from voidfindertk.box import Box

@pytest.fixture
def box_of_randoms(xyz_size=500,vxyz_size=220,m_size=12,n=1000):
    x = np.random.normal(0,xyz_size,size=n)
    y = np.random.normal(0,xyz_size,size=n)
    z = np.random.normal(0,xyz_size,size=n)
    vx = np.random.normal(0,vxyz_size,size=n)
    vy = np.random.normal(0,vxyz_size,size=n)
    vz = np.random.normal(0,vxyz_size,size=n)
    m = np.random.normal(0,m_size,size=n)
    return Box(x,y,z,vx,vy,vz,m)
    





## Verifica que un objeto sea clase box
def test_box(box_of_randoms):
    assert isinstance(box_of_randoms,Box)
     

# Nuestro objeto debe fallar cuando le entregamos un vector de distinto rango
# En este asumo que los inputs no tienen el mismo tamanio
def test_box_inputlength():
    cols = []
    for e in np.arange(0,stop=7,step=1):
        cols.append(np.random.normal(0,500,size=np.random.randint(250,size=1)))
    np.random.shuffle(cols)

    with pytest.raises(ValueError,match="Arrays should be of the same size"):
        Box(cols[0],cols[1],cols[2],cols[3],cols[4],cols[5],cols[6])

# Si todos los elementos del input del box son same size entonces su longitud deberia ser la misma que la caja
def test_box_length(box_of_randoms):
    box = box_of_randoms
    assert len(box) == list({len(e) for e in [box.x,box.y,box.z,box.vx,box.vy,box.vz,box.m]})[0]

#Pureba que todos las propiedades en box sean unidades de astropy
def test_all_float(box_of_randoms):
    box = box_of_randoms
    arr_box = [box.x,box.y,box.z,box.vx,box.vy,box.vz,box.m]
    bool_list = []
    for arr in arr_box:
        bool_list.append(all([isinstance(e,u.quantity.Quantity) for e in arr]))
    assert all(bool_list)