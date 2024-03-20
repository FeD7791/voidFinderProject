from ..voidfindertk import io
from ..voidfindertk.spherical import spherical_interface as sphi
from ..voidfindertk import models
from ..voidfindertk import data_box
from ..voidfindertk.spherical_vf import SphericalVF
box = io.read_table('../halos_fede.dat')

box2 = io.read_tableZovob('../halos_fede.dat')

data_box = data_box.DataBox(box)
model = SphericalVF() #params =
sparse = model.find(data_box)
print(sparse)
#voids_output = sphi.spherical_void_finder(box)
