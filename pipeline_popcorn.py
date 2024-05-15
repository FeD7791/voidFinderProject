from vftkproject2.voidFinderProject.voidfindertk.sphericalvf import spherical_tools
from vftkproject2.voidFinderProject.voidfindertk.zobovvf import zobov_tools
from voidfindertk import io, data_box
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


# box = io.read_table('./datasets/halos_cut_scorpio.dat')
# #box = io.read_table('./datasets/halos_ascii_1000_1024_npmin_10_z0.51.dat', names = ['m','x','y','z','vx','vy','vz'])
# d_box = data_box.DataBox(box)
# model = spherical_tools.SphericalVF() 
# #model = zobov_tools.ZobovVF() 
# sparse = model.find(d_box)
# print(sparse)


# svf con datos de Dante
box = io.read_table(
    './datasets/halos_ascii_1000_1024_npmin_10_z0.51.dat', 
    names = ['m','x','y','z','vx','vy','vz'])
d_box = data_box.DataBox(box)
model = spherical_tools.SphericalVF()
output = model.preprocess(d_box, m_min = 200)
ou = model.model_find(output)