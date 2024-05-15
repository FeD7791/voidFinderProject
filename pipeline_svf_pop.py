from vftkproject2.voidFinderProject.voidfindertk.popcorn import popcorn_tools
from voidfindertk import io, analysis_tools
import pandas as pd

path = './datasets/sphvds.dat'
svf_voids = pd.read_csv(path,sep='\s+')
svf_voids.columns = ['index','rad','x_void','y_void','z_void','delta']
svf = popcorn_tools.SVF(**svf_voids[['rad','x_void','y_void','z_void']].to_dict(orient='list'))

#box
box = io.read_table('./datasets/halos_ascii_1000_1024_npmin_10_z0.51.dat',names=['m','x','y','z','vx','vy','vz'])

#SPARSE COMBO
output = analysis_tools.join_box_void(box=box,voids=svf)

#VSF
vsf = analysis_tools.void_size_function_2(box=box,voids=svf)