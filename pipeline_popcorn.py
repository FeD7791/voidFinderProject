from voidfindertk.io_popcorn import read_popcorn_output_to_dictionary
from voidfindertk.io_popcorn import create_sparse_matrix_from_dict
from voidfindertk import io
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

filename='./clean_popvds.dat'
datos = read_popcorn_output_to_dictionary(filename)

#Define params
n_items = 500
param = 'void_volume'

#Define df_voids
voids_list = list(datos['voids'])
voids_all = []
for i in voids_list:
	voids_all.append(datos['voids'][f'{i}'])
df_voids = pd.DataFrame.from_dict(voids_all)

#Define sparse
sparse = create_sparse_matrix_from_dict(datos,700070).transpose()
	
#Define box
df_box = pd.DataFrame(io.read_table('./halos_fede.dat').__dict__)

#OPeration
df_filtered_param = df_voids.sort_values(by=param, ascending=False)
df_filtered_index_cut = df_filtered_param.index[:n_items]
velocity_voids = []
for i in df_filtered_index_cut:
	tracer_indexes = list(sparse.getcol(i).tocoo().row)
	velocity_voids.append({
	'void':i, 
	'velocity_norm':np.linalg.norm(df_box.loc[tracer_indexes,['vx','vy','vz']].sum(axis=0))/len(tracer_indexes),
  	'n_tracers':len(tracer_indexes) 
	})

velocity_voids_df = pd.DataFrame(velocity_voids)
sns.histplot(data=velocity_voids_df, x='velocity_norm', bins=50, kde=True).set(title=f'nvoids:{len(velocity_voids)}')
#sns.histplot(data=velocity_voids_df, x='n_tracers', bins=50, kde=True).set(title=f'nvoids:{len(velocity_voids)}')

plt.savefig('popcorn_velocity_hist.jpg')
#plt.savefig('popcorn_tracers_hist.jpg')
plt.show()
