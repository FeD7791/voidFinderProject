import pandas as pd
import numpy as np

def analysis_voids(void_out,n_items,param): #param: a param like 'rad', 'volume' property of a void

	df_voids = pd.DataFrame(void_out._get_voids())
	df_box = pd.DataFrame(void_out._get_box())
	sparse = void_out._sparse_matrix()
	df_filtered_param = df_voids.sort_values(by=param, ascending=False)
	df_filtered_index_cut = df_filtered_param.index[:n_items]
	parameter_voids = []
	for i in df_filtered_index_cut:
		tracer_indexes = list(sparse.getcol(i).tocoo().row)
		parameter_voids.append({
		'void':i, 
		'velocity_norm':np.linalg.norm(df_box.loc[tracer_indexes,['vx','vy','vz']].sum(axis=0))/len(tracer_indexes),
    'n_tracers':len(tracer_indexes)
		})
	return parameter_voids
