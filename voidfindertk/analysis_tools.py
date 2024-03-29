import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
from scipy.sparse import csr_matrix


def join_box_void(box,voids, **kwargs):
		
		params = {'tol':0.0}
		for key,value in kwargs.items():
			params[key] = value

		if type(voids).__name__ == 'SphericalVoids':

			tolerance = params['tol'] * np.ones(voids.rad.shape)
			rad = np.array(voids.rad) + tolerance
			pos_void = np.array([voids.x_void, voids.y_void, voids.z_void])
			pos_void = np.transpose(pos_void)
			arr = []
			for i in range(box._len):
				pos_box = np.array(
					[[box.x[i].value, box.y[i].value, box.z[i].value]]
				)
				d = scipy.spatial.distance.cdist(
					pos_box, pos_void, metric="euclidean"
				)  # Calculates the distance from a tracer to each void center
				rad_center_distance_comparing = np.greater(rad, d[0]).astype(
					float
				)  # Compares radius of a void and distance from the void center to the particle 1 if rad>d # noqa: E501
				arr.append(rad_center_distance_comparing)
			# Transform to sparse matrix
			output_sparse_matrix = scipy.sparse.csr_matrix(arr)

		if type(voids).__name__ == 'PopCornVoids':
			n_index =[]
			tracers = []

			for i in range(0,voids._void_len):
				n_index = n_index + list(i*np.ones(np.array(voids.ID_subhalo[i]).shape))
				tracers = tracers + voids.ID_subhalo[i]


			row = np.array(tracers)
			col = np.array(n_index)
			data = np.ones(row.shape)
			output_sparse_matrix = csr_matrix((data, (row, col)), shape=(box._len, voids._void_len))

		return {'box':box, 'voids':voids, 'sparse':output_sparse_matrix}



def analysis_voids(box,voids,sparse,n_items): #param: a param like 'rad', 'volume' property of a void
											#n_items: the most relevant n_items

	df_voids = pd.DataFrame(voids.__dict__)
	df_box = pd.DataFrame(box.__dict__)
	
	if type(voids).__name__ == 'SphericalVoids':
		df_filtered_param = df_voids.sort_values(by='rad', ascending=False)
	if type(voids).__name__  == 'PopCornVoids':
		df_filtered_param = df_voids.sort_values(by='void_volume', ascending=False)

	df_filtered_index_cut = df_filtered_param.index[:n_items]

	parameter_voids = {'void_number':[], 'velocity_norm':[], 'n_tracers':[]}

	for i in df_filtered_index_cut:
		tracer_indexes = list(sparse.getcol(i).tocoo().row) 
		parameter_voids['void_number'].append(i)
		parameter_voids['velocity_norm'].append(np.linalg.norm(df_box.loc[tracer_indexes,['vx','vy','vz']].sum(axis=0))/len(tracer_indexes))
		parameter_voids['n_tracers'].append(len(tracer_indexes))
	return parameter_voids


def void_size_function(box,voids, **kwargs):

	params = {
		'n_bins':10, 
		'save_fig':False,
		'plot':False
		}
	for key,value in kwargs.items():
		params[key] = value

	#Get radii
	if type(voids).__name__ == 'SphericalVoids':
		rad = voids.rad.value
	if type(voids).__name__ == 'PopCornVoids':
		popcorn_volume = voids.void_volume
		rad = ((3/(4*np.pi))*popcorn_volume)**(1/3)
	
	void_radii = np.array(rad)
	bins = np.linspace(min(void_radii), max(void_radii), params['n_bins'] )
	# Count voids in each bin
	counts, bins0 = np.histogram(void_radii, bins=bins)


	mids = (bins[1:] + bins[:-1]) / 2
	bin_width = bins[1:]- bins[:-1] 

	#Calculate survey volume 
	survey_volume = (math.ceil(np.max(box.x.value)))**3 

	# Normalize Counts
	normalized_counts = counts / (survey_volume * np.log10(bin_width))
	
	

	# Plot the void size function
	if params['plot']:
		plt.figure(figsize=(8, 8))
		plt.plot(mids, normalized_counts, marker='o', linestyle='-', label='Void Size Function')
		plt.xlabel(r'$log_{10}(R)$ '+ f' {box.x[0].unit}')
		plt.ylabel(r'$\frac{1}{V} \frac{dN_v}{dlnR_v}$', fontsize=20)  
		plt.xscale('log')
		plt.yscale('log')  
		plt.title('Void Size Function')
		plt.grid(True)
		plt.legend()
		if params['save_fig']:
			plt.savefig('void_size_function.jpg')
		plt.show()

	return [mids,normalized_counts]