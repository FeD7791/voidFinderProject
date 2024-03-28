import pandas as pd
import numpy as np
import math

def analysis_voids(void_out,n_items,param): #param: a param like 'rad', 'volume' property of a void
											#n_items: the most relevant n_items

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


def scl(rho_med,delta,n_tracers):
	rad = (3*n_tracers/(4*np.pi*(delta+1)*rho_med))**(1/3)
	return np.log10(rad)



def void_size_function(delta:float,void_out):
	vol = (math.ceil(np.max(void_out.x)))**3

	if void_out._get_type() == 'SphericalVoids':
		rad = void_out.rad.values
	
	rhomed = void_out._len()/vol 
	n1 = range(1,6,1)
	n2 = range(6,12,2)
	n3 = range(12,62,10)

	mxlg = math.log10(np.max(rad))

	br1 = [scl(rho_med=rhomed,delta=delta,n_tracers=n) for n in n1]
	br2 = [scl(rho_med=rhomed,delta=delta,n_tracers=n) for n in n2]
	br3 = [scl(rho_med=rhomed,delta=delta,n_tracers=n) for n in n3]
	br4 = list(np.linspace(scl(rho_med=rhomed,delta=delta,n_tracers=62), mxlg,50))
	br = br1 + br2 + br3 + br4

	#Build array of bins widths
	dlogr1 = np.array([0] + br )
	dlogr2 = np.array(br + [0])
	dlogr = list(dlogr2 - dlogr1)
	dlogr.pop(0)
	dlogr.pop(-1)

	#Middle elements of bins
	mids = list(np.array(br) + np.array(dlogr + [0])/2)
	mids.pop(-1)

	#Get bin counts for each bin
	counts = np.histogram(a=rad, bins=br)

	#Calculate Density
	density = counts/np.array(dlogr)/vol

	#Transform back mids
	x = 10**np.array(mids)
	y = density
	return [x,y]