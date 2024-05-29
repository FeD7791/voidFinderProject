import numpy as np
import ctypes
import pandas as pd
from . import box
from .models import DataBox
import math
import matplotlib.pyplot as plt
import scipy



def ctypes_input_params_builder(finder_type,**kwargs):
    if finder_type == 'spherical':
        # Input params
        params = {
        #Double Values
        "RadIncrement" : 0. , #cuanto incrementa la esfera a la hora de buscar con la random walk
        "BoxSize": 1000 ,
        "MaxRadiusSearch": 40.0, #Radio maximo que pienso encontrar en el catalogo de voids, muchas veces te frena la densidad. Tiene asociado un error que implica aumentar el maximo
        "ProxyGridSize": 5.0, # tamanio de la grilla para identificar vecinos el diametro del void debe ser masomenos el tamanio de la caja 
        "FracRadius" : 0.5 , # cuanto es el tamanio del salto de la caminata , no tiene un impacto muy fuerte sobre el algoritmo
        "DeltaThreshold": -0.9,  
        "DeltaSeed": -0.7,
        "OverlapTol":0,  # Cuanto se permite el overlap de tools , si le aumento la tol y se intersectan 2 los 2 se conservan
        "Redshift" : 0.99 ,
        "OmegaMatter" : 0.25,
        "OmegaLambda" : 0.75,
        "Hubble" : 0.73,
        "FidOmegaMatter": 0.2  ,
        "FidOmegaLambda" : 0.8 ,
        "FidHubble" : 0.7, 
        "MinProfileDist": 0.5,
        "MaxProfileDist" : 3.0,
        "ScalePos": 1,
        "ScaleVel" : 1,
        "InnerShellVel": 0.8,
        "OuterShellVel": 1.2,
        #Int Values
        "FormatTracers": 0,
        "NumFiles": 32,
        "NumRanWalk" : 75 , 
        "OMPcores": 8,
        "RSDist" : 0,
        "GDist": 0,
        "WriteProfiles" : 0 ,
        "NumProfileBins" : 100,   
            
        }
        #Replace params as needed
        for key,value in kwargs.items():
            params[key] = value
        
        # Load Input Params
        params1 = dict(list(params.items())[0:21]) # Params that are double type
        params2 = dict(list(params.items())[21:29]) # Params that are int type

        # Build ctypes structure for params
        p1 = [(key,ctypes.c_double) for key,value in params1.items()] 
        p2 = [(key,ctypes.c_int) for key,value in params2.items()]

        class InputParams(ctypes.Structure):
            _fields_ = p1 + p2
    
    return {'InputParams_class':InputParams,'params_dict':params}

def ctypes_output_builder(finder_type):
    if finder_type == 'spherical':
        # Create stucts --> class for output
        class voids(ctypes.Structure):
            _fields_ = [
                ("n_voids", ctypes.c_int),
                ("Rad", ctypes.c_float),
                # ("Rini", ctypes.c_float),
                # ("Ini", ctypes.c_float * 1),
                ("Pos", ctypes.c_float * 3),
                ("Vel", ctypes.c_float * 3),
                ("Dtype", ctypes.c_float),
                ("Delta", ctypes.c_float),
                ("Poisson", ctypes.c_float),
                # ("Dist4", ctypes.c_float),
                # ("ToF", ctypes.c_bool),
                ("Nran", ctypes.c_int),
            ]

        class voidArray(ctypes.Structure):
            _fields_ = [("voids", voids * 1)]
    return {'voidArray_class':voidArray}

def process_output_from_finder(finder_type, array_of_voids):
    va = array_of_voids
    if finder_type == 'spherical':
        va_arr = np.ctypeslib.as_array(va, shape=(10,))

        n_voids = va_arr[0][0][0][0] - 1

        full_arr = np.ctypeslib.as_array(va, shape=(n_voids,))

        radius = []
        x_coord = []
        y_coord = []
        z_coord = []
        vx_coord = []
        vy_coord = []
        vz_coord = []
        dtype = []
        delta = []
        poisson = []
        nran = []
        for i in range(0, n_voids, 1):
            radius.append(full_arr[i][0][0][1])
            x_coord.append(full_arr[i][0][0][2][0])
            y_coord.append(full_arr[i][0][0][2][1])
            z_coord.append(full_arr[i][0][0][2][2])
            vx_coord.append(full_arr[i][0][0][3][0])
            vy_coord.append(full_arr[i][0][0][3][1])
            vz_coord.append(full_arr[i][0][0][3][2])
            dtype.append(full_arr[i][0][0][4])
            delta.append(full_arr[i][0][0][5])
            poisson.append(full_arr[i][0][0][6])
            nran.append(full_arr[i][0][0][7])
    
    return {
        'rad':radius, 
        'x_void':x_coord, 'y_void':y_coord, 'z_void':z_coord,
        'vel_x_void':vx_coord, 'vel_y_void':vy_coord, 'vel_z_void':vz_coord,
        'delta':delta, 'dtype':dtype, 'poisson':poisson, 'nran':nran
        }

def preprocess_data_box(databox,**kwargs):
    kwargs.setdefault('m_min', 0)
    kwargs.setdefault('m_max', len(databox.box))
    kwargs.setdefault('round_to', 0)
    b = databox.box
    ##
    df = pd.DataFrame(b.__dict__)
    df = df[df['m'] >= kwargs['m_min']]
    df = df[df['m'] <= kwargs['m_max']]
    df.reset_index(drop=True,inplace=True)
    df.drop(columns=['_len'],axis=1,inplace=True)
    df = df.round(kwargs['round_to']) #Rounding numbers of dataframe for more or less precision
    df.drop_duplicates(inplace=True, ignore_index=True) #Drop any remaining duplicates
    box2 = box.Box(**df.to_dict(orient='list'))
    return DataBox(box2)

#ANALISIS TOOLS
def volume_radii_conversion(volume):
	#volume must be a np array
	if type(volume).__name__ != 'ndarray':
		raise TypeError('volume must be of type ndarray')
	r = ((3/(4*np.pi))*volume)**(1/3)
	return r


def distance_void_tracer_classification(
		x_tracers,y_tracers,z_tracers,
		x_voids,y_voids,z_voids,rad_voids
		):
	pos_void = np.array([
				x_voids, 
				y_voids, 
				z_voids
				])
	pos_void = np.transpose(pos_void)
	arr = []
	for i in range(len(x_tracers)):
		pos_box = np.array(
			[[x_tracers[i], y_tracers[i], z_tracers[i]]]
		)
		d = scipy.spatial.distance.cdist(
			pos_box, pos_void, metric="euclidean"
		)  # Calculates the distance from a tracer to each void center
			# d is voids._void_len	
		rad_center_distance_comparing = np.greater(rad_voids, d[0]).astype(
			float
		)  # Compares radius of a void and distance from the void center to the particle 1 if rad>d # noqa: E501
		arr.append(rad_center_distance_comparing)
	# Transform to sparse matrix
	return scipy.sparse.csr_matrix(arr)



def join_box_void(box,voids, **kwargs):
		
		params = {
			'tol':0.0,
			}
		for key,value in kwargs.items():
			params[key] = value
		
		if type(voids).__name__ in ['SphericalVoids','PopCornSVF']:
			
			tolerance = params['tol'] * np.ones(voids.rad.shape)
			rad = voids.rad.value + tolerance

			output_sparse_matrix = distance_void_tracer_classification(
				x_tracers=box.x.value,
				y_tracers=box.y.value,
				z_tracers=box.z.value,
				x_voids=voids.x_void.value,
				y_voids=voids.y_void.value, 
				z_voids=voids.z_void.value,
				rad_voids=rad

			)

		if type(voids).__name__ in ['PopCornVoids','ZobovVoids']:
			if type(voids).__name__ == 'PopCornVoids':
				tracers_in_voids = voids.ID_subhalo

			if type(voids).__name__ == 'ZobovVoids':
				tracers_in_voids = voids._tracers_in_void 

			n_index =[]
			tracers = []

			for i in range(0,voids._void_len):
				n_index = n_index + list(i*np.ones(np.array(tracers_in_voids[i]).shape))
				tracers = tracers + tracers_in_voids[i]


			row = np.array(tracers)
			col = np.array(n_index)
			data = np.ones(row.shape)
			output_sparse_matrix = scipy.sparse.csr_matrix((data, (row, col)), shape=(box._len, voids._void_len))
		
		# if type(voids).__name__ == 'ZobovVoids':

		# 	# rad = volume_radii_conversion(voids.VoidVol.value)
		# 	# centers = find_zobov_void_centers(box=box, zobov_voids=voids)
		# 	# output_sparse_matrix = distance_void_tracer_classification(
		# 	# 	x_tracers=box.x.value,
		# 	# 	y_tracers=box.y.value,
		# 	# 	z_tracers=box.z.value,
		# 	# 	x_voids=centers['x'],
		# 	# 	y_voids=centers['y'], 
		# 	# 	z_voids=centers['z'],
		# 	# 	rad_voids=rad

		# 	# )

		return output_sparse_matrix



def analysis_voids(box,voids,sparse,**kwargs): #n_items: the most relevant n_items

	params = {
		'n_items' : len(voids)
	}
	for key,value in kwargs.items():
		params[key] = value
										
	df_voids = pd.DataFrame(voids.__dict__)
	df_box = pd.DataFrame(box.__dict__)
	
	if type(voids).__name__ == 'SphericalVoids':
		df_filtered_param = df_voids.sort_values(by='rad', ascending=False)
	if type(voids).__name__  == 'PopCornVoids':
		df_filtered_param = df_voids.sort_values(by='void_volume', ascending=False)
	if type(voids).__name__ == 'ZobovVoids':
		df_filtered_param = df_voids.sort_values(by='VoidVol', ascending=False)


	df_filtered_index_cut = df_filtered_param.index[:params['n_items']]

	parameter_voids = {'void_number':[], 'velocity_norm':[], 'n_tracers':[], 'type':type(voids).__name__ }

	for i in df_filtered_index_cut:
		tracer_indexes = list(sparse.getcol(i).tocoo().row) 
		parameter_voids['void_number'].append(i)
		parameter_voids['velocity_norm'].append(
			np.linalg.norm(
				df_box.loc[
					tracer_indexes,['vx','vy','vz']].sum(axis=0))/len(tracer_indexes))
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
		rad = volume_radii_conversion(voids.void_volume.value)
	if type(voids).__name__ == 'ZobovVoids':
		rad = volume_radii_conversion(voids.VoidVol.value)
	
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


def find_zobov_void_centers(box,zobov_voids):
    box_df = pd.DataFrame(box.__dict__)
    zobov_df = pd.DataFrame(zobov_voids.__dict__)
    zobov_df = zobov_df.sort_values(by='CoreParticle')
    core_particle_list = list(zobov_df['CoreParticle'])
    val = box_df.index.isin(core_particle_list)
    x = list(box_df[val]['x'])
    y = list(box_df[val]['y'])
    z = list(box_df[val]['z'])
    return {'CoreParticle':core_particle_list,'x':x,'y':y,'z':z}



def void_size_function_2(box,voids, **kwargs):
	params = {
		'n_bins':50, 
		}
	for key,value in kwargs.items():
		params[key] = value

	#Get radii
	if type(voids).__name__ in ['SphericalVoids','SVF']:
		rad = voids.rad.value
	if type(voids).__name__ == 'PopCornVoids':
		rad = volume_radii_conversion(voids.void_volume.value)
	if type(voids).__name__ == 'ZobovVoids':
		rad = volume_radii_conversion(voids.VoidVol.value)

	#Vol simulation
	vol = (math.ceil(np.max(box.x.value)))**3 

	#rho med of tracers
	rho_med = len(box)/vol

	#Seq
	n1 = list(np.arange(1,6,1))
	n2 = list(np.arange(6,11,2))
	n3 = list(np.arange(12,62,10))
	n = np.array(n1 + n2 + n3)

	#Delta
	delta = -0.9

	#scl
	#scl = np.log10((3/(4*np.pi)*n/rho_med/(1+delta))**(1/3))
	scl = np.log10((3*n/((1+delta)*(4*np.pi)*rho_med))**(1/3))

	#max log rad
	mxlg = np.log10(max(rad))

	#breaker br
	br = np.array(list(scl[1:len(scl)-1]) + list(np.linspace(max(scl),mxlg, params['n_bins'])))

	#x
	x = rad
	#histogram
	counts, bins0 = np.histogram(np.log10(x), bins=br)

	#dlogr
	dlogr = br[1:len(br)] - br[0:len(br)-1]

	#mids
	mids = (bins0[1:] + bins0[:-1]) / 2

	#Normalized density
	#density = counts/dlogr/vol
	density = counts/(dlogr*vol) 
	return [10**mids,density]


#COMPARISON TOOLS
def find_radius(box,voids,sparse):
    radius = []
    tracers = [[box.x.value[i],box.y.value[i],box.z.value[i]] for i in range(len(box))]
    #DEPENDE DEL VOIDFINDER
    ##CONFIGURAR PARA VARIOS VOID FINDERS
    voids_centers = [
        [
            voids.x_void.value[i],
            voids.y_void.value[i],
            voids.z_void.value[i]] for i in len(voids)]
    ##########################################################
    sparse = sparse.toarray() #Transform sparse matrix to array [M,N]
    for i in range(sparse.shape[1]):
        c = sparse[:,i][:, None] * tracers # filter the tracers in voids
        c = c[np.any(c, axis=1)] #Removes [0,0,0] rows
        c = np.vstack([c,voids_centers[i]]) # Adding void center to points
        distances = scipy.spatial.distance.pdist(c) #distances between points
        radius.append(distances.sort(-1)[-1]) #appending biggest distance to radius as candidate
    return radius