import os
import pathlib
import tempfile

from ..utils import create_run_working_directory, save_file_from_box

class _Paths:
    CURRENT = pathlib.Path(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    )
    SCVF = CURRENT / 'SCVF'

# def _get_scvf_params(scvf_file_path):
#     """
#     Get the parameters and default values from a VF.param 
#     reference file. This is used to get all the params for use
#     of the SCVF based on the VF.param default input file. The 
#     function will read the file an map the name and values to 
#     a dictionary

#     Parameters:
#         scvf_file_path(str): It is the path to the VF.param 
#         file that is used as reference to get the parameters

#     Returns:
#         A dictionary with the keys of the input params in the 
#         VF.param file and values as the values of those params
#     """
#     with open(scvf_file_path,'r') as f:
#         a = f.readlines()
#     parameters = {}
#     for i in a:
#         k = i.split()
#         if len(k) > 1 and k[0] != '%':
#             parameters[f'{k[0]}'] = k[1]
#     return parameters

def _get_params(*,class_parameters:list):
    """
    This function is used to build the parameters for the 
    SCVF void finder.

    Parameters
    ----------
        class_parameters(list): A list of input parameters of
        the SCVF class

    Returns
    -------
        params(dictionary): A dictionary with the keys as the
        parameter properties that identify each input used in the 
        VF.param file
    """
    # List of parameters used in VF.param file
    get_vf_params = [
        'OMPcores', 'BoxSize','MaxRadiusSearch','ProxyGridSize',
        'DeltaThreshold','DeltaSeed','OverlapTol','FormatTracers',
        'NumFiles','FileTracers','FileVoids','ScalePos','ScaleVel',
        'FracRadius','NumRanWalk','RadIncrement','RandomSeed',
        'RSDist','Redshift','OmegaMatter','OmegaLambda','Hubble',
        'GDist','FidOmegaMatter','FidOmegaLambda','FidHubble',
        'WriteProfiles','MinProfileDist','MaxProfileDist',
        'NumProfileBins','PathProfiles','InnerShell','OuterShell'
        ]
    # get_vf_params = list(_get_scvf_params(SCVF / 'params' / VF.param).keys())
    # Get the param values from current class SCVF 
    scvf_params = [
        p for p in class_parameters 
        if not p.starts_with('__') and p not in
        [
            'model_find',
            'preprocess',
            'workdir',
            'vf_param_path',
            'main_executable_path'
            ]
        ]
    # Dictionary that holds the parameters used as input in the finder
    params = {}
    # Build dictionary from the previous param names and param values
    for param_name,param_value in zip(get_vf_params,scvf_params):
        params[f"{param_name}"] = param_value
    return params


def _build_vf_param_file(*,vf_param_path,**kwargs):
    """
    This function will open a file with name VF.param in some 
    location (vf_param_path) and will save any parameters in
    kwargs

    Parameters: 
        vf_param_path(pathlib.Path) : Is the path to the location
        directory of the file. 

    """
    with open(vf_param_path / 'VF.param', 'w') as conf:
        for key,value in kwargs.items():
            conf.writelines([key," ",value,"\n"])

class SCVF:
    def __init__(
        self,
        *,
        box_size=500, #box should get this parameter
        max_radius_search=40.0,
        proxy_grid_size=5.0,
        delta_threshold=-0.9,
        delta_seed=-0.7,
        overlap_tol=0.0,
        format_tracers=0,
        num_files=32,
        file_tracers=None,
        file_voids='test.dat',
        scale_pos=1.0,
        scale_vel=1.0,
        frac_radius=0.5,
        num_ran_walk=75,
        rad_increment=0.0,
        random_seed=1234,
        rsd_dist=0,
        red_shift=0.99,
        omega_matter=0.25,
        omega_lambda=0.75,
        hubble=0.73,
        g_dist=0,
        fid_omega_matter=0.2,
        fid_omega_lambda=0.8 ,
        fid_hubble=0.7,
        write_profile=0,
        min_profile_dist=0.5,
        max_profile_dist=3.0,
        num_profile_bins=100,
        path_profiles='data/profiles/',
        inner_shell=0.8,
        outher_shell=1.2,
        # Paths
        workdir = None,
        vf_param_path = None,
        main_executable_path = None
        ):
        self.box_size=box_size, #box should get this parameter
        self.max_radius_search=max_radius_search,
        self.proxy_grid_size=proxy_grid_size,
        self.delta_seed=delta_seed,
        self.delta_threshold=delta_threshold,
        self.overlap_tol=overlap_tol,
        self.format_tracers=format_tracers,
        self.num_files=num_files,
        self.file_tracers=file_tracers,
        self.file_voids=file_voids,
        self.scale_pos=scale_pos,
        self.scale_vel=scale_vel,
        self.frac_radius=frac_radius,
        self.num_ran_walk=num_ran_walk,
        self.rad_increment=rad_increment,
        self.random_seed=rad_increment,
        self.rsd_dist=rsd_dist,
        self.red_shift=red_shift,
        self.omega_matter=omega_matter,
        self.omega_lambda=omega_lambda,
        self.hubble=hubble,
        self.g_dist=g_dist,
        self.fid_omega_matter=fid_omega_matter,
        self.fid_omega_lambda=fid_omega_lambda ,
        self.fid_hubble=fid_hubble,
        self.write_profile=write_profile,
        self.min_profile_dist=min_profile_dist,
        self.max_profile_dist=max_profile_dist,
        self.num_profile_bins=num_profile_bins,
        self.path_profiles=path_profiles  ,
        self.inner_shell=inner_shell,
        self.outher_shell=outher_shell,
        self.random_seed = random_seed
        # paths
        # workdir holds the working directory
        self.workdir = pathlib.Path(
            tempfile.mkdtemp(prefix=f"vftk_{type(self).__name__}_")
            if workdir is None
            else workdir
        )
        # Location of the executable main.x that runs SCVF
        self.main_executable_path = pathlib.Path(
            _Paths.SCVF if main_executable_path is None else main_executable_path
        )
        # Path to save the VF.param file used to run the SCVF void finder
        self.vf_param_path = pathlib.Path(
            _Paths.SCVF /'params' if vf_param_path is None else vf_param_path
        )
        # Path to save the traces file
        self.file_tracers = pathlib.Path(
            file_tracers if file_tracers is not None else workdir
        )


    def preprocess(self,box):
        """
        This method is the previous step of preparation to run the SCVF
        It will first build the VF.params file needed to run the SCVF.
        The VF.params file is used to get the list of parameters  
        """
        
        params = _get_params(dir(self))

        # Build parameters file VF.param in the location vf_param_path using params
        _build_vf_param_file(vf_param_path=self.vf_param_path,**params)

        # Crete the input file if box was provided
        if self.file_tracers is not None:
            save_file_from_box.box_to_csv(box,self.file_tracers)
        else:
            save_file_from_box.box_to_csv(box,self.workdir)


    def model_find(self,box):
        # create sandbox for this run
        run_work_dir = create_run_working_directory()
        

        

    # def preprocess(self, databox, **kwargs):
    #     databox = preprocess_data_box(databox=databox, **kwargs)
    #     return databox

    # def model_find(self, llbox,**kwargs):
    #     sp_void = spherical_void_finder(llbox.box,**kwargs)
    #     return {'voids':sp_void}
    
    # def mk_vbox(self, voids,llbox):
    #     voids = voids['voids']
    #     # databox.box.__dict__.pop('_len')
    #     # voids.__dict__.pop('_void_len')
    #     #sparse_m = Classifier(**databox.box.__dict__, **voids.__dict__)._sparse_matrix(0.0)
    #     box_void_sparse = join_box_void(llbox.box, voids, tol=0.0)
    #     return box_void_sparse