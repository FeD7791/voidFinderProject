import uttr
import numpy as np
from ..models import ModelABC
from ..analysis_tools import join_box_void
from astropy import units as u
import os
import subprocess
import pandas as pd
from configparser import ConfigParser

@uttr.s(repr=False)
class PopCornVoids:
    void_id = uttr.ib()
    number_members = uttr.ib()
    void_volume = uttr.ib(converter=np.array,unit=u.Mpc**3)
    number_subhalos = uttr.ib() 
    internal_flag = uttr.ib()
    spheres = uttr.ib()
    ID_subhalo = uttr.ib()
    _void_len = uttr.ib(init=False)

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        lengths = set(())
        for e in (
            self.void_id,
            self.number_members,
            self.void_volume ,
            self.number_subhalos,
            self.internal_flag,
            self.spheres,
            self.ID_subhalo
        ):
            lengths.add(len(e))

        if len(lengths) != 1:
            raise ValueError("Arrays should be of the same size")

        super().__setattr__("_void_len", lengths.pop())

    def __len__(self):
        """Length method.

        Returns
        -------
            int
                the number of elements in SphericalVoids
        """
        return self._void_len
    
    def _slice(self,min,max):
        popvoid = {
            'void_id':self.void_id[min:max],
            'number_members':self.number_members[min:max],
            'void_volume':self.void_volume[min:max],
            'number_subhalos':self.number_subhalos[min:max],
            'internal_flag':self.internal_flag[min:max],
            'spheres':self.spheres[min:max],
            'ID_subhalo':self.ID_subhalo[min:max]
        }
        return PopCornVoids(**popvoid)
    
    def __repr__(self):
        """Representation method.

        Returns
        -------
            str
                Name plus number of points in the box
        """
        cls_name = type(self).__name__
        length = len(self)
        return f"<{cls_name} size={length}>"
    

class PopCornVF(ModelABC):
    def __init__(self):
        pass

    def preprocess(self, databox):
        return databox

    def model_find(self, llbox):
        None #Not implemented yet
        # sp_void = popcorn_void_finder(llbox.box) #This part calls the c method
        # return {'voids':sp_void}
    
    def mk_vbox(self, databox,voids,sparse,llbox):
        voids = voids['voids']
        sparse_m = join_box_void(databox,voids,tol=0.0)
        return sparse_m
    


@uttr.s(repr=False)
class PopCornSVF:
    radius = uttr.ib(converter=np.array,unit=u.Mpc)
    x_void = uttr.ib(converter=np.array,unit=u.Mpc)
    y_void = uttr.ib(converter=np.array,unit=u.Mpc)
    z_void = uttr.ib(converter=np.array,unit=u.Mpc)
    delta = uttr.ib(converter=np.array)

def popcorn_svf_input_data_builder(box, **kwargs):
    kwargs.setdefault("sep",'\t')
    kwargs.setdefault("index",False)
    kwargs.setdefault("header",False)
    path_file = os.path.dirname(os.path.realpath(__file__))
    df = pd.DataFrame(box.__dict__)
    df.drop(labels=['_len'], axis=1,inplace=True)
    df.to_csv(
        os.path.join(path_file,'popcorn','popcorn_void_finder','Source','box.dat'), 
        **kwargs) #Save the box in a file 
    
def config_file_maker(params, path):
    '''
    path: place where the input file is (always has vars.conf as name)
    '''
    params = dict([(key,str(value)) for key,value in params.items()]) #Transform non str items in strings
    config_element = ConfigParser()
    config_element.optionxform = str
    config_element['INPUT_PARAMS'] = params
    with open(os.path.join(path,'vars.conf'), 'w') as conf:
        config_element.write(conf)

def spherical_popcorn_void_finder():
    #Paths
    path = os.path.dirname(os.path.realpath(__file__))
    path_source = os.path.join(path,'popcorn','popcorn_void_finder','Source')
    os.chdir(path_source)
    subprocess.run(["mpirun", "-np", "1","--bind-to","none", "./svf", "./configuration/vars.conf"])

class SphericalPopCornVF(ModelABC):
    def __init__(self):
        pass

    def preprocess(self, databox, **kwargs):
        box = databox.box
        popcorn_svf_input_data_builder(box, **kwargs)
        return databox

    def model_find(self, llbox, **kwargs):
        #Paths
        path = os.path.dirname(os.path.realpath(__file__))
        path_config_file = os.path.join(
            path,
            'popcorn',
            'popcorn_void_finder',
            'Source',
            'configuration')
        #Default Params
        kwargs.setdefault('TRSFILE',
                          os.path.join(
                              path,
                              'popcorn',
                              'popcorn_void_finder',
                              'Source',
                              'box.dat'))
        kwargs.setdefault('FILEFMT','ASCII')
        kwargs.setdefault('NUM_FILE', 1)
        kwargs.setdefault('POPFILE', '../MWE/clean_popvds.dat')
        kwargs.setdefault('SPHFILE', '../MWE/sphvds.dat')
        kwargs.setdefault('AUXFILES' , 'true')
        kwargs.setdefault('RAWPOPFILE','../MWE/raw_popvds.dat')
        kwargs.setdefault('RAWSPHFILE', '../MWE/raw_sphv.dat')
        kwargs.setdefault('PAIRSFILE', '../MWE/pairs_raw_popvds.dat')
        kwargs.setdefault('BOXSIZE' , 1000.0)
        kwargs.setdefault('DENSTH', -0.9)
        kwargs.setdefault('MINRADIUS', 5)
        kwargs.setdefault('MAXRADIUS', 100)
        kwargs.setdefault('EPS', 0.00001)
        kwargs.setdefault('MASSMIN', 50)

        # Create config file
        config_file_maker(params=kwargs, path=path_config_file)
        # Use finder
        spherical_popcorn_void_finder()
        # Get output data
        sphvds = pd.read_csv(kwargs['SPHFILE'], sep='\s+')
        sphvds.columns = [
            'index',
            'radius',
            'x_center',
            'y_center',
            'z_center',
            'delta'
            ] #preguntar a dante cual es la salida
        # Tansform output to PopCornSVF class
        sphvds.drop(labels=['index'], axis=1, inplace=True)
        popcorn_svf = PopCornSVF(**sphvds.to_dict(orient='list'))
        return popcorn_svf
    
    def mk_vbox(self, databox,voids,sparse,llbox):
        voids = voids['voids']
        sparse_m = join_box_void(databox,voids,tol=0.0)
        return sparse_m



