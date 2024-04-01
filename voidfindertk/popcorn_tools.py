import uttr
import numpy as np
from .models import ModelABC
from .analysis_tools import join_box_void
from astropy import units as u

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


