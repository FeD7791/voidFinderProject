import uttr
import numpy as np
from astropy import units as u


@uttr.s(repr=False)
class ZobovVoids:
    Void_number = uttr.ib(converter=np.array)
    File_void_number = uttr.ib(converter=np.array)
    CoreParticle = uttr.ib(converter=np.array)
    CoreDens = uttr.ib(converter=np.array)
    ZoneVol = uttr.ib(converter=np.array)
    Zone_number_part = uttr.ib(converter=np.array)
    Void_number_Zones = uttr.ib(converter=np.array)
    VoidVol = uttr.ib(converter=np.array,unit=u.Mpc**3)
    Void_number_Part = uttr.ib(converter=np.array)
    VoidDensContrast = uttr.ib(converter=np.array)
    VoidProb = uttr.ib(converter=np.array)
    _void_len = uttr.ib(init=False)

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        lengths = set(())
        for e in (
            self.Void_number,
            self.File_void_number,
            self.CoreParticle,
            self.CoreDens,
            self.ZoneVol,
            self.Zone_number_part,
            self.Void_number_Zones,
            self.VoidVol,
            self.Void_number_Part,
            self.VoidDensContrast,
            self.VoidProb
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
        zobovvoid = {
            'Void_number':self.Void_number[min:max],
            'File_void_number':self.File_void_number[min:max],
            'CoreParticle':self.CoreParticle[min:max],
            'CoreDens':self.CoreDens[min:max],
            'ZoneVol':self.ZoneVol[min:max],
            'Zone_number_part':self.Zone_number_part[min:max],
            'Void_number_Zones':self.Void_number_Zones[min:max],
            'VoidVol':self.VoidVol[min:max],
            'Void_number_Part':self.Void_number_Part[min:max],
            'VoidDensContrast':self.VoidDensContrast[min:max],
            'VoidProb':self.VoidProb[min:max]
        }
        return ZobovVoids(**zobovvoid)
    
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