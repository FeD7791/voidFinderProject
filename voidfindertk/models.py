import abc
from attrs import define,field
from abc import abstractmethod, ABC
from . import data_box

@define 
class VoidMetrics():
    void_box_coordinates = field(init=False)
    void_mass = field(init=False)
    tracers = field(init=False)
    voids = field(init=False)


class ModelABC(ABC):

    def __init__(self):
        pass
    def find(self,box):
        llbox = self.preprocess(box)
        voids = self.model_find(llbox)
        vb = self.mk_vbox(voids,llbox)
        mass = self.get_void_mass(voids,llbox)
        metrics = VoidMetrics(**{
            'void_box_coordinates':vb,
            'void_mass':mass,
            'tracers':llbox.box,
            'voids':voids
        })
        return metrics

    @abstractmethod
    def preprocess(self, databox):
        pass
    @abstractmethod
    def model_find(self, llbox):
        pass
    @abstractmethod
    def mk_vbox(self,voids,llbox):
        pass
    @abstractmethod
    def get_void_mass(self,void_box,llbox):
        pass

