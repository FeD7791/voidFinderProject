import abc
from abc import abstractmethod, ABC
from . import data_box


class ModelABC(ABC):

    def __init__(self):
        pass
    def find(self,databox):
        llbox = self.preprocess(databox)
        voids = self.model_find(llbox)
        vb = self.mk_vbox(databox,voids,llbox)
        return vb

    @abstractmethod
    def preprocess(databox):
        pass
    @abstractmethod
    def model_find(llbox):
        pass
    @abstractmethod
    def mk_vbox(databox,voids,llbox):
        pass

