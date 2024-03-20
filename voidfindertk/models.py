import abc
from abc import abstractmethod, ABC


class ModelABC(ABC):

    def find(self,x,y,z):
        llbox = self.preprocess(x,y,z)
        voids = self.model_find(llbox)
        vb = self.mk_vbox(x,y,z,voids,llbox)
        return vb

    @abstractmethod
    def preprocess(x,y,z):
        pass
    @abstractmethod
    def model_find(llbox):
        pass
    @abstractmethod
    def mk_vbox(x,y,z,voids,llbox):
        pass

