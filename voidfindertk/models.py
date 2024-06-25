
from abc import ABC, abstractmethod

from attrs import define, field

from .box import DataBox



class ModelVoid(ABC):
    def __init__(self):
        pass

    @property
    @abstractmethod
    def tracers(self):
        pass 

    @property
    @abstractmethod
    def voids(self):
        pass

    @abstractmethod
    def voids_numbers(self):
        pass

    @abstractmethod
    def void_of(self,tracer):
        pass

class ModelABC(ABC):

    def __init__(self):
        pass

    def find(self, databox:DataBox):
        preprocess_parameters = self.preprocess(databox)
        model_find_parameters = self.model_find(preprocess_parameters)
        voids = self.build_voids(model_find_parameters)
        return voids

    @abstractmethod
    def preprocess(self, databox):
        pass

    @abstractmethod
    def model_find(self, preprocess_parameters):
        pass

    @abstractmethod
    def build_voids(self, model_find_parameters):
        pass


    # @abstractmethod
    # def get_void_mass(self, void_box, llbox):
    #     pass
