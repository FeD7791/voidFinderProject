from abc import ABC, abstractmethod

from .box import DataBox
from .voids import Voids


class ModelABC(ABC):

    def __init__(self):
        pass

    def find(self, databox: DataBox):
        preprocess_parameters = self.preprocess(databox)
        model_find_parameters = self.model_find(preprocess_parameters)
        voids_tuple, extra = self.build_voids(model_find_parameters)

        voids = Voids(
            method=type(self).__name__,
            tracers=databox.box,
            tracers_in_voids=voids_tuple,
            extra=extra,
        )

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
