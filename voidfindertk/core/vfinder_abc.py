from abc import ABC, abstractmethod

from .box import DataBox
from .voids import Voids


class VoidFinderABC(ABC):

    def __init__(self):
        pass

    def find(self, databox: DataBox):
        preprocess_parameters = self.preprocess(databox)
        model_find_parameters = self.model_find(preprocess_parameters)
        tracers_in_voids, centers, box, extra = self.build_voids(
            model_find_parameters
        )

        voids = Voids(
            method=type(self).__name__,
            box=databox.box,
            tracers_in_voids=tracers_in_voids,
            centers=centers,
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