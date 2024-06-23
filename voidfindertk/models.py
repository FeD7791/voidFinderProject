
from abc import ABC, abstractmethod

from attrs import define, field

from .box import DataBox


# class DataBox:
#     box = None

#     def __init__(self, aux_box):
#         self.box = aux_box

@define
class VoidMetrics:
    void_box_coordinates = field(init=False)
    void_mass = field(init=False)
    tracers = field(init=False)
    voids = field(init=False)


class ModelABC(ABC):

    def __init__(self):
        pass

    def find(self, databox:DataBox):
        preprocess_parameters = self.preprocess(databox)
        model_find_parameters = self.model_find(preprocess_parameters)
        build_void_parameters = self.build_void(model_find_parameters)
        # vb = self.mk_vbox(voids, llbox)
        # mass = self.get_void_mass(voids, llbox)
        # metrics = VoidMetrics(
        #     **{
        #         "void_box_coordinates": vb,
        #         "void_mass": mass,
        #         "tracers": llbox.box,
        #         "voids": voids,
        #     }
        # )
        # return metrics
        return voids

    @abstractmethod
    def preprocess(self, databox):
        pass

    @abstractmethod
    def model_find(self, prep_box):
        pass

    @abstractmethod
    def build_void(self, voids, llbox):
        pass

    # @abstractmethod
    # def get_void_mass(self, void_box, llbox):
    #     pass
