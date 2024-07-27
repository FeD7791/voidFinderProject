from . import _processing

from ..zobov import Files, ZobovVF



class DiveVF(ZobovVF):

    def __init__(self, *, dive_coso, **kwargs):
        super().__init__(**kwargs)
        self._dive_coso = dive_coso

    @property
    def dive_coso(self):
        return self._dive_coso

    def build_voids(self, model_find_parameters):
        particle_by_voids, extra = super().build_voids(model_find_parameters)

        # Get void properties
        void_properties = extra["zobov_voids"]
        # Get box
        box = extra["box"]
        # Get tracer volumes
        tracer_volumes = _processing.read_volume_file(
            filename=Files.VOLUME_RAW
            )
        # get center and radii
        centers, radii = _processing.get_center_and_radii(
            void_properties=void_properties,
            tracer_volumes=tracer_volumes,
            particle_by_voids=particle_by_voids,
            box=box,
        )
        return radii, centers
