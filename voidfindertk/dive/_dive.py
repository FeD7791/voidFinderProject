from . import _processing

from ..zobov import ZobovVF, Names

class _Files:
    VOL_FILE_RAW = f"vol{Names.OUTPUT_VOZINIT}.dat"
    ADJ_FILE_RAW = f"adj{Names.OUTPUT_VOZINIT}.dat"
    R_EFF_CENTER_FILE = "r_eff_center.dat"

class DiveVF(ZobovVF):

    def __init__(self, *, dive_coso, **kwargs):
        super().__init__(**kwargs)
        self._dive_coso = dive_coso

    @property
    def dive_coso(self):
        return self._dive_coso

    def build_voids(self, model_find_parameters):
        tracers_in_voids, extra = super().build_voids(model_find_parameters)
        # Get working directory
        workdir = extra["files_directory_path"]
        # Get void properties
        void_properties = extra["zobov_voids"]
        # Get box
        dbox = super().preprocess()
        box = dbox.box
        # Get tracer volumes
        tracer_volumes = _processing.read_volume_file(
            filename= workdir / _Files.VOL_FILE_RAW
        )
        # get center and radii
        centers, radii = _processing.get_center_and_radii(
            void_properties=void_properties,
            tracer_volumes=tracer_volumes,
            tracers_in_voids=tracers_in_voids,
            box=box,
        )
        return centers

