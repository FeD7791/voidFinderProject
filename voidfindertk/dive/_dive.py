from . import _processing

import attr

from ..zobov import ZobovVF, Names

class _Files:
    VOL_FILE_RAW = f"vol{Names.OUTPUT_VOZINIT}.dat"
    ADJ_FILE_RAW = f"adj{Names.OUTPUT_VOZINIT}.dat"
    XYZ_R_EFF_VOIDS_FILE = "xyz_r_eff_file.txt"
    XYZ_TRACERS_FILE = "xyz_tracers.txt"
    CLEANED_CATALOGUE = "cleaned_catalogue.txt"


class DiveVF(ZobovVF):

    def __init__(
            self,
            *,
            ratio,
            initial_radius,
            delta_r,
            threshold,
            **kwargs
            ):
        
        super().__init__(**kwargs)
        self._ratio = ratio
        self._initial_radius = initial_radius
        self._delta_r = delta_r
        self._threshold = threshold


    @property
    def ratio(self):
        return self._ratio
    @property
    def initial_radius(self):
        return self._initial_radius
    @property
    def delta_r(self):
        return self._delta_r
    @property
    def threshold(self):
        return self._threshold

    def build_voids(self, model_find_parameters):
        tracers_in_voids, extra = super().build_voids(model_find_parameters)
        # Get working directory
        run_work_dir = extra["files_directory_path"]
        # Get void properties
        void_properties = extra["zobov_voids"]
        # Get box
        box = model_find_parameters["box"]
        # Get tracer volumes
        tracer_volumes = _processing.read_volume_file(
            filename= run_work_dir / _Files.VOL_FILE_RAW
        )
        # get center and radii
        centers, radii = _processing.get_center_and_radii(
            void_properties=void_properties,
            tracer_volumes=tracer_volumes,
            tracers_in_voids=tracers_in_voids,
            box=box,
        )
        # Save to file
        # 1) center and radii
        _processing.save_r_eff_center(
            centers=centers,
            r_eff=radii,
            path=str(run_work_dir / _Files.XYZ_R_EFF_VOIDS_FILE)
            )
        # 2) xyz box
        _processing.save_xyz_tracers(
            box=box,
            path=str(run_work_dir / _Files.XYZ_TRACERS_FILE)
        )
        # Perform cleaning
        _processing.cbl_cleaner(
            file_voids=str(run_work_dir / _Files.XYZ_R_EFF_VOIDS_FILE),
            file_tracers=str(run_work_dir / _Files.XYZ_TRACERS_FILE),
            ratio=self._ratio,
            initial_radius=self._initial_radius,
            delta_r=self._delta_r,
            threshold=self._threshold,
            output_path=str(run_work_dir / _Files.CLEANED_CATALOGUE)
        )
        return centers

