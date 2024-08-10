
from . import _postprocessing
from ..zobov import Names, ZobovVF


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
        ratio=1.5,
        initial_radius=True,
        delta_r=[17.0, 150.0],
        threshold=0.3,
        overlap_criterion = True,
        **kwargs,
    ):

        super().__init__(**kwargs)
        self._ratio = ratio
        self._initial_radius = initial_radius
        self._delta_r = delta_r
        self._threshold = threshold
        self._overlap_criterion = overlap_criterion

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
    
    @property
    def threshold(self):
        return self._overlap_criterion

    def build_voids(self, model_find_parameters):
        tracers_in_voids, extra = super().build_voids(model_find_parameters)
        # Get working directory
        run_work_dir = extra["files_directory_path"]
        # Get void properties
        void_properties = extra["void_properties"]
        # Get box
        box = model_find_parameters["box"]
        # Get tracer volumes
        tracer_volumes = _postprocessing.read_volume_file(
            filename=run_work_dir / _Files.VOL_FILE_RAW
        )
        # get center and radii
        radii, centers = _postprocessing.get_center_and_radii(
            void_properties=void_properties,
            tracer_volumes=tracer_volumes,
            tracers_in_voids=tracers_in_voids,
            box=box,
        )
        # Save to file
        # 1) center and radii
        _postprocessing.save_r_eff_center(
            centers=centers,
            r_eff=radii,
            path=str(run_work_dir / _Files.XYZ_R_EFF_VOIDS_FILE),
        )
        # 2) xyz box
        _postprocessing.save_xyz_tracers(
            box=box, path=str(run_work_dir / _Files.XYZ_TRACERS_FILE)
        )
        # Perform cleaning
        _postprocessing.cbl_cleaner(
            file_voids=str(run_work_dir / _Files.XYZ_R_EFF_VOIDS_FILE),
            file_tracers=str(run_work_dir / _Files.XYZ_TRACERS_FILE),
            ratio=self._ratio,
            initial_radius=self._initial_radius,
            delta_r=self._delta_r,
            threshold=self._threshold,
            output_path=str(run_work_dir / _Files.CLEANED_CATALOGUE),
            ol_crit=self._overlap_criterion
        )
        # Get tracers in voids
        tinv_cleaned_catalogue = _postprocessing.get_tracers_in_voids(
            box=box, cbl_cleaned_path=str(
                run_work_dir / _Files.CLEANED_CATALOGUE
                )
        )
        # Updating extra with void_properties
        extra["void_properties"] = _postprocessing.get_dive_void_properties(
            cleaned_catalogue_path = run_work_dir / _Files.CLEANED_CATALOGUE
            )
        return tuple(tinv_cleaned_catalogue), extra
