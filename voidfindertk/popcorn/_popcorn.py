from ..svf_popcorn import PopCornVF, Paths, FileNames
from . import _wrapper

import attr


@attr.define
class PopCorn(PopCornVF):
    _shot_noise_threshold = attr.field(default=20)

    @property
    def shot_noise_threshold(self):
        return self._shot_noise_threshold

    def model_find(self, databox):
        parameters = super().model_find(databox)
        run_work_dir = parameters["run_work_dir"]
        box = parameters["box"]

        # Before continuing minradius must be re-configured
        _wrapper.read_and_modify_config(
            config_file_path=run_work_dir / FileNames.CONFIG,
            section="INPUT_PARAMS",
            parameter="MINRADIUS",
            new_value=str(self._shot_noise_threshold)
        )
        _wrapper.popcorn_void_finder(
            mpi_flags=self._mpi_flags,
            bin_path=Paths.SVF,
            conf_file_path=run_work_dir / FileNames.CONFIG,
            work_dir_path=run_work_dir
            )
        _wrapper.compute_intersects(
            bin_path=Paths.SVF,
            conf_file_path=run_work_dir / FileNames.CONFIG,
            work_dir_path=run_work_dir
        )
        _wrapper.clean_duplicates(
            bin_path=Paths.SVF,
            conf_file_path=run_work_dir / FileNames.CONFIG,
            work_dir_path=run_work_dir
        )
        return {"run_work_dir": run_work_dir, "box": box}
    
    def build_voids(self, model_find_parameters):
        return 1,2,3

















# @attrs.define
# class A:
#     a: int = 1
#     def howdy(self):
#         print('howdy neighbor?!')


# @attrs.define
# class B(A):
#     b: str = "default"
#     c: int = 56
#     def __attrs_post_init__(self):
#         self.a = self.c
#     def howdy_new(self):
#         super().howdy()
#         print("Hello neighbor")
