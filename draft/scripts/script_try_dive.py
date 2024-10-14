from voidfindertk import io
from voidfindertk import dive, zobov
from voidfindertk import core
import joblib
import time
import pathlib

# Run full finder

dataset_path = pathlib.Path(
    "/home/jorgefederico/updates/vftk_1109/voidFinderProject/datasets/"
)

# Get Box
box = io.read_table(
    dataset_path / "halos_ascii_1000_1024_npmin_10_z0.51.dat",
    names=["m", "x", "y", "z", "vx", "vy", "vz"],
)


workdir = pathlib.Path(
    "/home/jorgefederico/updates/vftk_1109/voidFinderProject/runz"
)


def build_zobov_model(workdir, *, delta=0.2):

    model = zobov.ZobovVF(
        density_threshold=0.2,  # zobov parameter
        box_size=1000,
        workdir=workdir,
    )
    model_find_parameters = {
        "run_work_dir": model.workdir
        / "tmpuo2jislx2024-09-18T23:23:15.845781+00:00",
        "box": box,
    }
    tinv, centers, extra = model.build_voids(model_find_parameters)
    void = core.Voids(
        method="ZOBOV",
        box=box,
        tracers_in_voids_=tinv,
        centers_=centers,
        extra_=extra,
    )
    joblib.dump(void, f"zobov_model_{delta}.pkl")


def build_dive_model(box, workdir, *, threshold=0.2, op_crit=False, ratio=1.5):
    # model
    dive_model = dive.DiveVF(
        density_threshold=0.2,  # zobov parameter
        box_size=1000,
        workdir=workdir,
        delta_r=[0.0, 100.0],
        threshold=threshold,  # This would be a density contrast of -0.8
        overlap_criterion=op_crit,  # False: The void with the higher central density is rejected
        ratio=ratio,
    )

    model_find_parameters = {
        "run_work_dir": dive_model.workdir
        / "tmpuo2jislx2024-09-18T23:23:15.845781+00:00",  # "tmpv4suzrwr2024-09-16T00:16:06.078226+00:00",
        "box": box,
    }
    tinv, centers, extra = dive_model.build_voids(model_find_parameters)
    void = core.Voids(
        method="Dive",
        box=box,
        tracers_in_voids_=tinv,
        centers_=centers,
        extra_=extra,
    )
    joblib.dump(void, f"dive_model_t_{threshold}_rat_{ratio}.pkl")
