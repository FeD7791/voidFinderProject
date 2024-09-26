from voidfindertk import io, popcorn, svf_popcorn, core
import os
import pathlib
import matplotlib.pyplot as plt
import joblib

dataset_path = pathlib.Path("./datasets")
workdir_path = pathlib.Path(os.path.abspath("./run_real_popcorn"))
dbox = io.read_table(
    dataset_path/"halos_ascii_1000_1024_npmin_10_z0.51.dat",
    names=["m","x", "y", "z", "vx", "vy", "vz"])
box = dbox.box

run_w_08 = "tmpipk1m12i2024-09-14T23:29:54.863443+00:00"

def build_svf_popcorn(box,workdir,run_w,*,densth=-0.8):
    model = svf_popcorn.SVFPopCorn(
    workdir=workdir,
    boxsize=1000,
    densth=densth,
    minradius = 10,
    maxradius = 100,
    )
    tracers,centers,extra = model.build_voids({
        "run_work_dir": model.workdir / run_w,
        "box": box})
    voids = core.Voids(
        method='SVF',
        box = box,
        tracers_in_voids_ = tracers,
        centers_ = centers,
        extra_ = extra
    )
    joblib.dump(voids, workdir / f"svf_{densth}.pkl")


def real_popcorn(box,workdir,run_w,*,densth=-0.8):

    model = popcorn.PopCorn(
        workdir=workdir,
        boxsize=1000,
        densth=densth,
        minradius = 10,
        maxradius = 100,
        shot_noise_threshold = 19
    )
    tracers,centers,extra = model.build_voids({
        "run_work_dir": model.workdir / run_w,
        "box": box})
    voids = core.Voids(
        method='POPCORN',
        box = dbox.box,
        tracers_in_voids_ = tracers,
        centers_ = centers,
        extra_ = extra
    )
    joblib.dump(voids, workdir / f"popcorn_{densth}.pkl")



