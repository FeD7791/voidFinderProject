from voidfindertk import zobov,io
import pathlib
import os
import joblib

dataset_path = pathlib.Path("./datasets")
workdir_path = pathlib.Path(os.path.abspath("./runz"))
box = io.read_table(
    dataset_path/"halos_ascii_1000_1024_npmin_10_z0.51.dat",
    names=["m","x", "y", "z", "vx", "vy", "vz"])


model = zobov.ZobovVF(
    density_threshold=0.1,
    box_size=1000,
    workdir=workdir_path
)
run_work_dir = pathlib.Path('/home/jorgefederico/updates/vftk_1109/voidFinderProject/runz/tmpqt51ba9_2024-09-23T19:33:00.395859+00:00')
model_find_parameters = {"run_work_dir": run_work_dir, "box": box}
void = model.build_voids(model_find_parameters=model_find_parameters)
# void = model.find(box)
# joblib.dump(voids, workdir / f"popcorn_{densth}.pkl")