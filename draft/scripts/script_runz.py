from voidfindertk import zobov, io
import pathlib
import os
import joblib

dataset_path = pathlib.Path("./datasets")
workdir_path = pathlib.Path(os.path.abspath("./runz"))
dbox = io.read_table(
    dataset_path / "halos_ascii_1000_1024_npmin_10_z0.51.dat",
    names=["m", "x", "y", "z", "vx", "vy", "vz"],
)
box = dbox.box

model = zobov.ZobovVF(
    density_threshold=0.1, box_size=1000, workdir=workdir_path
)
void = model.find(dbox)
# joblib.dump(voids, workdir / f"popcorn_{densth}.pkl")
