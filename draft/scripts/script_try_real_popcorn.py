from voidfindertk import io, popcorn, svf_popcorn, core
import os
import pathlib
import matplotlib.pyplot as plt
import joblib

dataset_path = pathlib.Path("./datasets")
workdir_path = pathlib.Path(os.path.abspath("./run_real_popcorn"))

dbox = io.read_table(
    dataset_path / "halos_ascii_1000_1024_npmin_10_z0.51.dat",
    names=["m", "x", "y", "z", "vx", "vy", "vz"],
)


# void09 = model_09.find(databox=dbox)
# void08 = model_08.find(databox=dbox)
# void07 = model_07.find(databox=dbox)

# solo hacer ultima parte
## Estudio radio de shot noise
# model07
model07 = svf_popcorn.SVFPopCorn(
    workdir=workdir_path,
    boxsize=1000,
    densth=-0.7,
    minradius=10,
    maxradius=100,
)
tracers07, centers07, extra07 = model07.build_voids(
    {
        "run_work_dir": model07.workdir
        / "tmplxhn4wcw2024-08-22T01:08:51.465883+00:00",
        "box": dbox.box,
    }
)
voids07 = core.Voids(
    method="svf_popcorn",
    box=dbox.box,
    tracers_in_voids_=tracers07,
    centers_=centers07,
    extra_=extra07,
)
joblib.dump(voids07, "svf_popcorn_model_07.pkl")
plot07 = core.VoidPlotter(voids=voids07)
vsf_kws = {
    "delta": -0.7,
    "n_neighbors": 600,
    "n_step1": 1,
    "n_step2": 3,
    "scale_1_num_samples": 8,
    "scale_2_num_samples": 1,
}
ax7 = plot07.void_size_function(vsf_kws=vsf_kws)

# fig2 = plt.figure()
# fig2.axes.append(ax7)

# fig2.savefig("vsf07.jpg")
ax7.figure.savefig("svf_popcorn_07.jpg")
fig, ax = plt.subplots()
#### FIn estudio


# model08
# model08 = svf_popcorn.SVFPopCorn(
#     workdir=workdir_path,
#     boxsize=1000,
#     densth=-0.8,
#     minradius = 10,
#     maxradius = 100,
# )
# tracers08,centers08,extra08 = model08.build_voids({
#     "run_work_dir": model08.workdir / "tmpbxqny8ro2024-09-14T23:59:51.396390+00:00",
#     "box": dbox.box})
# voids08 = core.Voids(
#     method='svf_popcorn',
#     box = dbox.box,
#     tracers_in_voids_ = tracers08,
#     centers_ = centers08,
#     extra_ = extra08
# )
# joblib.dump(voids08, "svf_popcorn_model_08.pkl")
# plot08 = core.VoidPlotter(
#     voids=voids08
# )
# vsf_kws={'delta':-0.8,'n_neighbors':400,'n_step1':1,'n_step2':3,'scale_1_num_samples':8,'scale_2_num_samples':1}
# ax8 = plot08.void_size_function(vsf_kws=vsf_kws)
# ax8.figure.savefig("vsf_popcorn_08.jpg")

# fig2 = plt.figure()
# fig2.axes.append(ax8)
# plt.savefig("vsf08.jpg")
# fig2.savefig("vsf08.jpg")
#### FIn estudio


# model09 for svf popcron
# model09 = svf_popcorn.SVFPopCorn(
#     workdir=workdir_path,
#     boxsize=1000,
#     densth=-0.9,
#     minradius = 10,
#     maxradius = 100,
# )
# tracers,centers,extra = model09.build_voids({
#     "run_work_dir": model09.workdir / "tmpipk1m12i2024-09-14T23:29:54.863443+00:00",
#     "box": dbox.box})
# voids = core.Voids(
#     method='svf_popcorn',
#     box = dbox.box,
#     tracers_in_voids_ = tracers,
#     centers_ = centers,
#     extra_ = extra
# )
# joblib.dump(voids, "svf_popcorn_model_09.pkl")
# plot09 = core.VoidPlotter(voids=voids)
# vsf_kws={'delta':-0.9,'n_neighbors':400,'n_step1':1,'n_step2':3,'scale_1_num_samples':8,'scale_2_num_samples':1}
# ax = plot09.void_size_function(vsf_kws=vsf_kws)
# ax.figure.savefig("svf_popcorn_09.jpg")
# fig,ax = plt.subplots()
# plot = core.VoidPlotter(
#     voids=voids
# )
# ax = plot.void_size_function()
# ax.figure.savefig("vsf09.jpg")
# #### Fin estudio

# model_09 = popcorn.PopCorn(
#     workdir=workdir_path,
#     boxsize=1000,
#     densth=-0.9,
#     minradius = 10,
#     maxradius = 100,
#     shot_noise_threshold = 19
# )


# tracers09,centers09,extra09 = model_09.build_voids({
#     "run_work_dir": model_09.workdir / "tmpipk1m12i2024-09-14T23:29:54.863443+00:00",
#     "box": dbox.box})

# voids = core.Voids(
#     method='real_popcorn',
#     box = dbox.box,
#     tracers_in_voids_ = tracers09,
#     centers_ = centers09,
#     extra_ = extra09
# )
# joblib.dump(voids, "real_popcorn_model_09.pkl")
# # get effetive radius
# # effective_rad = voids.effective_radius()
# plot09 = core.VoidPlotter(voids=voids)
# vsf_kws={'delta':-0.9,'n_neighbors':400,'n_step1':1,'n_step2':3,'scale_1_num_samples':8,'scale_2_num_samples':1}
# ax = plot09.void_size_function(vsf_kws=vsf_kws)
# ax.figure.savefig("real_popcorn_999.jpg")
# fig,ax = plt.subplots()


# model_08 = popcorn.PopCorn(
#     workdir=workdir_path,
#     boxsize=1000,
#     densth=-0.8,
#     minradius = 10,
#     maxradius = 100,
#     shot_noise_threshold = 20
# )

# tracers08,centers08,extra08 = model_08.build_voids({
#     "run_work_dir": model_08.workdir / "tmpbxqny8ro2024-09-14T23:59:51.396390+00:00",
#     "box": dbox.box})

# voids_real_08 = core.Voids(
#     method='real_popcorn',
#     box = dbox.box,
#     tracers_in_voids_ = tracers08,
#     centers_ = centers08,
#     extra_ = extra08
# )
# joblib.dump(voids_real_08, "real_popcorn_model_08.pkl")

# plot08 = core.VoidPlotter(voids=voids_real_08)
# vsf_kws={'delta':-0.8,'n_neighbors':400,'n_step1':1,'n_step2':3,'scale_1_num_samples':8,'scale_2_num_samples':1}
# ax = plot08.void_size_function(vsf_kws=vsf_kws)
# ax.figure.savefig("real_popcorn_08.jpg")
# fig,ax = plt.subplots()


model_07 = popcorn.PopCorn(
    workdir=workdir_path,
    boxsize=1000,
    densth=-0.7,
    minradius=10,
    maxradius=100,
    shot_noise_threshold=19.7,
)

tracers07, centers07, extra07 = model_07.build_voids(
    {
        "run_work_dir": model_07.workdir
        / "tmplxhn4wcw2024-08-22T01:08:51.465883+00:00",
        "box": dbox.box,
    }
)

voids07 = core.Voids(
    method="real_popcorn",
    box=dbox.box,
    tracers_in_voids_=tracers07,
    centers_=centers07,
    extra_=extra07,
)
joblib.dump(voids07, "real_popcorn_model_07.pkl")
plot08 = core.VoidPlotter(voids=voids07)
vsf_kws = {
    "delta": -0.7,
    "n_neighbors": 600,
    "n_step1": 1,
    "n_step2": 3,
    "scale_1_num_samples": 8,
    "scale_2_num_samples": 1,
}
ax = plot08.void_size_function(vsf_kws=vsf_kws)
ax.figure.savefig("real_popcorn_07.jpg")
fig, ax = plt.subplots()
