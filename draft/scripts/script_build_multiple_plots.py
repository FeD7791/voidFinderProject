from voidfindertk import io, popcorn, svf_popcorn, core
import os
import pathlib
import matplotlib.pyplot as plt
import joblib

svf = joblib.load("svf_-0.8.pkl")
zobov = joblib.load("zobov_model_08.pkl")
dive = joblib.load("dive_model_t_0.2_rat_1.5.pkl")
real_popcorn = joblib.load("popcorn_-0.8.pkl")

###
svf = joblib.load("svf_-0.9.pkl")
zobov = joblib.load("zobov_model_0.1.pkl")
dive = joblib.load("dive_model_t_0.1_rat_1.5.pkl")
real_popcorn = joblib.load("popcorn_-0.9.pkl")

plot = core.VoidPlotter(voids=zobov)
vsf_kws={'delta':-0.9,'n_neighbors':200,'n_step1':1,'n_step2':5,'scale_1_num_samples':4,'scale_2_num_samples':2}
ax = plot.void_size_function(vsf_kws=vsf_kws,label="ZOBOV")
ax.figure.savefig("models.jpg")

plot = core.VoidPlotter(voids=svf)
vsf_kws={'delta':-0.9,'n_neighbors':200,'n_step1':1,'n_step2':5,'scale_1_num_samples':4,'scale_2_num_samples':2}
ax = plot.void_size_function(vsf_kws=vsf_kws, label="Spherical")
ax.figure.savefig("models.jpg")

plot = core.VoidPlotter(voids=real_popcorn)
vsf_kws={
    'delta':-0.9,'n_neighbors':200,
    'n_step1':1,'n_step2':5,'scale_1_num_samples':4,
    'scale_2_num_samples':2,'r_eff_from_extra':True,
    'r_eff_from_extra_delta':-0.9
    }
ax = plot.void_size_function(vsf_kws=vsf_kws, label="Popcorn")
ax.figure.savefig("models.jpg")

plot = core.VoidPlotter(voids=dive)
vsf_kws={'delta':-0.9,'n_neighbors':200,'n_step1':1,'n_step2':3,'scale_1_num_samples':5,'scale_2_num_samples':1}
ax = plot.void_size_function(vsf_kws=vsf_kws,label="Dive")
ax.figure.savefig("models.jpg")
fig,ax = plt.subplots()