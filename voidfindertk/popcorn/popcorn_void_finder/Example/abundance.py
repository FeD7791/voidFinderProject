import numpy as np
import matplotlib.pyplot as plt

from read_voids import read_sph, read_pop

# Reading the spherical void catalogue:
sphvoids = read_sph("sphvoids.dat")

# Reading the popcorn void catalogue:
popvoids = read_pop("popvoids.dat")

# Calculating radii histograms:
hs, edges_s = np.histogram(sphvoids['r'])
hp, edges_p = np.histogram(popvoids['reff'])

# Midpoints of bins:
mids_s = 0.5 * (edges_s[1:] + edges_s[:-1])
mids_p = 0.5 * (edges_p[1:] + edges_p[:-1])

# Plotting radii distribution:
plt.loglog(mids_s, hs, label='spherical', linewidth=2, color='blue')
plt.loglog(mids_p, hp, label='popcorn', linewidth=2, color='red')

plt.xlabel(r'$R_v$ [Mpc/h]')
plt.ylabel('Void counts')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
