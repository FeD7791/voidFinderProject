import numpy as np

def compute_void_size_function(Delta, rsph1, n_tracers, box_size):
    pi = np.pi
    # Volume of the simulation
    # vol = 1000**3
    vol = box_size**3

    # Mean density of tracers
    rhomed = (n_tracers / vol)

    # N sequence
    N = np.concatenate([np.arange(1, 6), np.arange(6, 11, 2), np.arange(12, 53, 10)])

    # Scaling calculation
    scl = np.log10((3 / (4 * pi) * N / rhomed / (1 + Delta))**(1/3))

    mxlg = np.log10(max(rsph1))
    br = np.concatenate([scl[:-1], np.linspace(max(scl), mxlg, 50)])
    
    # Histogram calculation
    h, bin_edges = np.histogram(np.log10(rsph1), bins=br)

    # Compute dlogr
    dlogr = np.diff(br)
    mids = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    density = h / dlogr / vol

    return 10.0**mids, density

plt.figure(figsize=(8, 6))
plt.plot(x,y)
plt.xscale('log')
plt.yscale('log')
plt.show()