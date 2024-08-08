import numpy as np

def compute_void_size_function(delta, rsph1, n_tracers, box_size, final_n_bins):
    """
    Computes the void size function

    Parameters
    ----------
        delta : float
        Integrated density Contrast 
        rsph1 : array
        Array of radius of each void in Mpc
        n_tracers : int
        Number of tracers in the whole box.
        box_size : int
        Length in Mpc of the box
    """
    pi = np.pi
    # Volume of the simulation
    # vol = 1000**3
    vol = box_size**3

    # Mean density of tracers
    rhomed = (n_tracers / vol)

    # N sequence
    N = np.concatenate([np.arange(1, 6), np.arange(6, 11, 2), np.arange(12, 53, 10)])

    # Scaling calculation
    scl = np.log10((3 / (4 * pi) * N / rhomed / (1 + delta))**(1/3))

    mxlg = np.log10(max(rsph1))
    br = np.concatenate([scl[:-1], np.linspace(max(scl), mxlg, final_n_bins)])
    
    # Histogram calculation
    h, bin_edges = np.histogram(np.log10(rsph1), bins=br)

    # Compute dlogr
    dlogr = np.diff(br)
    mids = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    density = h / dlogr / vol

    return 10.0**mids, density


def calculate_void_size_function(radii, integrated_density, mean_density, bins):
    """
    Calculate the void size function from given void radii and densities.
    
    Parameters:
    - radii: Array of radii of the voids (in the same units).
    - integrated_density: Array of integrated densities of the voids.
    - mean_density: Mean density of tracers.
    - bins: Number of bins or edges for histogram.
    
    Returns:
    - bin_centers: Center of each bin.
    - void_counts: Number of voids in each bin.
    """
    # Calculate void volume from radii
    volumes = (4/3) * np.pi * radii**3
    
    # Calculate the density of each void
    void_density = integrated_density / volumes
    
    # Calculate the effective volume for each void
    effective_volumes = volumes * (mean_density / void_density)
    
    # Histogram the void sizes
    bin_edges = np.linspace(min(radii), max(radii), bins + 1)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    
    # Count voids in each bin
    void_counts, _ = np.histogram(radii, bins=bin_edges)
    
    # Normalize counts by bin width and total volume
    bin_widths = np.diff(bin_edges)
    total_volume = np.sum(volumes)  # Total volume for normalization
    void_counts = void_counts / (bin_widths * total_volume)
    
    return bin_centers, void_counts
