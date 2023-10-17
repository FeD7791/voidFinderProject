import numpy as np
import pandas as pd
from numpy import linalg as LA


def to_np_arr(element):
    """
    Transforms the input element in a numpy array
    """
    return np.array(element)


def distance_metric(x, y, z):
    """
    Calculates the Euclidean norm for an array of vectors.
    For each element of the vectors an Euclidean Norm is calculated
    using linalg norm.
    d[i] = linag.norm(x[i],y[i],z[i])
    Input vectors must have the same size.

    Parameters
    ----------
    x : array(float)
    y : array(float)
    z : array(float)

    Returns
    -------
    array(float)
        Array of Euclidean norms calculated from the input elements
    """
    d = np.array([])
    for i in range(x.size):
        d = np.append(d, LA.norm([x[i], y[i], z[i]]))
    return d


def slicer(x, y, z, n, cols):
    """
    Slices a set of points in a n-grid

    Parameters
    ----------
    x : array(float)
    y : array(float)
    z : array(float)
    n : scalar(int)
    cols : list(str)

    Returns
    -------
    list(dataframe)
        List of dataframes with the set of points inside of a grid
    """
    (x_min, x_max) = (x.min(), x.max())
    (y_min, y_max) = (y.min(), y.max())
    (z_min, z_max) = (z.min(), z.max())

    step_x = (x_max - x_min) / n
    step_y = (y_max - y_min) / n
    step_z = (z_max - z_min) / n

    xx = np.arange(x_min, x_max, step_x)
    yy = np.arange(y_min, y_max, step_y)
    zz = np.arange(z_min, z_max, step_z)

    grid_d = distance_metric(xx, yy, zz)
    arr = np.vstack((x, y, z, distance_metric(x, y, z)))
    arr = pd.DataFrame(arr.transpose())
    arr.columns = cols

    arreglo = []
    for i in range(len(grid_d)):
        if i == 0:
            arrayy = arr.loc[(arr[cols[3]]) < grid_d[i + 1]]
            arreglo.append(arrayy)
        elif i < len(grid_d) - 1:
            arrayy = arr.loc[
                (arr[cols[3]] > grid_d[i]) & (arr[cols[3]] < grid_d[i + 1])
            ]
            arreglo.append(arrayy)
        else:
            arrayy = arr.loc[(arr[cols[3]]) > grid_d[len(grid_d) - 1]]
            arreglo.append(arrayy)
    return arreglo
