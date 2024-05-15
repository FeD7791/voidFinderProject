import numpy as np
import scipy
def find_radius(box,voids,sparse):
    radius = []
    tracers = [[box.x.value[i],box.y.value[i],box.z.value[i]] for i in range(len(box))]
    #DEPENDE DEL VOIDFINDER
    ##CONFIGURAR PARA VARIOS VOID FINDERS
    voids_centers = [
        [
            voids.x_void.value[i],
            voids.y_void.value[i],
            voids.z_void.value[i]] for i in len(voids)]
    ##########################################################
    sparse = sparse.toarray() #Transform sparse matrix to array [M,N]
    for i in range(sparse.shape[1]):
        c = sparse[:,i][:, None] * tracers # filter the tracers in voids
        c = c[np.any(c, axis=1)] #Removes [0,0,0] rows
        c = np.vstack([c,voids_centers[i]]) # Adding void center to points
        distances = scipy.spatial.distance.pdist(c) #distances between points
        radius.append(distances.sort(-1)[-1]) #appending biggest distance to radius as candidate
    return radius