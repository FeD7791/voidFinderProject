import struct
import numpy as np

def read_volume_file(*, filename):
    with open(filename, 'rb') as f:
        # Read number of tracers
        number_voids = struct.unpack('i', f.read(4))[0]
        # Read volumes
        volumes = np.zeros(number_voids, dtype=np.float32)
        for i in range(number_voids):
            volume = struct.unpack('d', f.read(8))[0]
            volumes[i] = np.float32(volume)
    return volumes

def _get_tracers_xyz(*, box):
    tracer_x = box.x.value
    tracer_y = box.y.value
    tracer_z = box.z.value
    xyz_arr = np.stack([tracer_x,tracer_y,tracer_z],axis=1)
    return xyz_arr

def _calculate_barycentre(*, tracers_xyz,tracers, tracer_volumes):
    tracer_volumes = tracer_volumes[tracers]
    arr = tracers_xyz[tracers]
    center = np.average(arr,weights=tracer_volumes,axis=0)
    return center

def _get_volumes_from_properties(*,void_properties):
    void_volumes = []
    for void in void_properties:
        void_volumes.append(void.void_vol.value)
    void_volumes = np.array(void_volumes, dtype=np.float32)
    return void_volumes

def _calculate_r_eff(*, void_volumes):
    r_eff = ((3/(4*np.pi))*void_volumes)**(1/3)
    return r_eff

def get_center_and_radii(
        *,
        void_properties,
        tracer_volumes,
        particle_by_voids,
        box
        ):
    # Get the void_volumes
    void_volumes = _get_volumes_from_properties(
        void_properties=void_properties
        )
    # Get r_eff
    void_r_eff = _calculate_r_eff(void_volumes=void_volumes)

    # Get tracers xyz coords
    tracers_xyz = _get_tracers_xyz(box=box)

    # Get centers
    centers = []
    for tracers_in_void in particle_by_voids:
        center = _calculate_barycentre(
            tracers_xyz=tracers_xyz,
            tracers=tracers_in_void,
            tracer_volumes=tracer_volumes)
        centers.append(center)
    return void_r_eff, centers

