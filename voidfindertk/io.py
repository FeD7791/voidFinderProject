#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Rava Jorge Federico, Gualpa Sebastian
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
"""Table reader module."""
import pandas as pd

from . import box

# from .popcornvf import popcorn_tools
#
# from .sphericalvf import spherical_tools



def read_table(
    path_or_buffer, **kwargs
):  # names = ['x','y','z','vx','vy','vz','m'] to put column names
    """Input reader.

    Read a table from a file or buffer and returns a `box.Box` object.

    The table must contain 7 columns, with the following format:

    x y z vx vy vz m


    where:

    * `x`, `y`, and `z` are the coordinates of the box center
    * `vx`, `vy`, and `vz` are the velocities of the box center
    * `m` is the mass of the box

    Parameters
    ----------
        path_or_buffer: The path to the file or buffer containing the table.
        **kwargs: Keyword arguments to be passed to `pandas.read_csv()`.

    Returns
    -------
        A `box.Box` object containing the data from the table.

    """
    kwargs.setdefault("sep", r"\s+")
    kwargs.setdefault("usecols", [0, 1, 2, 3, 4, 5, 6])
    kwargs.setdefault("names", ["x", "y", "z", "vx", "vy", "vz", "m"])
    data = pd.read_csv(path_or_buffer, **kwargs, header=None)
    col_number = len(data.columns)

    if col_number != 7:
        raise ValueError(
            "There are not enough columns to create the coordinates of a box."
            f"Found {col_number} expected 7"
        )

    check_values = data.notnull().values.all()
    if not check_values:
        raise TypeError(
            f"There are: {data.isnull().sum().sum()}\
                  null or missing values"
        )
    #Clean duplicates
    data.drop_duplicates(ignore_index=True, inplace=True)

    the_box = box.Box(
        x=data.loc[:, "x"],
        y=data.loc[:, "y"],
        z=data.loc[:, "z"],
        vx=data.loc[:, "vx"],
        vy=data.loc[:, "vy"],
        vz=data.loc[:, "vz"],
        m=data.loc[:, "m"],
    )
    the_box = box.DataBox(the_box)
    return the_box


# def read_popcorn_output(filename):
#    """
#    Reads the output from a file named "popcorn" and stores the information in a dictionary.
#
#    Parameters
#    ----------
#    filename (str): The name of the file to read.
#
#    Returns:
#    --------
#    dict: A dictionary containing the information read from the file.
#
#    The dictionary structure is as follows:
#    {
#        'number_voids': int,
#        'voids': {
#            'void_1': {
#                'void_id': int,
#                'number_members': int,
#                'void_volume': float,
#                'number_subhalos': int,
#                'internal_flag': int,
#                'spheres': [
#                    {
#                        'x': float,'y': float,'z': float,
#                        'sphere_radius': float
#                    },
#                    ...
#                ],
#                'ID_subhalo': [int, ...]
#            },
#            'void_2': {
#                ...
#            },
#            ...
#        }
#    }
#
#    """
#    datos = {}
#
#    with open(filename, 'r') as file:
#        # Read the first line to get the number of voids
#        cantidad_voids = int(file.readline().strip())
#        datos['number_voids'] = cantidad_voids
#        datos['voids'] = {}
#
#        # Iterate over each void in the file
#        for _ in range(cantidad_voids):
#            # Read the information of each void
#            void_info = file.readline().strip().split(' ')
#            void_id = int(void_info[0])
#            num_members = int(void_info[1])
#            void_volume = float(void_info[2])
#            num_subhalos = int(void_info[3])
#            internal_flag = int(void_info[4])  # Cambiado a entero
#
#            # Create a dictionary for each void
#            void_dict = {
#                'void_id': void_id,
#                'number_members': num_members,
#                'void_volume': void_volume,
#                'number_subhalos': num_subhalos,
#                'internal_flag': internal_flag,
#                'spheres': [],
#                'ID_subhalo': []
#            }
#
#            # Read the list of spheres
#            for _ in range(num_members):
#                esfera_info = file.readline().strip().split()  # Dividir por espacios en blanco
#                x = float(esfera_info[0])
#                y = float(esfera_info[1])
#                z = float(esfera_info[2])
#                sphere_radius = float(esfera_info[3])
#                void_dict['spheres'].append({'x': x, 'y': y, 'z': z, 'sphere_radius': sphere_radius})
#
#            # Read the list of subhalo IDs
#            for _ in range(num_subhalos):
#                id_subhalo = int(file.readline().strip())
#                void_dict['ID_subhalo'].append(id_subhalo)
#
#            # Add the void dictionary to the main dictionary
#            datos['voids'][f'void_{void_id}'] = void_dict
#
#            #Define df_voids
#            voids_list = list(datos['voids'])
#            voids_all = []
#            for i in voids_list:
#                voids_all.append(datos['voids'][f'{i}'])
#            df_voids = pd.DataFrame.from_dict(voids_all)
#            #Popcorn object
#            popcorn_voids = popcorn_tools.PopCornVoids(**df_voids.to_dict(orient='list'))
#
#    return popcorn_voids

# def read_spherical_output(filename):
#    output = pd.read_csv(filename, sep='\s+')
#    output.columns = [
#        'rad',
#        'x_void',
#        'y_void',
#        'z_void',
#        'vel_x_void',
#        'vel_y_void',
#        'vel_z_void',
#        'delta',
#        'dtype',
#        'poisson',
#        'nran'
#        ]
#    spherical = spherical_tools.SphericalVoids(**output.to_dict(orient='list'))
#    return spherical
