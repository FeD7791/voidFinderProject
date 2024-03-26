#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Rava Jorge Federico, Gualpa Sebastian
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
"""Table reader module."""
import sys
import scipy.sparse
import numpy as np

def read_popcorn_output_to_dictionary(filename):
    """
    Reads the output from a file named "popcorn" and stores the information in a dictionary.

    Parameters
    ----------
    filename (str): The name of the file to read.

    Returns:
    --------
    dict: A dictionary containing the information read from the file.
    
    The dictionary structure is as follows:
    {
        'number_voids': int,
        'voids': {
            'void_1': {
                'void_id': int,
                'number_members': int,
                'void_volume': float,
                'number_subhalos': int,
                'internal_flag': int,
                'spheres': [
                    {
                        'x': float,'y': float,'z': float,
                        'sphere_radius': float
                    },
                    ...
                ],
                'ID_subhalo': [int, ...]
            },
            'void_2': {
                ...
            },
            ...
        }
    }

    """
    datos = {}
    
    with open(filename, 'r') as file:
        # Read the first line to get the number of voids
        cantidad_voids = int(file.readline().strip())
        datos['number_voids'] = cantidad_voids
        datos['voids'] = {}
        
        # Iterate over each void in the file
        for _ in range(cantidad_voids):
            # Read the information of each void
            void_info = file.readline().strip().split(' ')
            void_id = int(void_info[0])
            num_members = int(void_info[1])
            void_volume = float(void_info[2])
            num_subhalos = int(void_info[3])
            internal_flag = int(void_info[4])  # Cambiado a entero
            
            # Create a dictionary for each void
            void_dict = {
                'void_id': void_id,
                'number_members': num_members,
                'void_volume': void_volume,
                'number_subhalos': num_subhalos,
                'internal_flag': internal_flag,
                'spheres': [],
                'ID_subhalo': []
            }
            
            # Read the list of spheres
            for _ in range(num_members):
                esfera_info = file.readline().strip().split()  # Dividir por espacios en blanco
                x = float(esfera_info[0])
                y = float(esfera_info[1])
                z = float(esfera_info[2])
                sphere_radius = float(esfera_info[3])
                void_dict['spheres'].append({'x': x, 'y': y, 'z': z, 'sphere_radius': sphere_radius})
            
            # Read the list of subhalo IDs
            for _ in range(num_subhalos):
                id_subhalo = int(file.readline().strip())
                void_dict['ID_subhalo'].append(id_subhalo)
            
            # Add the void dictionary to the main dictionary
            datos['voids'][f'void_{void_id}'] = void_dict
    
    return datos

def create_sparse_matrix_from_dict(data, num_max_subhalo):
    """
    Crea una matriz sparse de scipy a partir de un diccionario.

    Parameters:
    data (dict): Un diccionario con el formato especificado.

    Returns:
    scipy.sparse.csr_matrix: Una matriz sparse de scipy.
    """
    # Extraer el número total de voids
    number_voids = data['number_voids']
    
    # Crear listas para almacenar los void_id y ID_subhalo
    void_ids = []
    id_subhalos = []
    
    # Iterar sobre cada void en el diccionario
    for void_key in data['voids']:
        void_info = data['voids'][void_key]
        
        # Obtener void_id
        void_id = void_info['void_id']
        void_ids.append(void_id)
        
        # Obtener ID_subhalo
        id_subhalo = void_info['ID_subhalo']

        if not isinstance(id_subhalo, list):
            var_aux = id_subhalos
            id_subhalo = []
            id_subhalo.append(var_aux)
        id_subhalos.append(id_subhalo)

    # Crear una matriz numpy con los valores de void_id y ID_subhalo
    matrix_data = np.zeros((len(id_subhalos), num_max_subhalo), dtype=float)
    #print(matrix_data)

    for i, subhalo_indices in enumerate(id_subhalos):
        matrix_data[i, subhalo_indices] = 1

    #for i, fila in enumerate(matrix_data):
    #    fila_filtrada = [valor for valor in fila if valor != 0]
    #    if fila_filtrada:
    #        print(f"Fila {i + 1}: {fila_filtrada}")


    # Crear una matriz sparse csr_matrix
    sparse_matrix = scipy.sparse.csr_matrix(matrix_data)

    return sparse_matrix