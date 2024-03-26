from io_popcorn import read_popcorn_output_to_dictionary
from io_popcorn import create_sparse_matrix_from_dict
import numpy as np

filename='../../../popcorn_202403/pocorn4/popcorn/MWE/clean_popvds.dat'
datos = read_popcorn_output_to_dictionary(filename)
sparse_matrix = create_sparse_matrix_from_dict(datos,700070)
print(sparse_matrix)
