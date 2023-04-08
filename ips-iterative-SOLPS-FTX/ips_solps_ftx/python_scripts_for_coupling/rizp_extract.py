#!/usr/bin/env python
import numpy as np
from collections import defaultdict

def rizp_extract(species):
    # Read ionization potential data from file
    filename = 'ionization_potentials.txt'

    with open(filename, 'r') as file:
        lines = file.readlines()

    # Read the lower diagonal data
    n = len(lines)
    ionization_potentials = np.zeros((n, n))
    for i, line in enumerate(lines):
        line = line.replace('d', 'e')  # Replace 'd' with 'e' for correct exponent notation
        values = list(map(float, line.split()))
        ionization_potentials[i, :i+1] = values

    # Define a mapping from species to atomic numbers
    species_to_atomic_number = dict(zip(
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
         'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
         'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
         'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
         'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La',
         'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
         'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
         'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',
         'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
         'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg',
         'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og', 'ZERO', 'D', 'T'],
        list(range(1, 119)) + [0, 1, 1]))

    # Split the species string and extract atomic numbers
    species_list = species.split('+')
    atomic_numbers = [species_to_atomic_number[x] for x in species_list]


    # Extract the data for the specified species
    extracted_data = []
    for idx in atomic_numbers:
        extracted_row = ionization_potentials[idx - 1, :idx]
        extracted_data.extend([0] + extracted_row.tolist())

    extracted_data = np.array(extracted_data)
    return extracted_data
