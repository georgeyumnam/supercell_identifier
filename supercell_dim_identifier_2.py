"""
__author__ = George Yumnam
Code written to choose the supercell for different materials with different geometries.
"""

import numpy as np
from ase.geometry import cell as cl
from spglib import find_primitive
from phonopy.interface import vasp as v


# From the input POSCAR file
poscar = 'POSCAR'
cell = v.read_vasp(poscar)
prim_lat, prim_pos, prim_num = find_primitive(cell, symprec=1e-5)


def number_of_atoms_in_cell(position):
    """
    Count the number of atoms in unitcell
    """
    num_atoms_std = np.shape(position)[0]
    return num_atoms_std


def zero_checker(var):
    if var < 0.001:
        return 0
    else:
        return var


def func_gcd(*numbers):
    from fractions import gcd
    return reduce(gcd, numbers)


def func_lcm(*numbers):
    def lcm(a, b):
        return (a * b) // func_gcd(a, b)
    return reduce(lcm, numbers, 1)


def supercell_lattice(lattice, nx, ny, nz):
    """
    This function will create a supercell lattice of the given unit cell
    lattice with a dimension of [nx, ny, nz]
    """
    a, b, c, alpha, beta, gamma = cl.cell_to_cellpar(lattice)
    T = [nx, ny, nz]
    sup_A = round((a * T[0]), 5)
    sup_B = round((b * T[1]), 5)
    sup_C = round((c * T[2]), 5)
    SUP_cellpar = sup_A, sup_B, sup_C, alpha, beta, gamma
    SUP_cell = cl.cellpar_to_cell(SUP_cellpar,
                                  ab_normal=(0, 0, 1), a_direction=None)
    print " "
    print "Unitcell  -- a, b, c ==", a, " - ", b, " - ", c
    print "Supercell -- a, b, c ==", sup_A, " - ", sup_B, " - ", sup_C
    return SUP_cell


def supercell_lattice_smoothen(checker_lattice):
    """
    This function smoothens the created supercell lattice from numerical error
    """
    k = np.zeros(np.shape(checker_lattice))
    for i in range(3):
        for j in range(3):
            k[i][j] = zero_checker(checker_lattice[i][j])
    return k


def matrix_convertor(array):
    """
    This function caclulates the supercell lattice value on the go
    """
    array[0][2] = array[0][0] * array[0][1]
    array[1][2] = array[1][0] * array[1][1]
    array[2][2] = array[2][0] * array[2][1]
    # The converted array
    return array


def round_off_of_num(x):
    n = int(x)
    if n <= x < (n + 0.5):
        return n
    else:
        return (n+1)


def main_function(lattice, positions):
    a, b, c, alpha, beta, gamma = cl.cell_to_cellpar(lattice)
    num_atoms_in_unit_cell = number_of_atoms_in_cell(positions)

    max_allowed = 120
    min_allowed = 100

    tot_xyz_allowed_max = max_allowed / num_atoms_in_unit_cell
    tot_xyz_allowed_min = min_allowed / num_atoms_in_unit_cell

    sup_nx, sup_ny, sup_nz = 1, 1, 1  # Starting from the unitcell
    sup_lat_a, sup_lat_b, sup_lat_c = a, b, c  # Initially

    abc_nxyz = np.array([[a, sup_nx, sup_lat_a, 0],
                         [b, sup_ny, sup_lat_b, 1],
                         [c, sup_nz, sup_lat_c, 2]])

    while (sup_nx * sup_ny * sup_nz) < tot_xyz_allowed_max:
        if (sup_nx * sup_ny * sup_nz) < tot_xyz_allowed_min:
            abc_nxyz = abc_nxyz[abc_nxyz[:, 2].argsort()]
            abc_nxyz[0][1] = abc_nxyz[0][1] + 1
            abc_nxyz = matrix_convertor(abc_nxyz)
            sup_nx, sup_ny, sup_nz = abc_nxyz[0][1], \
                abc_nxyz[1][1], abc_nxyz[2][1]
            roonum0 = round_off_of_num(abc_nxyz[0][2])
            roonum1 = round_off_of_num(abc_nxyz[1][2])
            roonum2 = round_off_of_num(abc_nxyz[2][2])
            sup_lat_a, sup_lat_b, sup_lat_c = roonum0, roonum1, roonum2
        else:
            sup_nx, sup_ny, sup_nz = 10, 10, 10  # Just to break the while loop

    sup_nx, sup_ny, sup_nz = abc_nxyz[0][1], abc_nxyz[1][1], abc_nxyz[2][1]

    if (sup_nx * sup_ny * sup_nz) < tot_xyz_allowed_min:
        abc_nxyz[0][1] = abc_nxyz[0][1] + 1

    abc_nxyz = abc_nxyz[abc_nxyz[:, 3].argsort()]
    sup_nx, sup_ny, sup_nz = abc_nxyz[0][1], abc_nxyz[1][1], abc_nxyz[2][1]
    num_atoms_in_sup_cell = sup_nx * sup_ny * sup_nz * num_atoms_in_unit_cell
    SUP_LATT = supercell_lattice(lattice, sup_nx, sup_ny, sup_nz)
    SMOOTH_SUP_LATT = supercell_lattice_smoothen(SUP_LATT)

    print " "
    print "The required supercell size for this case is : ", \
        sup_nx, "x", sup_ny, "x", sup_nz
    print "The number of atoms in the unitcell is       : ", \
        num_atoms_in_unit_cell
    print "The number of atoms with this supercell is   : ", \
        num_atoms_in_sup_cell
    print " "
    print " The unit cell lattice is : "
    print lattice
    print " "
    print " The supercell lattice is : "
    print SMOOTH_SUP_LATT
    print " "

main_function(prim_lat, prim_pos)
