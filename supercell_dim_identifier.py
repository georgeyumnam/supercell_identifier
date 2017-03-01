import sys
import os.path
from ase.geometry import cell as cl
from spglib import find_primitive
from phonopy.interface import vasp as v
from phonopy.structure.cells import *
import numpy as np


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


def round_off_of_num(x):
    n = int(x)
    return n if n <= x < (n + 0.5) else n+1


def func_gcd(*numbers):
    from fractions import gcd
    return reduce(gcd, numbers)


def func_lcm(*numbers):
    def lcm(a, b):
        return (a * b) // func_gcd(a, b)
    return reduce(lcm, numbers, 1)


def zero_checker(var):
    if var < 0.001:
        return 0
    else:
        return var


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
    SUP_cell = cl.cellpar_to_cell(SUP_cellpar, ab_normal=(0, 0, 1),
                                  a_direction=None)
    return SUP_cell


def supercell_lattice_smoothen(checker_lattice):
    """
    This function smoothens the created supercell lattice from numerical errors
    """
    k = np.zeros(np.shape(checker_lattice))
    for i in range(3):
        for j in range(3):
            k[i][j] = zero_checker(checker_lattice[i][j])
    return k


def main_function(lattice, positions):
    a, b, c, alpha, beta, gamma = cl.cell_to_cellpar(lattice)
    aa = round_off_of_num(a)
    bb = round_off_of_num(b)
    cc = round_off_of_num(c)

    LCM = func_lcm(aa, bb, cc)
    sup_nx_i = LCM / aa
    sup_ny_i = LCM / bb
    sup_nz_i = LCM / cc

    num_atoms_in_unit_cell = number_of_atoms_in_cell(positions)
    num_atoms_in_sup_cell = sup_nx_i * sup_ny_i * sup_nz_i \
        * num_atoms_in_unit_cell
    sup_nx = sup_nx_i
    sup_ny = sup_ny_i
    sup_nz = sup_nz_i
    mult = 1
    while num_atoms_in_sup_cell < 105:
        if num_atoms_in_sup_cell < 90:
            mult = mult + 1
            sup_nx = mult * sup_nx
            sup_ny = mult * sup_ny
            sup_nz = mult * sup_nz
            num_atoms_in_sup_cell = sup_nx * sup_ny * sup_nz \
                * num_atoms_in_unit_cell
    if num_atoms_in_sup_cell > 113:
        mult = 1
        MAX_N = max(sup_nx, sup_ny, sup_nz)
        MIN_N = min(sup_nx, sup_ny, sup_nz)
        new_sup_nx = sup_nx / MIN_N
        new_sup_ny = sup_ny / MIN_N
        new_sup_nz = sup_nz / MIN_N
        num_atoms_in_sup_cell = new_sup_nx * new_sup_ny * new_sup_nz \
            * num_atoms_in_unit_cell

        if num_atoms_in_sup_cell > 113:
            sup_nx, sup_ny, sup_nz = new_sup_nx, new_sup_ny, new_sup_nz
        elif num_atoms_in_sup_cell <= 113:
            if num_atoms_in_sup_cell > 88:
                sup_nx, sup_ny, sup_nz = new_sup_nx, new_sup_ny, new_sup_nz
            else:
                while num_atoms_in_sup_cell < 88:
                    mult = mult + 1
                    sup_nx, sup_ny, sup_nz = new_sup_nx * mult, new_sup_ny \
                        * mult, new_sup_nz * mult
                    num_atoms_in_sup_cell = sup_nx * sup_ny * sup_nz \
                        * num_atoms_in_unit_cell
        num_atoms_in_sup_cell = sup_nx * sup_ny * sup_nz \
            * num_atoms_in_unit_cell

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
