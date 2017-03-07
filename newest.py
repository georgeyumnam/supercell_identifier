from ase.geometry import cell as cl
from spglib import find_primitive, standardize_cell
from phonopy.interface import vasp as v
from phonopy.structure.cells import *
import numpy as np


# From the input POSCAR file
poscar = 'POSCAR'
cell = v.read_vasp(poscar)
prim_lat, prim_pos, prim_num = find_primitive(cell, symprec=1e-1)
stnd_lat, stnd_pos, stnd_num = standardize_cell(cell, symprec=1e-1)


def number_of_atoms_in_cell(position):
    """
    Count the number of atoms for given cell -- positions
    """
    num_atoms = np.shape(position)[0]
    return num_atoms


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
    if np.absolute(var) < 0.000001:
        return 0.0
    else:
        return var


def supercell_lattice(lattice, nx, ny, nz):
    """
    This function will create a supercell lattice of the given unit cell
    lattice with a dimension of [nx, ny, nz]
    """
    a, b, c, alpha, beta, gamma = cl.cell_to_cellpar(lattice)
    sup_A = a * nx
    sup_B = b * ny
    sup_C = c * nz
    SUP_cellpar = sup_A, sup_B, sup_C, alpha, beta, gamma
    SUP_cell = cl.cellpar_to_cell(SUP_cellpar, ab_normal=(0, 0, 1),
                                  a_direction=(1, 0, 0))
    return SUP_cell


def supercell_lattice_smoothen(checker_lattice):
    """
    This function smoothens the created supercell lattice from numerical errors
    """
    k = checker_lattice
    for i in range(3):
        for j in range(3):
            k[i][j] = zero_checker(checker_lattice[i][j])
    return k

def permutation_definer(num1, num2, num3):
    """
    This function returns the total permutations possible for all 
    given numbers for computing the LCM
    """
    matrix = np.zeros((27, 3))
    index = 0
    for i in range(3):
        x = num1 - 1 + i
        for j in range(3):
            y = num2 - 1 + j
            for k in range(3):
                z = num3 - 1 + k
                matrix[index][:] = [x, y, z]
                index = index + 1
    return matrix


def best_combo_with_least_LCM(num1, num2, num3):
    """
    Determines the best permutations possible with least LCM
    """
    mat = permutation_definer(num1, num2, num3)
    mat = np.array(mat)
    LCM = 1000 # Arbitrary large LCM value for initialization
    for i in range(27):
        LCM_new = func_lcm(mat[i][0], mat[i][1], mat[i][2])
        if LCM_new <= LCM :
            LCM = LCM_new
            req_mat = mat[i]
    a, b, c = req_mat
    return a, b, c, LCM

def main_function(lattice, positions):
    a, b, c, alpha, beta, gamma = cl.cell_to_cellpar(lattice)
    aa = round_off_of_num(a)
    bb = round_off_of_num(b)
    cc = round_off_of_num(c)

    aaa, bbb, ccc, LCM = best_combo_with_least_LCM(aa, bb, cc)

    sup_nx = LCM / aaa
    sup_ny = LCM / bbb
    sup_nz = LCM / ccc

    num_atoms_in_unit_cell = number_of_atoms_in_cell(positions)

    given = 100
    given_hard = 140
    
    num_at_sup = num_atoms_in_unit_cell * sup_nx * sup_ny * sup_nz
    check = given / num_at_sup
   
    checker = "unacceptable"
 
    while checker != "accept" :
      print 1
      if check < 1 :
        print 2
        if num_at_sup > given_hard :
            print 3
            sup_nx = sup_nx - 1
            sup_ny = sup_ny - 1
            sup_nz = sup_nz - 1
            num_at_sup = num_atoms_in_unit_cell * sup_nx * sup_ny * sup_nz
            check = given / num_at_sup
        else :
            print 4
            nx, ny, nz = sup_nx, sup_ny, sup_nz
            checker = "accept"
      elif check > 1 :
        print 5
        sup_nx = sup_nx + 1
        sup_ny = sup_ny + 1
        sup_nz = sup_nz + 1
        num_at_sup = num_atoms_in_unit_cell * sup_nx * sup_ny * sup_nz
        check = given / num_at_sup
        if check < 1 :
          print 6
          if num_at_sup > given_hard :
            print 7
            sup_nx = sup_nx - 1
            sup_ny = sup_ny - 1
            sup_nz = sup_nz - 1
            num_at_sup = num_atoms_in_unit_cell * sup_nx * sup_ny * sup_nz
            check = given / num_at_sup
            checker = "accept"        

    SUP_LATT = supercell_lattice(lattice, sup_nx, sup_ny, sup_nz)
    SMOOTH_SUP_LATT = supercell_lattice_smoothen(SUP_LATT)
    print " "
    print "The required supercell size for this case is : ", \
        sup_nx, "x", sup_ny, "x", sup_nz
    print "The number of atoms in the unitcell is       : ", \
        num_atoms_in_unit_cell
    print "The number of atoms with this supercell is   : ", \
        num_at_sup
    print " "
    print " The unit cell lattice is : "
    print lattice
    print " "
    print " The supercell lattice is : "
    print SMOOTH_SUP_LATT
    print " "

main_function(stnd_lat, stnd_pos)
