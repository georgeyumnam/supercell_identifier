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


def sup_zero_conv(sup_nx, sup_ny, sup_nz):
    sup = [sup_nx, sup_ny, sup_nz]
    for i in range(3):
        if sup[i] == 0:
	    sup[i] = 1
	    for j in range(3):
	        if j != i :
		    sup[j] = sup[j] - 1
    return sup[0], sup[1], sup[2]


def final_sup_check(sup_nx, sup_ny, sup_nz, aaa, bbb, ccc, num_atoms_in_unit_cell, given):
    supp = [sup_nx, sup_ny, sup_nz]
    n1 = supp[0] * aaa
    n2 = supp[1] * bbb
    n3 = supp[2] * ccc
    sup = np.array([[sup_nx, n1, 0], [sup_ny, n2, 1], [sup_nz, n3, 2]])
    mmm = sup[sup[:, 1].argsort()]
    min_n = mmm[0][1]
    max_n = mmm[2][1]
    num_at_sup = num_atoms_in_unit_cell * sup_nx * sup_ny * sup_nz
    if num_at_sup < given :
        if mmm[2][1] / mmm[0][1] > 1.8 :
	    mmm[0][0] = mmm[0][0] + 1
	if mmm[2][1] / mmm[1][1] > 1.8 :
	    mmm[1][0] = mmm[1][0] + 1 
    mmm = np.array(mmm)
    nnn = mmm[mmm[:, 2].argsort()]
    supx, supy, supz = nnn[0][0], nnn[1][0], nnn[2][0]
    return supx, supy, supz
    




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
      if check < 1 :
        if num_at_sup > given_hard :
            sup_nx = sup_nx - 1
            sup_ny = sup_ny - 1
            sup_nz = sup_nz - 1
	    sup_nx, sup_ny, sup_nz = sup_zero_conv(sup_nx, sup_ny, sup_nz)
            num_at_sup = num_atoms_in_unit_cell * sup_nx * sup_ny * sup_nz
            check = given / num_at_sup
        else :
            nx, ny, nz = sup_nx, sup_ny, sup_nz
            checker = "accept"
      elif check > 1 :
        sup_nx = sup_nx + 1
        sup_ny = sup_ny + 1
        sup_nz = sup_nz + 1
        num_at_sup = num_atoms_in_unit_cell * sup_nx * sup_ny * sup_nz
        check = given / num_at_sup
        if check < 1 :
          if num_at_sup > given_hard :
            sup_nx = sup_nx - 1
            sup_ny = sup_ny - 1
            sup_nz = sup_nz - 1
	    sup_nx, sup_ny, sup_nz = sup_zero_conv(sup_nx, sup_ny, sup_nz)
            num_at_sup = num_atoms_in_unit_cell * sup_nx * sup_ny * sup_nz
            check = given / num_at_sup
            checker = "accept"        
    

    sup_nx, sup_ny, sup_nz = final_sup_check(sup_nx, sup_ny, sup_nz, aaa, bbb, ccc, num_atoms_in_unit_cell, given)
    num_at_sup = num_atoms_in_unit_cell * sup_nx * sup_ny * sup_nz

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
