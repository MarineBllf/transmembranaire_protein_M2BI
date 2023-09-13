"""  This is the main script where all the program functions are executed """

# import 

import argparse
import sys
import protein_pdb
import hydrophobicity
import sphere
import start_plane
import search_best_plane
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
import numpy as np




if __name__ == "__main__":

    # use argparse to have the user provide arguments on the command line
    parser = argparse.ArgumentParser(description="Determine a membrane design that maximizes relative hydrophobicity. .")

    parser.add_argument('file_name', type=str, help='specify pdb file name (../data/file.pdb).')
    parser.add_argument('-n', '--number-of-points', type=int, required=True,
                        help='number of points on the half-sphere.')
    parser.add_argument('-o', '--output-file', type=str, required=True,
                        help='output pdb file name.')

    args = parser.parse_args()
    file_name = args.file_name
    num_points = args.number_of_points
    output_file = args.output_file
    
    # read pdb file and obtain center of mass, solvent accessibility, and CA list
    print("RUN")

    pdb_output = protein_pdb.PDB_manipulation(file_name)
    center_of_mass = pdb_output[0]
    all_CA = pdb_output[1]

    # using the hydrophobicity module

    # list of accessible aa
    accessible_list = hydrophobicity.sep_aa_accessibility(all_CA)

    # total number of hydrophobic protein residues
    all_hydrophobe = hydrophobicity.hydrophobe_count(all_CA)

    total = hydrophobicity.hydro_count(accessible_list)

    # using the sphere module for half-sphere construction

    radius = 10.0  # half-sphere radius
    sphere_points = sphere.saff_kuijlaars_half_sphere(num_points*2, radius, center_of_mass)

    # using star_plane module

    # calculate normal vector
    normal_vectors = start_plane.calc_vecteur_normal(sphere_points, center_of_mass)

    # calculating plane equations from normal vectors
    cartesian_eq = start_plane.d_planes_coords(normal_vectors, center_of_mass)

    # search for the best membrane position

    # function to find for each axis the membrane maximizing relative hydrophobicity
    each_best = search_best_plane.max_hydroph(normal_vectors, cartesian_eq, accessible_list, total, all_hydrophobe)

    # searches for the membrane(s) with the best hydrophobicity value
    search_best_val = search_best_plane.best_plane(each_best, normal_vectors)

    # enlarged membrane thickness for the planes for above
    optimize_best = search_best_plane.optimize_hydrophobicity(search_best_val, total, accessible_list, all_hydrophobe)
    print(f"Best membrane plane, relative hydrophobicity value {optimize_best[0][5]}")
    
    # This function produces a pdb file of the best membrane plane
    protein_pdb.build_memb(optimize_best, file_name, output_file)
    print(f"the pdb file {output_file} has been created")
