import argparse
import sys
import protein_pdb
import hydrophobicity
import sphere
import start_plane
import search_best_plane

def main():
    parser = argparse.ArgumentParser(description="Description de votre script.")

    parser.add_argument('file_name', type=str, help='Nom du fichier PDB.')
    parser.add_argument('-n', '--nombre-de-points', type=int, required=True,
                        help='Le nombre de points sur la demi-sphère.')
    parser.add_argument('-o', '--output-file', type=str, required=True,
                        help='Nom du fichier de sortie.')

    args = parser.parse_args()
    file_name = args.file_name
    num_points = args.nombre_de_points
    output_file = args.output_file

    # read pdb file and obtain center of mass, solvent accessibility, and CA list
    pdb_output = protein_pdb.PDB_manipulation(file_name)
    center_of_mass = pdb_output[0]
    all_CA = pdb_output[1]

    # using the hydrophobicity module

    # list of accessible aa
    accessible_list = hydrophobicity.sep_aa_accessibility(all_CA)
    print("accessible_list")

    # total number of hydrophobic protein residues
    all_hydrophobe = hydrophobicity.hydrophobe_count(all_CA)
    print("all_hydrophobe")

    #
    total = hydrophobicity.hydro_count(accessible_list)

    # using the sphere module for half-sphere construction

    radius = 10.0  # half-sphere radius
    sphere_points = sphere.saff_kuijlaars_half_sphere(num_points, radius, center_of_mass)

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
    print(search_best_val)

    # enlarged membrane thickness for the planes for above
    optimize_best = search_best_plane.optimize_hydrophobicity(search_best_val, total, accessible_list, all_hydrophobe)

    # This function produces a pdb file of the best membrane plane
    protein_pdb.build_memb(optimize_best, file_name, output_file)
    print("hein")

if __name__ == "__main__":
    main()
