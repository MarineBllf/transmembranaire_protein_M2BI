from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
import numpy as np


masses = {"N": 14, "C": 12, "CA": 12,"O": 16,"S": 12, "H":1}

def PDB_manipulation(pdb_name): 
    """ From a pdb file, this function allows you to file to extract CA 
        calculate the protein's center of mass 
        and solvent accessibility for each atom 

    Parameters
    ----------
    pdb_name : name of the protein pdb file in the directory 

    Returns
    ----------
    The program returns a list of lists.
    outpout_pdb[0] = center of mass cordinates [x , y, z]
    outpout_pdb[1] = list of CA, including for each :
        outpout_pdb[1][0]   : resname
        outpout_pdb[1][1]   : residue id 
        outpout_pdb[1][2:4] : residue coordinates
        outpout_pdb[1][5]   : accessible (A) or buried(E)
    """ 

    parser = PDBParser() 
    structure = parser.get_structure("PHA-L", pdb_name)

    # variable initialization

    all_CA = []  # protein CA lists

    # mass center coordinates, numerator (num) denominator (den) in the formula
    xnum, xden = 0, 0
    ynum, yden = 0, 0 
    znum, zden = 0, 0

    # extracts information from each residue and atom 
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    xi = f"{atom.coord[0]:.3f}" 
                    yi = f"{atom.coord[1]:.3f}" 
                    zi = f"{atom.coord[2]:.3f}" 
                    # recovers the atomic mass of the corresponding atom
                    mi = masses[atom.element]
                    
                    # sum of terms in center-of-mass formula
                    xnum += (mi * float(xi))
                    xden += mi
                    ynum += (mi * float(yi))
                    yden += mi
                    znum += (mi * float(zi))
                    zden += mi

                    # if an atom is a CA, it is added
                    if atom.name == 'CA' :
                        x_coord = atom.coord[0]
                        y_coord = atom.coord[1]
                        z_coord = atom.coord[2]
                        all_CA.append([residue.resname,residue.id[1],x_coord,y_coord,z_coord])
    
    #  final center-of-mass calculation
    xg = xnum /xden
    yg = ynum /yden
    zg = znum /zden
    center_of_mass = [xg, yg, zg]

    # calculating solvent accessibility with DSSP
    dssp = DSSP(structure[0],pdb_name, dssp='mkdssp')
    
    for i in range (0,len(all_CA)):
        print(dssp.keys())
        # extraction of the solvent accessibility value of the residue 
        a_key = list(dssp.keys())[i]
        # if value > 0.3  accessible (A) otherwise buried (E)
        if dssp[a_key][3] > 0.3: 
            all_CA[i].append("A")
            all_CA[i].append(dssp[a_key][1])  # extraction of one-letter residue code
        else: 
            all_CA[i].append("E")
            all_CA[i].append(dssp[a_key][1])
            
    outpout_pdb = []
    outpout_pdb.append(center_of_mass)
    outpout_pdb.append(all_CA)
    return (outpout_pdb)


def build_memb(planes_eq,pdb_file,outpout_name): 
    """ This function constructs a pdb with a 3D representation of the membrane found.

        Parameters
        ----------
        planes_eq : equation of the best plans found
        pdb_name : name of the protein pdb file in the directory 
        outpout_name : name of the outpout file (.pdb)
    
        Returns
        ----------
        The program returns a PDB files with x,y,z coordinates of the membranes atoms
    """ 
    
    # initialization of the list of x and y coordinates of membrane atoms
    x_list = []
    y_list = []
    # the minimum and maximum x and y values of the protein atoms are determined 
    # so that the plane is desined at the protein level
    parser = PDBParser() 
    structure = parser.get_structure("PHA-L", pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x_coord = atom.coord[0]
                    y_coord = atom.coord[1]
                    x_list.append(x_coord)
                    y_list.append(y_coord)
    xmin = min(x_list)
    xmax = max(x_list)
    ymin = min(y_list)
    ymax = max(y_list)
    
    # initialize the list of atom coordinates for each plane
    plane1 = []
    plane2 = []
    
    # calculates z coordinates for x and y values
    for x in np.arange(xmin-10 , xmax+10,1 ) :
        for y in np.arange(ymin-10  , ymax+10,1 ):
            z1 = (-planes_eq[0][0] * x - planes_eq[0][1] * y - planes_eq[0][3]) / planes_eq[0][2]
            z2 = (-planes_eq[0][0] * x - planes_eq[0][1] * y - planes_eq[0][4]) / planes_eq[0][2]
            plane1.append((x, y, z1))
            plane2.append((x, y, z2))
    
    # writing atomic coordinates to a pdb file
    with open(outpout_name, 'w') as f:
            for i in range (0,len(plane1)):
                    f.write(f"ATOM  {i+1:4}  CA  MEM M   1    {plane1[i][0]:8.3f}{plane1[i][1]:8.3f}{plane1[i][2]:8.3f}\n")
            for i in range (0,len(plane1)):
                    f.write(f"ATOM  {i+1:4}  CA  MEM M   1    {plane2[i][0]:8.3f}{plane2[i][1]:8.3f}{plane2[i][2]:8.3f}\n")             