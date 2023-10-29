BAILLIF MARINE
Project : ASSIGNATION ET DETECTION DES PARTIES TRANSMEMBRANAIRES D'UNE PROTEINE


This program is used to detect the transmembrane segments of a protein.
It is used on a protein chain only. 

## Installation 

Several tools are required for the program to run smoothly.

### Dependencies
DSSP : 
This program uses DSSP software to calculate the solvent accessibility 
of amino acids. 
DSSP will be used in combination with the BioPython package. 

2 ways to install DSSP : 
- with conda, which is done directly with the environment.yml file (https://github.com/MarineBllf/transmembranaire_protein_M2BI.git)
- in the event of a problem, DSSP can be installed by 
`sudo apt-get install dssp`

Pymol : 
Pymol is not required to run the program, but it is necessary to visualize the results 
(or any other molecule visualization software).
`sudo apt-get install pymol`


For the program to work properly, the zip directory must be downloaded, unzipped.
Source : 
Change your current directory to the download directory (project) and then run the following 
command line to install the necessary environment and activate it :

```
conda env create -f  ./environment_TM.yml
conda activate env_TM
```

### Virtual environment
Now the virtual environment can be use.


### Usage 
Go to the 'src' directory to run the commands.
```
cd src/
```
You can now use the program with the following line, which requires you to enter certain arguments : 

```
python TM_main [-h] [-n N] [-o outpout]  path/to/filename

positional arguments:
  filename         A PDB file
  -n N             Number of points to place on the sphere
  -o outpout       name of output pdb file with atom coordinates of membrane planes.
```

Test data are located in the :  '../data/'

For example : 
```
python TM_main.py -n 200 -o membrane_1G90.pdb ../data/1G90.pdb
```

## Structure

The code structured in a TM_main.py where all the functions needed to run the program are imported.
- [protein_pdb.py] contains the functions required to manipulate the PDB file.
- [sphere.py] Determination of points on the surface of the half-sphere using the Saaf and Kuijlaars algorithm
- [hydrophobicity.py] module for characterizing hydrophobic residues (separation, counting)
- [start_plane.py] characterization of normal vectors and Cartesian equations of planes
- [search_best_plane.py] module for finding the best membrane plan



