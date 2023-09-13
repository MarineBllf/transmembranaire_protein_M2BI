""" Module for calculating hydrophobicity """


# list of aa considered hydrophilic or hydrophobic
hydrophobic = ["F" , "G" , "I" , "L" , "M" , "V" , "W" , "Y" ]
hydrophile = [ "A" , "C" , "D" , "E" , "H" , "K" , "N" , "P" , "Q" , "R" , "S" , "T" ]


def sep_aa_accessibility(list_CA):
    """ This function recovers the CA of residues accessible to the solvent
    and defines their hydrophilic and hydrophobic character. 
    
        Parameters
    ----------
    list_CA : the list of protein CAs

    Returns
    ----------
    The program returns a list of lists.
    Each list contains for a residue i 
        accessible_aa[i][0]   : resname
        accessible_aa[i][1]   : residue id 
        accessible_aa[i][2:4] : residue coordinates
        accessible_aa[i][5]   : accessible (A) or buried (E)
        accessible_aa[i][5]   : one-letter aa code 
        accessible_aa[i][6]   : hydrophile (hphile) or hydrophobic (hphobic)
    """
   
    # initialization of accessible aa list
    accessible_aa = []
    
    for i in range (0,len(list_CA)):
        print(list_CA)
    #  caracterization of accessible hydrophilic and hydrophobic aa
        if list_CA[i][5] == "A" and (list_CA[i][6] in hydrophile) :
            list_CA[i].append("hphile")
            accessible_aa.append(list_CA[i])
        
        elif list_CA[i][5] == "A" and list_CA[i][6] in hydrophobic :
            list_CA[i].append("hphobic")
            accessible_aa.append(list_CA[i])
    
    return(accessible_aa)



def hydrophobe_count(liste_CA): 
    """ A function to determine the total number of hydrophobic aa in the protein.

    Parameters
    ----------
    list_CA : the list of protein CAs with their caracteristic.

    Returns
    ----------
    The total number of hydrophobic residues in the protein. 
    
    """

    # initialization of hydrophobic count
    hydrophobe_count = 0 

    for i in range (0,len(liste_CA)):
        if liste_CA[i][6] in hydrophobic : 
            hydrophobe_count += 1 

    return(hydrophobe_count)

def hydro_count(list_accessible):
    """ A function to determine the total number of hydrophobic and hydrophile accessible aad.

    Parameters
    ----------
    list of accessible aa 

    Returns
    ----------
    The total number of hydrophobic and hydrophile accessible residues in the protein. 
    
    """
     # initiatization of hydrophile and hydrophobe count 
    tot_hydrophile_exp = 0
    tot_hydrophob_exp = 0 

    for i in range (0,len (list_accessible)):
        if list_accessible[i][7] == "hphobic" : 
            tot_hydrophob_exp += 1
        else : 
            tot_hydrophile_exp += 1 
    totaux = []
    totaux.append(tot_hydrophile_exp)
    totaux.append(tot_hydrophob_exp)
    
    return totaux