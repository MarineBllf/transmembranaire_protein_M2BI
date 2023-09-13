""" These functions can be used to search for the best membrane planes based on relative hydrophobicity."""


def membrane_position(normal_vector,accessible_aa,d_up,d_down):

    """ This function lets you know whether the planes are still positioned on the protein.

     Parameters
    ----------
    normal vector coordinates
    list of solvent-accessible aa
    d values: for each plane delimiting the membrane
    
    Returns
    ----------
    if the membrane is still positioned on the protein 
    the list of accessible aa included in the plan in question

    """
    # initialization of the list of aa accessible in between the two plans
    id_aa = []
    # initialization of the position as FALSE, if at least one accessible aa is in the plane, then the position becomes TRUE
    memb_position = "FALSE"

    for i in range (0,len(accessible_aa)):
     # the sign of the values obtained shows whether the aa is still between the 2 planes
        high_val = normal_vector[0]*accessible_aa[i][2]+normal_vector[1]*accessible_aa[i][3]+normal_vector[2]*accessible_aa[i][3]+d_up
        down_val = normal_vector[0]*accessible_aa[i][2]+normal_vector[1]*accessible_aa[i][3]+normal_vector[2]*accessible_aa[i][3]+d_down
        if high_val <= 0 and down_val >= 0 : 
            memb_position = "TRUE"
            id_aa.append(i)
        elif high_val >= 0 and down_val <= 0 : 
            memb_position = "TRUE"
            id_aa.append(i)

    return (memb_position,id_aa)




def calcul_hydrophobicite_rel(list_accessible_aa,id_aa,totaux,all_hydrophobe):
    """ This function calculates the relative hydrophobicity between 2 planes.

         Parameters
        ----------
        list of accessible aa 
        list of accessible amino acid identifiers between the two planes
        total of accessible hydrophobes and hydrophiles
    
         Returns
        ----------
        relative hydrophobicity value

    """
    
    # initialization of hydrophobic and hydrophilic membrane counts
    hydrophobe_memb = 0
    hydrophile_memb = 0

    for i in id_aa : 
        if list_accessible_aa[i][7] == "hphobic" : 
            hydrophobe_memb += 1
        elif list_accessible_aa[i][7] == "hphilie" :
            hydrophile_memb += 1
    # calculation of the number of hydrophiles exposed outside the membrane
    hphile_out_memb = (totaux[0] - hydrophile_memb)
    relative_hydrophobicity =  ((hphile_out_memb /totaux[0]) + (hydrophobe_memb/all_hydrophobe))
    
    return relative_hydrophobicity


def max_hydroph(list_norm_vect,eq_cartesian1,list_accessible,tot,all_hydrophobe):
    """ For each axis, this function determines the position of the 15-angstrom membrane
        maximizing relative hydrophobicity.
        
        On each axis, the membrane is slid from top to bottom on the axis 
        to achieve maximum relative hydrophobicity.
        
    Parameters
    ----------
    coordinates of normal vectors and values of d 
    list of solvent-accessible aa
    total of accessible hydrophobes and hydrophiles
    
    Returns
    ----------
    a list of lists 
    each_best_plane[0]: for each axis, returns the starting d values of the plane equation 
    each_best_plane[1][0:2]: d values for each axis where relative hydrophobicity is high 
    each_best_plane[1][2]: the maximum hydrophobicity values found between planes for each axis
        
        """
    # this list will be updated for each axis if a higher relative hydrophobicity value is found
    each_best_plane = eq_cartesian1.copy()
    
    for i in range (0 , len(list_norm_vect)):
        # at the start, every shot is on the protein, so TRUE
        top_slice = "TRUE"
        down_slice  = "TRUE"
        # for a given axis, the top (m1) and bottom (m2) of the axis will be studied starting from the center of mass 
        d_top_m1 = each_best_plane[i][0][0]
        d_top_m2 = each_best_plane[i][0][0]
        d_dwn_m1 = each_best_plane[i][0][1]
        d_dwn_m2 = each_best_plane[i][0][1]
        while top_slice == "TRUE" : 
            print(top_slice)
            # verification of membrane position (on or off the protein) 
            top_plane = membrane_position(list_norm_vect[i],list_accessible,d_top_m1,d_dwn_m1)
            top_slice = top_plane[0]  # TRUE or FALSE
            id_aa_memb = top_plane[1]  # id of accessible residuals between the two planes 
            if top_slice == "TRUE" : 
                # calculation of relative hydrophobicity
                rel_hydrophobicity  = calcul_hydrophobicite_rel(list_accessible,id_aa_memb,tot,all_hydrophobe)
                # comparison of the relative hydrophobicity value with the best found so far for this axis
                if rel_hydrophobicity > each_best_plane[i][1][2]:
                    # if the value is better, it is saved and the position of the corresponding plane is also saved
                    each_best_plane[i][1][2] = rel_hydrophobicity
                    each_best_plane[i][1][0] = d_top_m1
                    each_best_plane[i][1][1] = d_dwn_m1
                    # the membrane is shifted by 1 on axis
                    d_top_m1 += 1 
                    d_dwn_m1 += 1
                else : 
                    # if the hydrophobicity value found is no higher than the current best
                    d_top_m1 += 1 
                    d_dwn_m1 += 1          
        # the bottom of the same axis is explored after the top
        while down_slice == "TRUE":
            dwn_plane = membrane_position(list_norm_vect[i],list_accessible,d_top_m2,d_dwn_m2)
            down_slice = dwn_plane[0]  # TRUE or FALSE
            id_aa_memb = dwn_plane[1]
            if down_slice == "TRUE":  # accessible aa are between the two planes
                rel_hydrophobicity_dwn  = calcul_hydrophobicite_rel(list_accessible,id_aa_memb,tot,all_hydrophobe)
                if rel_hydrophobicity_dwn > each_best_plane[i][1][2]: 
                    each_best_plane[i][1][2] = rel_hydrophobicity_dwn
                    each_best_plane[i][1][0] = d_top_m2
                    each_best_plane[i][1][1] = d_dwn_m2 
                    d_top_m2 -= 1 
                    d_dwn_m2 -= 1
                else :
                    d_top_m2 -= 1 
                    d_dwn_m2 -= 1
    # for each planes returns the position which maximise relative hydrophobicity values

    return (each_best_plane)


def max_hydroph2(list_norm_vect,eq_cartesian1,list_accessible,tot,all_hydrophobe,list_CA):
    """ For each axis, this function determines the position of the 15-angstrom membrane
        maximizing relative hydrophobicity.
        
        On each axis, the membrane is slid from top to bottom on the axis 
        to achieve maximum relative hydrophobicity.
        
    Parameters
    ----------
    coordinates of normal vectors and values of d 
    list of solvent-accessible aa
    total of accessible hydrophobes and hydrophiles
    
    Returns
    ----------
    a list of lists 
    each_best_plane[0]: for each axis, returns the starting d values of the plane equation 
    each_best_plane[1][0:2]: d values for each axis where relative hydrophobicity is high 
    each_best_plane[1][2]: the maximum hydrophobicity values found between planes for each axis
        
        """
    # this list will be updated for each axis if a higher relative hydrophobicity value is found
    each_best_plane = eq_cartesian1.copy()
    
    for i in range (0 , len(list_norm_vect)):
        print(i)
        # at the start, every shot is on the protein, so TRUE
        slice_protein_top = "TRUE" 
        slice_protein_down = "TRUE" 
        top_slice_acc = "TRUE"
        down_slice_acc  = "TRUE"
        # for a given axis, the top (m1) and bottom (m2) of the axis will be studied starting from the center of mass 
        d_top_m1 = each_best_plane[i][0][0]
        d_top_m2 = each_best_plane[i][0][0]
        d_dwn_m1 = each_best_plane[i][0][1]
        d_dwn_m2 = each_best_plane[i][0][1]
        while slice_protein_top == "TRUE" : 
            # verification of membrane position (on or off the protein) 
            top_plane = membrane_position(list_norm_vect[i],list_accessible,d_top_m1,d_dwn_m1,list_CA)
            slice_protein_top = top_plane[0]
            print("top")
            top_slice_acc = top_plane[1]  # TRUE or FALSE
            id_aa_memb = top_plane[2]  # id of accessible residuals between the two planes 
            if slice_protein_top == "TRUE":
                if top_slice_acc == "TRUE" : 
                    # calculation of relative hydrophobicity
                    rel_hydrophobicity  = calcul_hydrophobicite_rel(list_accessible,id_aa_memb,tot,all_hydrophobe)
                    # comparison of the relative hydrophobicity value with the best found so far for this axis
                    if rel_hydrophobicity > each_best_plane[i][1][2]:
                    # if the value is better, it is saved and the position of the corresponding plane is also saved
                        each_best_plane[i][1][2] = rel_hydrophobicity
                        each_best_plane[i][1][0] = d_top_m1
                        each_best_plane[i][1][1] = d_dwn_m1
                        # the membrane is shifted by 1 on axis
                        d_top_m1 += 1 
                        d_dwn_m1 += 1
                else : 
                    # if the hydrophobicity value found is no higher than the current best
                    top_slice_acc == "TRUE"
                    d_top_m1 += 1 
                    d_dwn_m1 += 1          
        # the bottom of the same axis is explored after the top

        while slice_protein_down == "TRUE":
            dwn_plane = membrane_position(list_norm_vect[i],list_accessible,d_top_m2,d_dwn_m2,list_CA)
            slice_protein_down = dwn_plane[0]
            down_slice_acc = dwn_plane[1]  # TRUE or FALSE
            id_aa_memb = dwn_plane[2]
            if slice_protein_down == "TRUE" : 
                if down_slice_acc == "TRUE":  # accessible aa are between the two planes
                    rel_hydrophobicity_dwn  = calcul_hydrophobicite_rel(list_accessible,id_aa_memb,tot,all_hydrophobe)
                    if rel_hydrophobicity_dwn > each_best_plane[i][1][2]: 
                        print(rel_hydrophobicity_dwn)
                        each_best_plane[i][1][2] = rel_hydrophobicity_dwn
                        each_best_plane[i][1][0] = d_top_m2
                        each_best_plane[i][1][1] = d_dwn_m2 
                        d_top_m2 -= 1 
                        d_dwn_m2 -= 1
                else :
                    down_slice_acc = "TRUE"
                    d_top_m2 -= 1 
                    d_dwn_m2 -= 1
    # for each planes returns the position which maximise relative hydrophobicity values

    return (each_best_plane)



def best_plane(list_best_4_each,list_nvector):
    """ This function lets you know which plane has the highest relative hydrophobicity.
        If several planes have the same value, they are taken into account.
        
    Parameters
    ----------
    the list containing for each plane the max relative hydrophobicity value found 
    and the corresponding membrane position
    list of normal vectors
    
    Returns
    ----------
    List of coordinates of the plane(s) with the best relative hydrophobicity 
    and the best value.
        
    """
        
    # initialization of the list of best hydrophobicity values for each plane   
    best_hphobe_values = []

    # variable initialization for the best value among all axes
    best_value = 0 
    # id of planes with the best value
    id_best = []
    for i in range(0,len(list_best_4_each)):
        best_hphobe_values.append(list_best_4_each[i][1][2])
   
    for i in range(0,len(best_hphobe_values)):
        if best_hphobe_values[i] >= best_value :
            best_value = best_hphobe_values[i]
    print(best_value)
    for i in range(0,len(best_hphobe_values)):
        if best_hphobe_values[i] >= best_value :
            id_best.append(i)
    print(id_best)
    # equation for the best plan or plans
    best_planes = []
    for i in range (0,len(id_best)):
        memb_id = id_best[i]
        # coords ax,by,cz,d for the best membrane or membranes
        ax= list_nvector[memb_id][0]
        by= list_nvector[memb_id][1]
        cz= list_nvector[memb_id][2]
        dh= list_best_4_each[memb_id][1][0]
        db= list_best_4_each[memb_id][1][1]
        # best relative hydrophobicity value
        rel_hydrophobicity = list_best_4_each[memb_id][1][2]
        best_planes.append([ax, by, cz, dh, db, rel_hydrophobicity])
        print(best_planes)
    return best_planes




def optimize_hydrophobicity2(list_best_planes1,tot,accessible_aa,all_hydrophobe): 
    """ for the plane(s) with the best relative hydrophobicity value, this function increases
      the membrane thickness to optimize the hydrophobicity value

    If the relative hydrophobicity value decreases too much or converges, the increase is stopped.

    Parameters
    ----------
    coordinates of best plane(s) and best relative hydrophobicity value 
    list of total hydrophilic and hydrophobic numbers
    list of solvent-accessible aa
    
    Returns
    ----------
    best membrane coordinates (two planes)
        
    """

    # this list will be updated for each axis if a higher relative hydrophobicity value is found
    list_best_planes = list_best_planes1.copy()
    # if only one membrane plane has the best hydrophobicity value
    if len(list_best_planes) == 1 : 
        dtop = list_best_planes[0][3] 
        ddown  = list_best_planes[0][4]
        # planes are displaced in 0.5 increments
        dtop += 0.5
        ddown -= 0.5
        norm_vector = list_best_planes[0][0:3] 
        # the limit value is the number of iterations made before stopping to increase the membrane thickness 
        # the number decreases if the relative hydrophobicity decreases too much or converges
        limit  = 10
        
        while limit != 0 : 
            memb_slice = membrane_position(norm_vector,accessible_aa,dtop,ddown)
            id_aa_memb = memb_slice[1]
            rel_hydrophobicity  = calcul_hydrophobicite_rel(accessible_aa,id_aa_memb,tot,all_hydrophobe)
            # if the relative hydrophobicity is better than the current one, the list is updated
            if rel_hydrophobicity > list_best_planes[0][5]:
                limit = 10 
                list_best_planes[0][5] = rel_hydrophobicity
                list_best_planes[0][3] = dtop
                list_best_planes[0][4] = ddown
                dtop += 0.5
                ddown -= 0.5
            else :
                limit -= 1
                dtop  += 0.5
                ddown -= 0.5

    # if several planes have the same (best) relative hydrophocity value
    else : 
        for i in range(0,len(list_best_planes)):
            dtop  = list_best_planes[i][3] 
            ddown   = list_best_planes[i][4]
            dtop += 0.5
            ddown  -= 0.5
            norm_vector = list_best_planes[i][0:3] 
            limit = 10
           
            while limit != 0 : 
                memb_slice = membrane_position(norm_vector,accessible_aa,dtop,ddown)
                id_aa_memb = memb_slice[1]
                print(id_aa_memb)
                rel_hydrophobicity  = calcul_hydrophobicite_rel(accessible_aa,id_aa_memb,tot,all_hydrophobe)
                # if a better relative hydrophobicity value is found, it is retained
                if rel_hydrophobicity > list_best_planes[i][5]:
                    limit = 10 
                    list_best_planes[i][5] = rel_hydrophobicity
                    list_best_planes[i][3] = dtop
                    list_best_planes[i][4] = ddown
                    dtop += 0.5
                    ddown -= 0.5
                else :
                    limit -= 1
                    dtop += 0.5
                    ddown -= 0.5
        print("OH")
        # select the best if initially several had the same hydrophobicity value
        best_phobe_values = []
        for i in range(0,len(list_best_planes)):
            best_phobe_values.append((list_best_planes[i][5]))
        best_num = 0
        list_best = [] 
        for i in range(0,len(best_phobe_values)):
            if best_phobe_values[i] >= best_num :
                best_num = best_phobe_values[i]

        for i in range(0,len(best_phobe_values)):
            if best_phobe_values[i] >= best_num :
                list_best.append(i)
        # recuperer l'equation des meilleurs plan et la valeur d'hydrophobicite associe 
        best_planes = []
        for i  in range (0, len(list_best)) :
            id_best = list_best[i]
            best_planes.append(list_best_planes[id_best])
        list_best_planes = best_planes

     # if now only one plane has the best value only this membrane is kept otherwise 
     # all planes with the same value are kept
    return list_best_planes


def optimize_hydrophobicity(list_best_planes1,tot,accessible_aa,all_hydrophobe): 
    """ for the plane(s) with the best relative hydrophobicity value, this function increases
      the membrane thickness to optimize the hydrophobicity value

    If the relative hydrophobicity value decreases too much or converges, the increase is stopped.

    Parameters
    ----------
    coordinates of best plane(s) and best relative hydrophobicity value 
    list of total hydrophilic and hydrophobic numbers
    list of solvent-accessible aa
    
    Returns
    ----------
    best membrane coordinates (two planes)
        
    """

    # this list will be updated for each axis if a higher relative hydrophobicity value is found
    list_best_planes = list_best_planes1.copy()
    # if only one membrane plane has the best hydrophobicity value
    if len(list_best_planes) == 1 : 
        dtop = list_best_planes[0][3] 
        ddown  = list_best_planes[0][4]
        # planes are displaced in 0.5 increments
        dtop += 0.5
        ddown -= 0.5
        norm_vector = list_best_planes[0][0:3] 
        # the limit value is the number of iterations made before stopping to increase the membrane thickness 
        # the number decreases if the relative hydrophobicity decreases too much or converges
        limit  = 20
        
        while limit != 0 : 
            memb_slice = membrane_position(norm_vector,accessible_aa,dtop,ddown)
            id_aa_memb = memb_slice[1]
            rel_hydrophobicity  = calcul_hydrophobicite_rel(accessible_aa,id_aa_memb,tot,all_hydrophobe)
            # if the relative hydrophobicity is better than the current one, the list is updated
            if rel_hydrophobicity > list_best_planes[0][5]:
                limit = 20 
                list_best_planes[0][5] = rel_hydrophobicity
                list_best_planes[0][3] = dtop
                list_best_planes[0][4] = ddown
                dtop += 0.5
                ddown -= 0.5
            else :
                limit -= 1
                dtop  += 0.5
                ddown -= 0.5

    # if several planes have the same (best) relative hydrophocity value
    else : 
        for i in range(0,len(list_best_planes)):
            dtop  = list_best_planes[i][3] 
            ddown   = list_best_planes[i][4]
            dtop += 0.5
            ddown  -= 0.5
            norm_vector = list_best_planes[i][0:3] 
            limit = 20
           
            while limit != 0 : 
                memb_slice = membrane_position(norm_vector,accessible_aa,dtop,ddown)
                id_aa_memb = memb_slice[1]
                print(id_aa_memb)
                rel_hydrophobicity  = calcul_hydrophobicite_rel(accessible_aa,id_aa_memb,tot,all_hydrophobe)
                # if a better relative hydrophobicity value is found, it is retained
                if rel_hydrophobicity > list_best_planes[i][5]:
                    limit = 20 
                    list_best_planes[i][5] = rel_hydrophobicity
                    list_best_planes[i][3] = dtop
                    list_best_planes[i][4] = ddown
                    dtop += 0.5
                    ddown -= 0.5
                else :
                    limit -= 1
                    dtop += 0.5
                    ddown -= 0.5
        print("OH")
        # select the best if initially several had the same hydrophobicity value
        best_phobe_values = []
        for i in range(0,len(list_best_planes)):
            best_phobe_values.append((list_best_planes[i][5]))
        best_num = 0
        list_best = [] 
        for i in range(0,len(best_phobe_values)):
            if best_phobe_values[i] >= best_num :
                best_num = best_phobe_values[i]

        for i in range(0,len(best_phobe_values)):
            if best_phobe_values[i] >= best_num :
                list_best.append(i)
        # recuperer l'equation des meilleurs plan et la valeur d'hydrophobicite associe 
        best_planes = []
        for i  in range (0, len(list_best)) :
            id_best = list_best[i]
            best_planes.append(list_best_planes[id_best])
        list_best_planes = best_planes

     # if now only one plane has the best value only this membrane is kept otherwise 
     # all planes with the same value are kept
    return list_best_planes