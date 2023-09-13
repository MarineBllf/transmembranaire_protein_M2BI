""" Determining eigenvectors and starting planes . """



def calc_vecteur_normal(array_coord_pts,mass_center):
    """ The function is used to determine the coordinates of the normal vectors .

    Parameters
    ----------
    the coordinates of the points on the half-sphere and the center of mass.
    The vectors are determined between each point of the half-sphere and 
    the center of mass.

    Returns
    ----------
    The total number of hydrophobic residues in the protein. 
    
    """

    # initialization of the list of normal vectors 
    vnorm_coords=[]

    for i in range(0,len(array_coord_pts)):
    # calculation of the x ,y ,z coordinates of the vector 
        x_vnor = mass_center[0] - array_coord_pts[i][0]
        y_vnor = mass_center[1] - array_coord_pts[i][1]
        z_vnor = mass_center[2] - array_coord_pts[i][2]
        niem_vnorm = [x_vnor,y_vnor, z_vnor]
        #  list of coordinates of each vector added to the list of all vectors
        vnorm_coords.append(niem_vnorm)           
   
    return vnorm_coords



def d_planes_coords(norm_vec,mass_center):
    """ the function is used to determine the d value of plane equations .

    Parameters
    ----------
    The normal vector whose coordinates correspond to a, b, c in the equation 
    of the plane to which it is normal (equation: ax+by+cz+d=0).
    The center of mass 

    Returns
    ----------
    A list of list 
    cartesian_equation[0] : d values of the plane delimiting each hypothetical
                            membrane for each normal vector
    cartesian_equation[1] : empty list that will contain for each plane : 
        cartesian_equation[1][0:2] : the d values of the membrane with the best hydrophobicity value 
        cartesian_equation[1][2] :the hydrophobicity value corresponding to that plane """
    cartesian_equation = []
    
    for i in range (0,len(norm_vec)): 
        # calculation of d in the equation ax+by+cz+d, taking the center of mass as the point
        d=norm_vec[i][0]*mass_center[0]+norm_vec[i][1]*mass_center[1]+norm_vec[i][2]*mass_center[2]
        d=round((-1*d),3)
        # the two planes delimiting the membrane are 15 angstroms apart
        cartesian_equation.append([[d+7.5,d-7.5],[0,0,0]])

    return cartesian_equation


