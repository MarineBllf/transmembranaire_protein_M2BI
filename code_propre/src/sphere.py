""" Calculating points on the surface of a half-sphere."""

import numpy as np

def saff_kuijlaars_half_sphere(num_points, radius,mass_center):
    """This function calculates the points on the surface of a half-sphere using the Saaf and Kuijlarrs algorithm
    
    Parameters
    ----------
    num_points :the number of points on the surface of the half-sphere
    radius : sphere radius 
    mass_center : protein center of mass
    
    Returns
    ----------
    x ,y, z coordinates of points on the surface of the half sphere 
    center on the center of mass 
    
    """

    # point list initialization
    points = []

    phi_0 = 0
    for k in range(2, num_points//2):  # total number of points // 2 to obtain a half sphere
        # hk, theta and phi are the parameters of the Saaf and Kuijlaars algorithm
        hk = -1 + ((2 * (k - 1)) / (num_points - 1))
        theta = np.arccos(hk)
        phi = (phi_0 + ((3.6) / (num_points) ** (1 / 2)) * (1 / (1 - hk ** 2) ** (1 / 2)))% (2*3.14)
        phi_0 = phi
        x = radius * np.sin(theta) * np.cos(phi)
        y = radius * np.sin(theta) * np.sin(phi)
        z = radius * np.cos(theta)
        points.append([x, y, z])

    # the half-sphere is centered on the center of mass
    half_sphere_points_centered = [np.array(point) + np.array(mass_center) for point in points]
    
    return half_sphere_points_centered

