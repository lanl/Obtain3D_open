# This file contains some of the functions from the Obtain3D_open script that weren't being used (and one that was.)
import math

#%%#################################### Functions ################################
def rot(l,m,n,ang,x,y,z):
    l_unit = l/(l**2+m**2+n**2)**0.5
    n_unit = n/(l**2+m**2+n**2)**0.5
    m_unit = m/(l**2+m**2+n**2)**0.5
#
    c = math.cos(ang)
    s = math.sin(ang)
#
    a11 = l_unit*l_unit*(1-c) + c
    a12 = m_unit*l_unit*(1-c) - n_unit*s
    a13 = n_unit*l_unit*(1-c)+m_unit*s
    a21 = l_unit*m_unit*(1-c)+n_unit*s
    a22 = m_unit*m_unit*(1-c)+c
    a23 = n_unit*m_unit*(1-c)-l_unit*s
    a31 = l_unit*n_unit*(1-c)-m_unit*s
    a32 = m_unit*n_unit*(1-c)+l_unit*s
    a33 = n_unit*n_unit*(1-c)+c
    cord = [a11*x+a12*y+a13*z, a21*x+a22*y+a23*z, a31*x+a32*y+a33*z]
    return cord

# Not used in program
def cross_product(x1,y1,z1,x2,y2,z2):
    cp = [y1*z2-y2*z1, z1*x2-z2*x1, x1*y2-x2*y1]
    return cp

# Not used in program
def unit_vector(h,k,l):
    unit = [h/(h**2+k**2+l**2)**0.5, k/(h**2+k**2+l**2)**0.5, l/(h**2+k**2+l**2)**0.5]
    return unit

# Not used in program
def angle_two_vectors(h1,k1,l1,h2,k2,l2):
    a_dot_b = (h1*h2+k1*k2+l1*l2)
    mag_a = (h1**2+k1**2+l1**2)**0.5
    mag_b = (h2**2+k2**2+l2**2)**0.5
    cos_theta = a_dot_b/(mag_a*mag_b+1e-33) #replaced 0.000000000000000000000000000000001
    theta = math.acos(cos_theta)
    return theta
