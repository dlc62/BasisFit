#!/usr/bin/python

import sys
import numpy
import math

def overlap(e1,e2,l):
    s = math.pow(2.0*e1,0.75) * math.pow(2.0*e2,0.75) * math.pow(e1+e2,-(1.5+l))
    return s

def large_large_matrix(shells,angmom):
# note: should be abstracted to generalized function
    S_matrix = []
    for i in range(0,len(shells)):
        S_vec = []
        s1 = shells[i]
	for j in range(0,len(shells)):
	    s2 = shells[j]
	    S_element = 0
	    for k in range(0,len(s1)):
                e1 = s1[k][0]
		c1 = s1[k][1]
		for l in range(0,len(s2)):
		    e2 = s2[l][0]
		    c2 = s2[l][1]
		    S_int = overlap(e1,e2,angmom)
		    S_element += c1*c2*S_int
	    S_vec.append(S_element)
	S_matrix.append(S_vec)
    return S_matrix

def small_large_vector(s1,shells,angmom):
    T_vec = []
    for j in range(0,len(shells)):
	s2 = shells[j]
	T_element = 0
	for k in range(0,len(s1)):
            e1 = s1[k][0]
            c1 = s1[k][1]
	    for l in range(0,len(s2)):
		e2 = s2[l][0]
		c2 = s2[l][1]
		T_int = overlap(e1,e2,angmom)
		T_element += c1*c2*T_int
	T_vec.append([T_element])
    return T_vec

def small_large_fit(S_matrix,T_vec):
# convert S and T to numpy matrix format
    A_matrix = numpy.matrix(S_matrix)
    B_vec = numpy.matrix(T_vec)
# find coefficient vector, append to s_coeffs 
    C_vec = numpy.linalg.solve(A_matrix,B_vec)
    U_vec = numpy.array(C_vec).flatten().tolist()
    return U_vec

def do(min_bas_dat,lrg_bas_dat):

# unpack basis set data shell by shell
    [s_min,p_min,d_min,f_min,g_min,h_min] = min_bas_dat
    [s_lrg,p_lrg,d_lrg,f_lrg,g_lrg,h_lrg] = lrg_bas_dat

# do large-to-small fitting on shell-by-shell basis
#------------------------------------------------------
# s functions
#------------------------------------------------------
    if ((s_min != []) and (s_lrg != [])):
        angmom = 0.0
        shells = s_lrg
        s_fit_coeffs = []
# construct overlap matrix for pairs of contracted functions from large basis
        S_matrix = large_large_matrix(shells,angmom)
# construct overlap vector between each minimal contracted basis function 
# and contracted large basis functions
        for i in range(0,len(s_min)):
            s1 = s_min[i]
	    T_vec = small_large_vector(s1,shells,angmom)
	    U_vec = small_large_fit(S_matrix,T_vec) 
	    s_fit_coeffs.append(U_vec)
    else:
	s_fit_coeffs = []

#------------------------------------------------------
# p functions
#------------------------------------------------------
    if ((p_min != []) and (p_lrg != [])):
        angmom = 1.0
        shells = p_lrg
        p_fit_coeffs = []
# construct overlap matrix for pairs of contracted functions from large basis
        S_matrix = large_large_matrix(shells,angmom)
# construct overlap vector between each minimal contracted basis function 
# and contracted large basis functions
        for i in range(0,len(p_min)):
            s1 = p_min[i]
	    T_vec = small_large_vector(s1,shells,angmom)
	    U_vec = small_large_fit(S_matrix,T_vec) 
	    p_fit_coeffs.append(U_vec)
    else:
	p_fit_coeffs = []

#------------------------------------------------------
# d functions
#------------------------------------------------------
    if ((d_min != []) and (d_lrg != [])):
        angmom = 2.0
        shells = d_lrg
        d_fit_coeffs = []
# construct overlap matrix for pairs of contracted functions from large basis
        S_matrix = large_large_matrix(shells,angmom)
# construct overlap vector between each minimal contracted basis function 
# and contracted large basis functions
        for i in range(0,len(d_min)):
            s1 = d_min[i]
	    T_vec = small_large_vector(s1,shells,angmom)
	    U_vec = small_large_fit(S_matrix,T_vec) 
	    d_fit_coeffs.append(U_vec)
    else:
	d_fit_coeffs = []

#------------------------------------------------------
# f functions
#------------------------------------------------------
    if ((f_min != []) and (f_lrg != [])):
        angmom = 3.0
        shells = f_lrg
        f_fit_coeffs = []
# construct overlap matrix for pairs of contracted functions from large basis
        S_matrix = large_large_matrix(shells,angmom)
# construct overlap vector between each minimal contracted basis function 
# and contracted large basis functions
        for i in range(0,len(f_min)):
            s1 = f_min[i]
	    T_vec = small_large_vector(s1,shells,angmom)
	    U_vec = small_large_fit(S_matrix,T_vec) 
	    f_fit_coeffs.append(U_vec)
    else:
	f_fit_coeffs = []

#------------------------------------------------------
# g functions
#------------------------------------------------------
    if ((g_min != []) and (g_lrg != [])):
        angmom = 4.0
        shells = g_lrg
        g_fit_coeffs = []
# construct overlap matrix for pairs of contracted functions from large basis
        S_matrix = large_large_matrix(shells,angmom)
# construct overlap vector between each minimal contracted basis function 
# and contracted large basis functions
        for i in range(0,len(g_min)):
            s1 = g_min[i]
	    T_vec = small_large_vector(s1,shells,angmom)
	    U_vec = small_large_fit(S_matrix,T_vec) 
	    g_fit_coeffs.append(U_vec)
    else:
	g_fit_coeffs = []

#------------------------------------------------------
# g functions
#------------------------------------------------------
    if ((h_min != []) and (h_lrg != [])):
        angmom = 5.0
        shells = h_lrg
        h_fit_coeffs = []
# construct overlap matrix for pairs of contracted functions from large basis
        S_matrix = large_large_matrix(shells,angmom)
# construct overlap vector between each minimal contracted basis function 
# and contracted large basis functions
        for i in range(0,len(h_min)):
            s1 = h_min[i]
	    T_vec = small_large_vector(s1,shells,angmom)
	    U_vec = small_large_fit(S_matrix,T_vec) 
	    h_fit_coeffs.append(U_vec)
    else:
	h_fit_coeffs = []

    return [s_fit_coeffs,p_fit_coeffs,d_fit_coeffs,f_fit_coeffs,g_fit_coeffs,h_fit_coeffs]
