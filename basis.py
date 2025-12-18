#!/usr/bin/python

import sys
import string

def get(atom,basis):

    #------------------------------------------------------------
    # Whole pile of set up stuff
    #------------------------------------------------------------

    # open appropriate library file and read in data
#    f = open("BasisFitPath/BasisData/Gamess/"+basis+".Gamess",'r')
    f = open("/BasisData/Gamess/"+basis+".Gamess",'r')
    all_data = f.readlines()
    f.close()

    # remove first and last lines
    all_data.pop(0)
    all_data.pop()
    
    # extract chunk of file pertaining to 'atom'
    got_atom = 0
    data = []
    for i in range(0,len(all_data)):
	if string.find(all_data[i],atom) != -1:
	    got_atom = 1
	    continue
        if len(all_data[i].split()) != 0 and got_atom == 1:
	    data.append(all_data[i])
	else:
	    got_atom = 0

    # exit if atom is not found in basis file 
    if len(data) == 0:
	print 'Error: ' + atom + ' data not found in basis data file'
	sys.exit()

    # set up arrays to catch output
    s_vec = []; p_vec = []; d_vec = []; f_vec = []; g_vec = []; h_vec = []

    #------------------------------------------------------------
    # Process GAMESS-style basis set data 
    #------------------------------------------------------------
    # Gaussian / QChem and Gamess
    #------------------------------------------------------------

    i = 0
    ii = 1
    ordering = []
    while i < len(data):
    # catch the last line and exit while loop if required
    # shouldn't be required though
        if i == len(data)-1:
            break
        type = data[i].split()[0]
        num = int(data[i].split()[1])
        if type == 'S': 
	    s_function = []
            for j in range(i+1,i+num+1):
	        exponent = float(data[j].split()[0+ii])
	        coefficient = float(data[j].split()[1+ii])
	        s_function.append([exponent,coefficient])
	    s_vec.append(s_function)
	    i += num+1
	    ordering.append('S')
        elif type == 'SP' or type == 'L':
	    s_function = []
	    p_function = []
            for j in range(i+1,i+num+1):
	        exponent = float(data[j].split()[0+ii])
	        coefficient_s = float(data[j].split()[1+ii])
	        coefficient_p = float(data[j].split()[2+ii])
	        s_function.append([exponent,coefficient_s])
	        p_function.append([exponent,coefficient_p])
	    s_vec.append(s_function)
	    p_vec.append(p_function)
	    i += num+1
	    ordering.append('S')
	    ordering.append('P')
        elif type == 'P': 
	    p_function = []
            for j in range(i+1,i+num+1):
	        exponent = float(data[j].split()[0+ii])
	        coefficient = float(data[j].split()[1+ii])
	        p_function.append([exponent,coefficient])
	    p_vec.append(p_function)
	    i += num+1
	    ordering.append('P')
        elif type == 'D': 
	    d_function = []
            for j in range(i+1,i+num+1):
	        exponent = float(data[j].split()[0+ii])
	        coefficient = float(data[j].split()[1+ii])
	        d_function.append([exponent,coefficient])
	    d_vec.append(d_function)
	    i += num+1
	    ordering.append('D')
        elif type == 'F': 
	    f_function = []
            for j in range(i+1,i+num+1):
	        exponent = float(data[j].split()[0+ii])
	        coefficient = float(data[j].split()[1+ii])
	        f_function.append([exponent,coefficient])
	    f_vec.append(f_function)
	    i += num+1
	    ordering.append('F')
        elif type == 'G': 
	    g_function = []
            for j in range(i+1,i+num+1):
	        exponent = float(data[j].split()[0+ii])
	        coefficient = float(data[j].split()[1+ii])
	        g_function.append([exponent,coefficient])
	    g_vec.append(g_function)
	    i += num+1
	    ordering.append('G')
        elif type == 'H': 
	    h_function = []
            for j in range(i+1,i+num+1):
	        exponent = float(data[j].split()[0+ii])
	        coefficient = float(data[j].split()[1+ii])
	        h_function.append([exponent,coefficient])
	    h_vec.append(h_function)
	    i += num+1
	    ordering.append('H')
        else:
	    print 'Error: exceeded maximum angular momentum (h)'
	    print data[i]
	    sys.exit()

    #------------------------------------------------------------
    # Return processed results
    #------------------------------------------------------------

    return [s_vec,p_vec,d_vec,f_vec,g_vec,h_vec,ordering]
