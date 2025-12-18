#!/usr/bin/python

import sys
import string
import basis
import fit
import math

periodic_table = {1:'HYDROGEN', 2:'HELIUM', 3:'LITHIUM', 4:'BERYLLIUM',
		  5:'BORON', 6:'CARBON', 7:'NITROGEN', 8:'OXYGEN',
		  9:'FLUORINE', 10:'NEON', 11:'SODIUM', 12:'MAGNESIUM',
		  13:'ALUMINUM', 14:'SILICON', 15:'PHOSPHOROUS',
		  16:'SULFUR', 17:'CHLORINE', 18:'ARGON', 19:'POTASSIUM',
		  20:'CALCIUM', 21:'SCANDIUM', 22:'TITANIUM', 23:'VANADIUM',
		  24:'CHROMIUM', 25:'MANGANESE',26:'IRON',27:'COBALT',
		  28:'NICKEL',29:'COPPER',30:'ZINC',31:'GALLIUM',
		  32:'GERMANIUM',33:'ARSENIC',34:'SELENIUM',35:'BROMINE',
		  36:'KRYPTON',37:'RUBIDIUM',38:'STRONTIUM',39:'YTTRIUM',
		  40:'ZIRCONIUM',41:'NIOBIUM',42:'MOLYBDENUM',43:'TECHNETIUM',
		  44:'RUTHENIUM',45:'RHODIUM',46:'PALLADIUM',47:'SILVER',
		  48:'CADMIUM',49:'INDIUM',50:'TIN',51:'ANTIMONY',
		  52:'TELLURIUM',53:'IODINE',54:'XENON'}

periodic_symbols = {'H':'HYDROGEN', 'HE':'HELIUM', 'LI':'LITHIUM', 'BE':'BERYLLIUM',
		  'B':'BORON', 'C':'CARBON', 'N':'NITROGEN', 'O':'OXYGEN',
		  'F':'FLUORINE', 'NE':'NEON', 'NA':'SODIUM', 'MG':'MAGNESIUM',
		  'AL':'ALUMINUM', 'SI':'SILICON', 'P':'PHOSPHOROUS',
		  'S':'SULFUR', 'CL':'CHLORINE', 'AR':'ARGON', 'K':'POTASSIUM',
		  'CA':'CALCIUM', 'SC':'SCANDIUM', 'TI':'TITANIUM', 'V':'VANADIUM',
		  'CR':'CHROMIUM', 'MN':'MANGANESE','FE':'IRON','CO':'COBALT',
		  'NI':'NICKEL','CU':'COPPER','ZN':'ZINC','GA':'GALLIUM',
		  'GE':'GERMANIUM','AS':'ARSENIC','SE':'SELENIUM','BR':'BROMINE',
		  'KR':'KRYPTON','RB':'RUBIDIUM','SR':'STRONTIUM','Y':'YTTRIUM',
		  'ZR':'ZIRCONIUM','NB':'NIOBIUM','MO':'MOLYBDENUM','TC':'TECHNETIUM',
		  'RU':'RUTHENIUM','RH':'RHODIUM','PD':'PALLADIUM','AG':'SILVER',
		  'CD':'CADMIUM','IN':'INDIUM','SN':'TIN','SB':'ANTIMONY',
		  'TE':'TELLURIUM','I':'IODINE','XE':'XENON'}
# Add 'CARBON':'CARBON' etc to periodic_symbols
for element in periodic_table.values():
                  periodic_symbols[element] = element

basis_list = ['3-21G','6-31G','6-31Gs','6-31Gss','6-31+G','6-31+Gs','6-31+Gss',
              '6-31++G','6-31++Gs','6-31++Gss','6-311G','6-311Gs','6-311Gss',
	      '6-311+G','6-311+Gs','6-311+Gss','6-311++G','6-311++Gs','6-311++Gss',
	      '6-31G_2d,2p','6-31+G_2d,2p','6-31++G_2d,2p',
	      '6-31G_3df,3pd','6-31+G_3df,3pd','6-31++G_3df,3pd',
	      '6-311G_2d,2p','6-311+G_2d,2p','6-311++G_2d,2p',
	      '6-311G_3df,3pd','6-311+G_3df,3pd','6-311++G_3df,3pd',
	      'cc-pVDZ','cc-pVTZ','cc-pVQZ','cc-pV5Z',
	      'aug-cc-pVDZ','aug-cc-pVTZ','aug-cc-pVQZ','aug-cc-pV5Z',
	      'ccd','cct','ccq','cc5','accd','acct','accq','acc5',
	      'aug-pc-0','aug-pc-1','aug-pc-2','aug-pc-3','aug-pc-4',
              'pc-0','pc-1','pc-2','pc-3','pc-4',
	      'pc1min','pc2min','pc3min','pc4min',
	      'apc0','apc1','apc2','apc3','apc4']
bare_basis_list = []
for i in range(0,len(basis_list)):
              bare = basis_list[i].replace("-","").replace(",","").replace("_","").upper()
	      bare_basis_list.append(bare)

small_basis_list = ['STO-3G','STO-3G-EMSL','STO-6G','MINI']
bare_small_basis_list = []
for i in range(0,len(small_basis_list)):
              bare = small_basis_list[i].replace("-","").replace(",","").replace("_","").upper()
	      bare_small_basis_list.append(bare)
basis_list.extend(small_basis_list)
bare_basis_list.extend(bare_small_basis_list)
#print bare_basis_list

#--------------------------------------------------------------------------
# Handy subroutines
#--------------------------------------------------------------------------

def get_contents(file):
    f = open(file,'r')
    fc = f.readlines()
    f.close()
    return fc

def extract_atom_basis(data,small_basis):
    geom = []
    for i in range(0,len(data)):
	if string.find(data[i],'END') == -1:
	    geom.append(data[i])
        else:
	    break
    symmetry = data[2].split()[0]
    geom.pop(0)
    geom.pop(0)
    geom.pop(0)
    atom_basis = []
    if small_basis != None:
        basis_string = small_basis
        atoms = []
        for i in range(0,len(data)):
            if string.find(data[i],'POPULATION ANALYSIS') != -1:
                atoms_start = i+1
            if string.find(data[i],'MOMENTS AT POINT') != -1:
	        atoms_stop = i 
        for i in range(atoms_start,atoms_stop):
            atom = data[i].split()[0].upper()
	    if atom not in list(periodic_symbols):
		print 'Error: atomic symbol ', atom, ' not recognised'
		print 'Please edit the POPULATION ANALYSIS section of your GAMESS'
		print 'dat file and $DATA section of your GAMESS input files'
		print 'Acceptable symbols are (case-insensitive)'
		print list(periodic_symbols)
		sys.exit()
	    else:
	        atom_basis.append([periodic_symbols[atom],basis_string])
    elif symmetry == 'C1':
        natom = len(geom)/3
        for i in range(0,natom):
	    atom_num = int(float(geom[3*i].split()[1]))
            basis = geom[3*i+1].split()
	    atom = periodic_table[atom_num]
	    if small_basis == None:
	        if len(basis) == 2:
	            basis_string = basis[0] + basis[1] + "G"
	            if basis_string != ('STO3G' or 'STO6G'):
		        print 'Warning: unless you are fitting to a minimal'
		        print 'basis your virtual orbitals will be rubbish'
	        elif len(basis) == 1:
	            basis_string = basis[0]
	            if basis_string != 'MINI':
		        print 'Warning: unless you are fitting to a minimal'
		        print 'basis your virtual orbitals will be rubbish'
	        else:
	            print 'Basis set not recognised from GAMESS dat file'
		    print 'If you want to use a Pople basis with extra diffuse'
		    print 'or polarization functions as your non-minimal small basis'
		    print 'you will need to specify its name on the command line'
	            sys.exit()
	    else:
	        basis_string = small_basis
	    atom_basis.append([atom,basis_string])
    else:
        print 'Sorry: I am not smart enough to detect both molecular symmetry'
        print 'and atomic basis sets concurrently. Please copy your full'
        print 'geometry from your output file into your dat file in C1 format'
        print 'or specify your universal small basis on the command line'
	sys.exit()
    return atom_basis

def extract_vec(data):
    have_uhf = 0
    for i in range(0,len(data)):
	if string.find(data[i],'UHF') != -1:
	    have_uhf = 1
    invec = 0
    for i in range(0,len(data)):
	if string.find(data[i],'VEC') != -1:
	    invec = 1
	    # keep last VEC section from dat file (clear earlier ones)
	    vec_raw_data = []
	    continue
	if string.find(data[i],'END') == -1 and invec == 1:
	    vec_raw_data.append(data[i])
	else:
	    invec = 0
    ivec = ' 0'
    allvecs = []
    vec = []
# need to insert end flag so processing code below
# knows when final MO has ended
    vec_raw_data.append('END')
    for i in range(0,len(vec_raw_data)):
	line = vec_raw_data[i]
	ii = line[0:2]
	c1 = line[5:20]
	c2 = line[20:35]
	c3 = line[35:50]
	c4 = line[50:65]
	c5 = line[65:80]
	if ii != ivec:
# append last vec if we have moved on (appends empty vec initially)
	    allvecs.append(vec)
	    ivec = ii
	    vec = []
        if c1 != '\n':
	    if c1 != '':
                vec.append(c1)
        if c2 != '\n':
	    if c2 != '':
                vec.append(c2)
        if c3 != '\n':
	    if c3 != '':
                vec.append(c3)
        if c4 != '\n':
	    if c4 != '':
                vec.append(c4)
        if c5 != '\n':
	    if c5 != '':
                vec.append(c5)
    allvecs.pop(0)
# check that all vecs are the same length
# not necessarily the case if nvec > 99
# shouldn't be able to get here anymore
# now that GAMESS prints identifiers modulo 100 rather than **
    vec_length = len(allvecs[0])
    for i in range(1,len(allvecs)):
	if len(allvecs[i]) != vec_length:
	    print 'Error: not all MO coefficient vectors the same length'
	    print 'Most likely due to having more than 99 MOs'
	    print 'Please edit GAMESS dat file to assign unique identifiers'
	    print 'to orbitals 100 and above, then re-run'
	    sys.exit()
    if have_uhf:
        allvecs1 = allvecs[:len(allvecs)/2]
	allvecs2 = allvecs[len(allvecs)/2:]
	allvecs_set = [allvecs1,allvecs2]
    else:
	allvecs_set = [allvecs]
    return allvecs_set

#--------------------------------------------------------------------------
# Main program
#--------------------------------------------------------------------------

if len(sys.argv) == 3:
    dat_file = sys.argv[1]
    lrg_basis = sys.argv[2]
    large_basis = lrg_basis.replace("-","").replace(",","").replace("_","").upper()
    small_basis = ''
    stem = dat_file.split("/")[-1].split('.')[0]
    if large_basis not in bare_basis_list:
	print 'Error: you have chosen an unsupported large basis set (' + lrg_basis + ')'
	print 'Available basis sets are:'
	print basis_list
	sys.exit()
elif len(sys.argv) == 4:
    dat_file = sys.argv[1]
    lrg_basis = sys.argv[2]
    sml_basis = sys.argv[3]
    large_basis = lrg_basis.replace("-","").replace(",","").replace("_","").upper()
    small_basis = sml_basis.replace("-","").replace(",","").replace("_","").upper()
    stem = dat_file.split("/")[-1].split('.')[0]
    if large_basis not in bare_basis_list:
	print 'Error: you have chosen an unsupported large basis set (' + lrg_basis + ')'
	print 'Available basis sets are:'
	print basis_list
	sys.exit()
    if small_basis not in bare_small_basis_list:
	print 'Warning: unless you are fitting to a minimal basis set'
	print '         or MCSCF optimized orbitals'
	print '         your virtual orbitals will be rubbish'
        if small_basis not in bare_basis_list:
	    print 'Error: you have chosen an unsupported small basis set (' + sml_basis + ')'
	    print 'Available minimal basis sets are:'
	    print small_basis_list
	    print 'All available basis sets are:'
	    print basis_list
	    sys.exit()
else:
    print 'Usage: process_datfile.py dat_file_name large_basis (small_basis)'
    sys.exit()

data = get_contents(dat_file)
if small_basis == '':
    atom_basis = extract_atom_basis(data,None)
else:
    atom_basis = extract_atom_basis(data,small_basis)
vecs_set = extract_vec(data)

# Initialize array to hold sets of output MO vectors (1 set for RHF, 2 for UHF)
out_vecs_set = []

for n in range(0,len(vecs_set)):
    vecs = vecs_set[n]
    # initialize output movec array
    out_vecs = []
    for i in range(0,len(vecs)):
        out_vecs.append([])

    # sanity check, number of MO coefficients should be equal to number of basis functions
    min_vec = vecs[0] 
    n_bas = 0
    for i in range(0,len(atom_basis)):
	atom = atom_basis[i][0]
	small_basis = atom_basis[i][1]
	min_bas_dat = basis.get(atom,small_basis)
	min_bas_order = min_bas_dat.pop()
	for j in range(0,len(min_bas_order)):
	    type = min_bas_order[j]
	    if type == 'S':
		n_bas += 1
	    elif type == 'P':
		n_bas += 3
            elif type == 'D':
		n_bas += 6
	    elif type == 'F':
		n_bas += 10
            elif type == 'G':
		n_bas += 15
	    elif type == 'H':
		n_bas += 21
	    else:
	        print 'Error: your minimal basis set appears to contain'
	        print 'basis functions of higher angular momentum than h.'
	        sys.exit()
    if n_bas != len(min_vec):
        print 'Error: number of molecular orbital coefficients does not'
        print 'match number of basis functions. Please check that your'
        print 'basis set is correctly defined, and/or that you have supplied the'
	print 'correct dat file for the minimal basis set you have nominated'
	sys.exit()
	    
    # expand out MOs from small to large basis atom by atom
    for i in range(0,len(atom_basis)):
        atom = atom_basis[i][0]
        small_basis = atom_basis[i][1]
        min_bas_dat = basis.get(atom,small_basis)
        lrg_bas_dat = basis.get(atom,large_basis)
        min_bas_order = min_bas_dat.pop()
        lrg_bas_order = lrg_bas_dat.pop()
	if len(min_bas_order) > len(lrg_bas_order):
            print 'Warning: you are fitting a larger basis with a smaller one'
	    print '         for atom: ', atom
    # do fitting of large basis to small (RMS fit)
        [s_fit,p_fit,d_fit,f_fit,g_fit,h_fit] = fit.do(min_bas_dat,lrg_bas_dat)
    # for each minimal basis molecular orbital
        for j in range(0,len(vecs)):
	    min_vec = vecs[j]
    # initialize MO coefficients for large basis 
    # (treating each orbital type separately for now)
            lrg_vec = []
            for k in range(0,len(lrg_bas_dat)):
                if len(lrg_bas_dat[k]) != 0:
                    shell_vec = []
		    for l in range(0,len(lrg_bas_dat[k])):
# s orbital case
		        if k == 0:
			    shell_vec.append([0.0])
# p orbital case
		        elif k == 1:
			    shell_vec.append([0.0,0.0,0.0])
## d orbital case (Spherical)
#		        elif k == 2:
#			    shell_vec.append([0.0,0.0,0.0,0.0,0.0])
## f orbital case (Spherical)
#		        elif k == 3:
#			    shell_vec.append([0.0,0.0,0.0,0.0,0.0,0.0,0.0])
## g orbital case (Spherical)
#		        elif k == 4:
#			    shell_vec.append([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
## h orbital case (Spherical)
#		        elif k == 5:
#			    shell_vec.append([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
# d orbital case (Cartesian)
		        elif k == 2:
			    shell_vec.append([0.0,0.0,0.0,0.0,0.0,0.0])
# f orbital case (Cartesian)
		        elif k == 3:
			    shell_vec.append([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
# g orbital case (Cartesian)
		        elif k == 4:
			    shell_vec.append([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
# h orbital case (Cartesian)
		        elif k == 5:
			    shell_vec.append([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    	              0.0,0.0,0.0,0.0,0.0,0.0])
	        lrg_vec.append(shell_vec)
# make a copy of fitting coefficients that we can sequentially pop off 
# one copy required for each MO
	    s_fit_cp = s_fit[:]
	    p_fit_cp = p_fit[:]
	    d_fit_cp = d_fit[:]
	    f_fit_cp = f_fit[:]
	    g_fit_cp = g_fit[:]
	    h_fit_cp = h_fit[:]
# expand out MOs
# note that MOs are ordered by atom type then atomic orbital 1s,2s,2p,3s,3p etc
# as specified in min_bas_order array
# accumulate new contraction coefficients as product of old cc's & fitting coeffs
            for k in range(0,len(min_bas_order)):
                type = min_bas_order[k]
	        if type == 'S':
		    cs = float(min_vec.pop(0))
		    if (s_fit != []):
		        sf = s_fit_cp.pop(0)
		        for l in range(0,len(sf)):
		            lrg_vec[0][l][0] += cs * sf[l]
#	            else:
#			print 'Warning: you are fitting a larger basis with a smaller one'
	        elif type == 'P':
		    cx = float(min_vec.pop(0))
		    cy = float(min_vec.pop(0))
		    cz = float(min_vec.pop(0))
		    if (p_fit != []):
		        pf = p_fit_cp.pop(0)
		        for l in range(0,len(pf)):
		            lrg_vec[1][l][0] += cx * pf[l]
		            lrg_vec[1][l][1] += cy * pf[l]
		            lrg_vec[1][l][2] += cz * pf[l]
#	            else:
#			print 'Warning: you are fitting a larger basis with a smaller one'
	        elif type == 'D':
		    cd1 = float(min_vec.pop(0))
		    cd2 = float(min_vec.pop(0))
		    cd3 = float(min_vec.pop(0))
		    cd4 = float(min_vec.pop(0))
		    cd5 = float(min_vec.pop(0))
# extra Cartesian function/s
		    cd6 = float(min_vec.pop(0))
		    if (d_fit != []):
		        df = d_fit_cp.pop(0)
		        for l in range(0,len(df)):
		            lrg_vec[2][l][0] += cd1 * df[l]
		            lrg_vec[2][l][1] += cd2 * df[l]
		            lrg_vec[2][l][2] += cd3 * df[l]
		            lrg_vec[2][l][3] += cd4 * df[l]
		            lrg_vec[2][l][4] += cd5 * df[l]
# extra Cartesian function/s
		            lrg_vec[2][l][5] += cd6 * df[l]
#	            else:
#			print 'Warning: you are fitting a larger basis with a smaller one'
	        elif type == 'F':
		    cf1 = float(min_vec.pop(0))
		    cf2 = float(min_vec.pop(0))
		    cf3 = float(min_vec.pop(0))
		    cf4 = float(min_vec.pop(0))
		    cf5 = float(min_vec.pop(0))
		    cf6 = float(min_vec.pop(0))
		    cf7 = float(min_vec.pop(0))
# extra Cartesian function/s
		    cf8 = float(min_vec.pop(0))
		    cf9 = float(min_vec.pop(0))
		    cf10 = float(min_vec.pop(0))
		    if (f_fit != []):
		        ff = f_fit_cp.pop(0)
		        for l in range(0,len(ff)):
		            lrg_vec[3][l][0] += cf1 * ff[l]
		            lrg_vec[3][l][1] += cf2 * ff[l]
		            lrg_vec[3][l][2] += cf3 * ff[l]
		            lrg_vec[3][l][3] += cf4 * ff[l]
		            lrg_vec[3][l][4] += cf5 * ff[l]
		            lrg_vec[3][l][5] += cf6 * ff[l]
		            lrg_vec[3][l][6] += cf7 * ff[l]
		            lrg_vec[3][l][7] += cf8 * ff[l]
		            lrg_vec[3][l][8] += cf9 * ff[l]
		            lrg_vec[3][l][9] += cf10 * ff[l]
#	            else:
#			print 'Warning: you are fitting a larger basis with a smaller one'
	        elif type == 'G':
		    cg1 = float(min_vec.pop(0))
		    cg2 = float(min_vec.pop(0))
		    cg3 = float(min_vec.pop(0))
		    cg4 = float(min_vec.pop(0))
		    cg5 = float(min_vec.pop(0))
		    cg6 = float(min_vec.pop(0))
		    cg7 = float(min_vec.pop(0))
		    cg8 = float(min_vec.pop(0))
		    cg9 = float(min_vec.pop(0))
# extra Cartesian function/s
		    cg10 = float(min_vec.pop(0))
		    cg11 = float(min_vec.pop(0))
		    cg12 = float(min_vec.pop(0))
		    cg13 = float(min_vec.pop(0))
		    cg14 = float(min_vec.pop(0))
		    cg15 = float(min_vec.pop(0))
		    if (g_fit != []):
		        gf = g_fit_cp.pop(0)
		        for l in range(0,len(gf)):
		            lrg_vec[4][l][0] += cg1 * gf[l]
		            lrg_vec[4][l][1] += cg2 * gf[l]
		            lrg_vec[4][l][2] += cg3 * gf[l]
		            lrg_vec[4][l][3] += cg4 * gf[l]
		            lrg_vec[4][l][4] += cg5 * gf[l]
		            lrg_vec[4][l][5] += cg6 * gf[l]
		            lrg_vec[4][l][6] += cg7 * gf[l]
		            lrg_vec[4][l][7] += cg8 * gf[l]
		            lrg_vec[4][l][8] += cg9 * gf[l]
		            lrg_vec[4][l][9] += cg10 * gf[l]
		            lrg_vec[4][l][10] += cg11 * gf[l]
		            lrg_vec[4][l][11] += cg12 * gf[l]
		            lrg_vec[4][l][12] += cg13 * gf[l]
		            lrg_vec[4][l][13] += cg14 * gf[l]
		            lrg_vec[4][l][14] += cg15 * gf[l]
#	            else:
#			print 'Warning: you are fitting a larger basis with a smaller one'
	        elif type == 'H':
		    ch1 = float(min_vec.pop(0))
		    ch2 = float(min_vec.pop(0))
		    ch3 = float(min_vec.pop(0))
		    ch4 = float(min_vec.pop(0))
		    ch5 = float(min_vec.pop(0))
		    ch6 = float(min_vec.pop(0))
		    ch7 = float(min_vec.pop(0))
		    ch8 = float(min_vec.pop(0))
		    ch9 = float(min_vec.pop(0))
# extra Cartesian function/s
		    ch10 = float(min_vec.pop(0))
		    ch11 = float(min_vec.pop(0))
		    ch12 = float(min_vec.pop(0))
		    ch13 = float(min_vec.pop(0))
		    ch14 = float(min_vec.pop(0))
		    ch15 = float(min_vec.pop(0))
		    ch16 = float(min_vec.pop(0))
		    ch17 = float(min_vec.pop(0))
		    ch18 = float(min_vec.pop(0))
		    ch19 = float(min_vec.pop(0))
		    ch20 = float(min_vec.pop(0))
		    ch21 = float(min_vec.pop(0))
		    if (h_fit != []):
		        hf = h_fit_cp.pop(0)
		        for l in range(0,len(hf)):
		            lrg_vec[5][l][0] += ch1 * hf[l]
		            lrg_vec[5][l][1] += ch2 * hf[l]
		            lrg_vec[5][l][2] += ch3 * hf[l]
		            lrg_vec[5][l][3] += ch4 * hf[l]
		            lrg_vec[5][l][4] += ch5 * hf[l]
		            lrg_vec[5][l][5] += ch6 * hf[l]
		            lrg_vec[5][l][6] += ch7 * hf[l]
		            lrg_vec[5][l][7] += ch8 * hf[l]
		            lrg_vec[5][l][8] += ch9 * hf[l]
		            lrg_vec[5][l][9] += ch10 * hf[l]
		            lrg_vec[5][l][10] += ch11 * hf[l]
		            lrg_vec[5][l][11] += ch12 * hf[l]
		            lrg_vec[5][l][12] += ch13 * hf[l]
		            lrg_vec[5][l][13] += ch14 * hf[l]
		            lrg_vec[5][l][14] += ch15 * hf[l]
		            lrg_vec[5][l][15] += ch16 * hf[l]
		            lrg_vec[5][l][16] += ch17 * hf[l]
		            lrg_vec[5][l][17] += ch18 * hf[l]
		            lrg_vec[5][l][18] += ch19 * hf[l]
		            lrg_vec[5][l][19] += ch20 * hf[l]
		            lrg_vec[5][l][20] += ch21 * hf[l]
#	            else:
#			print 'Warning: you are fitting a larger basis with a smaller one'
	        else:
		    print 'Error: your minimal basis set appears to contain'
		    print 'basis functions of higher angular momentum than h.'
		    sys.exit()
# sort new contraction coefficients by orbital order, remove structure
# and append to out_vecs 
            out_vec = []
            for k in range(0,len(lrg_bas_order)):
                type = lrg_bas_order[k]
	        if type == 'S':
		    if len(lrg_vec[0]) == 1:
		        vec = lrg_vec[0][0]
		    else:
		        vec = lrg_vec[0].pop(0)
	        elif type == 'P':
		    if len(lrg_vec[1]) == 1:
		        vec = lrg_vec[1][0]
		    else:
		        vec = lrg_vec[1].pop(0)
	        elif type == 'D':
		    if len(lrg_vec[2]) == 1:
		        vec = lrg_vec[2][0]
		    else:
		        vec = lrg_vec[2].pop(0)
	        elif type == 'F':
		    if len(lrg_vec[3]) == 1:
		        vec = lrg_vec[3][0]
		    else:
		        vec = lrg_vec[3].pop(0)
	        elif type == 'G':
		    if len(lrg_vec[4]) == 1:
		        vec = lrg_vec[4][0]
		    else:
		        vec = lrg_vec[4].pop(0)
	        elif type == 'H':
		    if len(lrg_vec[5]) == 1:
		        vec = lrg_vec[5][0]
		    else:
		        vec = lrg_vec[5].pop(0)
	        out_vec.extend(vec)
            out_vecs[j].extend(out_vec)
    out_vecs_set.append(out_vecs)

# now have new fitted vecs in out_vecs, in correct atom and shell order
# format out_vecs in GAMESS style and print
output_set = []
for n in range(0,len(out_vecs_set)):
    out_vecs = out_vecs_set[n]
    output = []
    for i in range(0,len(out_vecs)):
#for i in range(0,1):
        vec = out_vecs[i]
        vec_len = float(len(vec))
        nlines = int(math.ceil(vec_len/5.0))
        for j in range(0,nlines):
	    outline = ''
	    i_out = (i+1) % 100
	    j_out = j+1
	    outline += '%2i' % i_out
	    outline += '%3i' % j_out
	    if j < nlines-1:
                ccs = vec[5*j:5*j+5]
	    if j == nlines-1:
	        ccs = vec[5*j:]
            for k in range(0,len(ccs)):
	        outline += '% 15.8E' % ccs[k]
	    outline += '\n'
	    output.append(outline)
    output_set.append(output)

if len(output_set) == 1:
    output = output_set[0]
    rhf_output = []
    rhf_output.append(' $VEC\n')
    rhf_output.extend(output)
    rhf_output.append(' $END\n')
    uhf_output = []
    uhf_output.append(' $VEC\n')
    uhf_output.extend(output)
    uhf_output.extend(output)
    uhf_output.append(' $END\n')
elif len(output_set) == 2:
    output1 = output_set[0]
    output2 = output_set[1]
    rhf_output = []
    rhf_output.append(' $VEC\n')
    rhf_output.extend(output1)
    rhf_output.append(' $END\n')
    uhf_output = []
    uhf_output.append(' $VEC\n')
    uhf_output.extend(output1)
    uhf_output.extend(output2)
    uhf_output.append(' $END\n')
else:
# Should not be able to get here
    print 'Error: you appear to have more than two different electron spins'
    sys.exit()

f = open(stem+'.vec_rhf_'+large_basis,'w')
f.writelines(' $GUESS GUESS=MOREAD NORB=' + str(len(vecs_set[0])) + ' $END\n')
f.writelines(rhf_output)
f.close()
f = open(stem+'.vec_uhf_'+large_basis,'w')
f.writelines(' $GUESS GUESS=MOREAD NORB=' + str(len(vecs_set[0])) + ' $END\n')
f.writelines(uhf_output)
f.close()
