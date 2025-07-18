import numpy as np 
import itertools as it
from numpy.linalg import inv
import struct as st 
import matplotlib.pyplot as plt 
from csv import writer
import os

def Calc_pseudo(arrayL, arrayH):
    """calc vacuum expectaion value for each time step - pseudoscalar""" 
    pscalar = complex(0,0)
    for x,y,z,c1,c2 in it.product(range(Ns), range(Ns), range(Ns), range(3), range(3)):
        l = np.matrix(arrayL[x,y,z,:2,c1,:2,c2], dtype=complex)
        h = np.matrix(arrayH[x,y,z,:,c1,:,c2], dtype=complex)
        S = l*h.H
        pscalar += np.real(S[0,0] + S[1,1]) 
        
    return pscalar

def Calc_vector(arrayL, arrayH):
    """calc vacuum expectation value for each time step - vector """
    vector_val = 0
    for x,y,z,c1,c2 in it.product(range(Ns), range(Ns), range(Ns), range(3), range(3)):

        l = np.matrix(arrayL[x,y,z,:2,c1,:2,c2])
        h = np.matrix(arrayH[x,y,z,:,c1,:,c2])
        matr_1 = -S15*h*S15*l.H
        vector_val += np.real(matr_1[0,0] + matr_1[1,1])
        matr_2 = -S25*h*S25*l.H
        vector_val += np.real(matr_2[0,0] + matr_2[1,1])
        matr_3 = -S35*h*S35*l.H
        vector_val += np.real(matr_3[0,0] + matr_3[1,1])

    return vector_val

def Calc_correlator(filename_l, filename_h):
 ##import the propagators
    numdoub_l = (Ns**3)*4*4*3*3*2 
    numbyt_l = 8*numdoub_l

    numdoub_h = (Ns**3)*2*2*3*3*2
    numbyt_h = 8*numdoub_h

    with open(filename_l, 'rb') as light, open(filename_h, 'rb') as heavy:#read file and extract prop    
        #save 4 ints at start of file (light)  
        data = light.read(16)
        ints = st.unpack("iiii", data[0:16])
        #save norm
        data = light.read(8)
        norm = st.unpack("d", data[0:8])[0]
        
        #save filesize (heavy)
        data = heavy.read(4)
        size = st.unpack("i", data[0:4])
        
        prop_l = np.zeros((Ns, Ns, Ns, 4, 3, 4, 3), dtype=complex)
        prop_tran = np.zeros((Ns, Ns, Ns, 4, 3, 4, 3), dtype=complex)
        prop_h = np.zeros((Ns, Ns, Ns, 2, 3, 2, 3), dtype=complex)
        pseudo = np.zeros((Nt), dtype=complex)
        vector = np.zeros((Nt), dtype=complex)

        for t in range(Nt):

            data_l = light.read(numbyt_l)
            pointer_l = 0
            
            data_h = heavy.read(numbyt_h)
            pointer_h = 0
            
            #read in light prop
            for x,y,z,s1,c1,s2,c2 in it.product(range(Ns), range(Ns), range(Ns), range(4), range(3), range(4), range(3)):
                prop_l[x,y,z,s1,c1,s2,c2] = (st.unpack("d", data_l[pointer_l:pointer_l+8])[0] + st.unpack("d", data_l[pointer_l+8:pointer_l+16])[0]*1j)

                pointer_l += 16
            
            #change light prop chiral to NR rep
            for x,y,z,c1,c2 in it.product(range(Ns), range(Ns), range(Ns), range(3), range(3)):
                #prop_tran[x,y,z,:,c1,:,c2] = trans_matrix_inv*prop_l[x,y,z,:,c1,:,c2]*trans_matrix*(1/2)
                prop_tran[x,y,z,:,c1,:,c2] = trans_matrix*prop_l[x,y,z,:,c1,:,c2]*trans_matrix_inv


            #read in heavy prop
            for s1,s2,c1,c2,z,y,x in it.product(range(2), range(2), range(3), range(3), range(Ns), range(Ns), range(Ns)):
                prop_h[x,y,z,s1,c1,s2,c2] = (st.unpack("d", data_h[pointer_h:pointer_h+8])[0] + st.unpack("d", data_h[pointer_h+8:pointer_h+16])[0]*1j)

                pointer_h += 16

            pseudo[t] = Calc_pseudo(prop_tran, prop_h)
            vector[t] = Calc_vector(prop_tran, prop_h)

    light.close()
    heavy.close()

    return pseudo, vector

#Define matrices
trans_matrix = np.matrix([[1, 0, -1, 0],
                          [0, 1, 0, -1],
                          [1, 0, 1, 0],
                          [0, 1, 0, 1]], dtype = complex)
trans_matrix_inv = np.linalg.inv(trans_matrix)

Sigma0 = np.matrix([[1,0],
                    [0,1]], dtype=complex)
Sigma1 = np.matrix([[0,1],
                    [1,0]], dtype=complex)
Sigma2 = np.matrix([[0,-1j],
                    [1j,0]], dtype=complex)
Sigma3 = np.matrix([[1,0],
                    [0,-1]], dtype=complex)
Sigma5 = Sigma1*Sigma2*Sigma3*Sigma0
S15 = Sigma1*Sigma5
S25 = Sigma2*Sigma5
S35 = Sigma3*Sigma5

Gamma0 = np.matrix([[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, -1, 0],
                    [0, 0, 0, -1]], dtype=complex)
Gamma1 = np.matrix([[0, 0, 0, -1j],
                    [0, 0, -1j, 0],
                    [0, 1j, 0, 0],
                    [1j, 0, 0, 0]], dtype=complex)
Gamma2 = np.matrix([[0, 0, 0, -1],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [-1, 0, 0, 0]], dtype=complex)
Gamma3 = np.matrix([[0, 0, -1j, 0],
                    [0, 0, 0, 1j],
                    [1j, 0, 0, 0],
                    [0, -1j, 0, 0]], dtype=complex)
Gamma5 = Gamma1*Gamma2*Gamma3*Gamma0
G15 = Gamma1*Gamma5
G25 = Gamma2*Gamma5
G35 = Gamma3*Gamma5


######################################################################################

Ns=32
Nt=32

#light_dir = f'/home/postgrad/rachelhd/Bmeson/rqcd/propagators/{Ns}x{Nt}/strange/prop'
light_dir = f'/data/rachel/Bmeson/rqcd/propagators/{Ns}x{Nt}/s0'
heavy_dir = f'/data/rachel/Bmeson/nrqcd/propagators/{Ns}x{Nt}/prop/sprop'
file_dir =  f'/home/postgrad/rachelhd/Bmeson/hl_corrs/Gen2'

l_files = os.listdir(light_dir)    #list files in directory
h_files = os.listdir(heavy_dir)
l_files.sort()                     #sort because python doesn't
h_files.sort()
#print(l_files, h_files)
if len(l_files) > len(h_files):
    number_of_cfgs = len(l_files)
else:
    number_of_cfgs = len(h_files)
print(number_of_cfgs)
pseudo = []
vector = []

numcfgs=0
for i in range(6010,6020,10):
    if any(str(i) in x for x in l_files) and any(str(i) in x for x in h_files):
        ps, v = Calc_correlator(light_dir+f'/Gen2_{Ns}x{Nt}cfgn{i}.s0.m0',  heavy_dir+f'/sprop.{Ns}x{Nt}_{i}')
        cfg = np.array([i])
        ps = np.hstack((cfg, ps))
        v = np.hstack((cfg, v))
        pseudo.append(np.real(ps))
        vector.append(np.real(v))
        with open(file_dir+f'/hl_ps_{Ns}x{Nt}.csv', 'a') as file:
            write = writer(file)
            write.writerow(pseudo[numcfgs])
            file.close()

        with open(file_dir+f'/hl_v_{Ns}x{Nt}.csv', 'a') as file:
            write = writer(file)
            write.writerow(vector[numcfgs])
            file.close()
        numcfgs+=1
    else:
        print(f'no {i}')
<<<<<<< HEAD






=======
>>>>>>> c
