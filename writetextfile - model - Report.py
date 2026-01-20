# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 22:01:39 2023

@author: lizad
"""

import numpy as np

n_runs = 4
constants = np.zeros((2, n_runs))

'''3 arms'''
#high peak at 10^4, low receptor depletion
constants[0,0] = 1e-12
constants[1,0] = 5e-3

#low peak at 10^4, high receptor depletion
constants[0,1] = 5e-8
constants[1,1] = 1e-4

#high peak at 10^5, low receptor depletion
constants[0,2] = 1e-11
constants[1,2] = 3e-4

# low peak at 10^5, high receptor depletion
constants[0,3] = 1e-7
constants[1,3] = 1e-5


'''6 arms'''
#high peak at 10^4, low receptor depletion
constants[0,0] = 1e-15
constants[1,0] = 5e-4

#low peak at 10^4, high receptor depletion (already done)
constants[0,1] = 1e-8
constants[1,1] = 5e-5

#high peak at 10^5, low receptor depletion
constants[0,2] = 1e-13
constants[1,2] = 4e-5

#low peak at 10^5, high receptor depletion
constants[0,3] = 1e-8
constants[1,3] = 5e-6

for i_run in range(n_runs):
    koffA = 1
    konsurfA = constants[1,i_run]
    konsolA = constants[0,i_run]

    koffB = 1
    konsurfB = 1e-5
    konsolB = 1e-9
        
    n_armsA = 6
    n_armsB = 0

    regel1 = '// Binding arm typa A'
    regel2 = 'Binding n00 receptA'
    regel3 = 'n00 -1 receptA 1 n01 1'
    regel4 = ''

    file_name = f'{n_armsA} arms konsol{konsolA: .2e} konsur{konsurfA: .2e}_model.txt'
    f = open(f'.\\BEP runs\\txt files\\{file_name}', "w")
    f.write("// Whitespace at the end of line (especially change in reactants) causes strange behaviour")
    f.write("// Make sure to remove all whitespace!!!\n")
    f.write(f"nanostar {n_armsA} arms model, konsol: {konsolA}, konsur: {konsurfA}, koff: {koffA}\n\n\
> Species, species dimension and diffuion constants\n\
// On membrane\n\
receptA         2       0\n\
receptB         2       0\n")
                    
    for ii in range(n_armsA+1):
        for jj in range(n_armsB+1):
            if ii==0 and jj==0:
                a = 1
            else:
                f.write("n"+str(ii)+str(jj)+"\t\t 2\t\t 0\n")

    f.write("\n\n\
// In solution\n\
n00\t\t 3\t 0\n\
receptora_tot\t 3\t 0\n\
receptorb_tot\t 3\t 0\n\
nanostars_tot\t 3\t 0\n\
coverage \t 3\t 0 \n\n\
> Reactions and rate\n")

    # Receptor A binding

    for ii in range(n_armsA):
        konrate = konsurfA*(n_armsA-ii)
        print("ii = "+ str(ii))
        print("konsurfA = "+ str(konsurfA))
        print("n_armsA = " + str(n_armsA))
        print("konsrufA*n_arms = " + str(konsurfA*n_armsA))
        print("konsurfA*ii = "+ str(konsurfA*ii))
        print("konsurfA*(n_armsA-ii) = " + str(konsurfA*(n_armsA-ii)))
        #print(konsurfA*ii-konsurfA*ii)
        for jj in range(n_armsB+1):
            f.write("// binding arm with type A\n\
binding n"+str(ii)+str(jj)+"\n\
n"+str(ii)+str(jj)+" -1 receptA -1 n"+str(ii+1)+str(jj)+" 1\n")
            if ii==0 and jj==0:
                f.write(str(konsolA*n_armsA)+"\n\n")
            else:
                f.write(str(konsurfA*(n_armsA-ii))+"\n\n")
                
            f.write("// unbinding arm with type A\n\
unbinding n"+str(ii+1)+str(jj)+"\n\
n"+str(ii+1)+str(jj)+" -1 receptA 1 n"+str(ii)+str(jj)+" 1\n"\
+str(koffA*(ii+1))+"\n\n")



    # Receptor B binding

    for ii in range(n_armsA+1):
        for jj in range(n_armsB):
            f.write("// binding arm with type B\n\
binding n"+str(ii)+str(jj)+"\n\
n"+str(ii)+str(jj)+" -1 receptB -1 n"+str(ii)+str(jj+1)+" 1\n")
            if ii==0 and jj==0:
                f.write(str(konsolB*n_armsB)+"\n\n")
            else:
                f.write(str(konsurfB*(n_armsB-jj))+"\n\n")

            f.write("// unbinding arm with type B\n\
unbinding n"+str(ii)+str(jj+1)+"\n\
n"+str(ii)+str(jj+1)+" -1 receptB 1 n"+str(ii)+str(jj)+" 1\n"\
+str(koffB*(jj+1))+"\n\n")

    f.close()