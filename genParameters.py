#!/usr/bin/env python
import csv, subprocess
import numpy as np
import time
# import matplotlib.pyplot as plt

floatOut='%.4f'
intOut='%i'
outputFormat=[floatOut,floatOut,floatOut,floatOut,floatOut,floatOut,floatOut,intOut,intOut,intOut]
 








## variable parameter

wc_seed = [1.42]
inveps_seed = [1]


eps_ss_seed = [8]
wc_ss_seed = [1]
eps_sm_seed = [4]
wc_sm_seed = [1]
eps_bm_seed = [4]
wc_bm_seed = [1]


len_out = len(wc_seed)*len(inveps_seed)*len(eps_ss_seed)*len(wc_ss_seed)*len(eps_sm_seed)*len(wc_sm_seed)*len(eps_bm_seed)*len(wc_bm_seed)


realsMax = 1

len_out_master = len_out*realsMax

outarray=np.zeros(shape=(len_out,len(outputFormat)))
outarrayMaster=np.zeros(shape=(len_out_master,len(outputFormat)))

counterMaster = 0
for real in range(0,realsMax):
	counter = 0
	for i in range(len(wc_seed)):
		for j in range(len(inveps_seed)):
			for k in range(len(eps_ss_seed)):
				for l in range(len(wc_ss_seed)):
					for m in range(len(eps_sm_seed)):
						for n in range(len(wc_sm_seed)):
							for o in range(len(eps_bm_seed)):
								for p in range(len(wc_bm_seed)):

									seedBuf = np.random.randint(1000000)
									assert(seedBuf!=0)
									outarray[counter,0] = wc_seed[i]
									outarray[counter,1] = inveps_seed[j]
									outarray[counter,2] = eps_ss_seed[k]
									outarray[counter,3] = wc_ss_seed[l]
									outarray[counter,4] = eps_sm_seed[m]
									outarray[counter,5] = wc_sm_seed[n]
									outarray[counter,6] = eps_bm_seed[o]
									outarray[counter,7] = wc_bm_seed[p]
									outarray[counter,8] = seedBuf
									outarray[counter,9] = real

									outarrayMaster[counterMaster,0] = outarray[counter,0]
									outarrayMaster[counterMaster,1] = outarray[counter,1]
									outarrayMaster[counterMaster,2] = outarray[counter,2]
									outarrayMaster[counterMaster,3] = outarray[counter,3]

									outarrayMaster[counterMaster,4] = outarray[counter,4]
									outarrayMaster[counterMaster,5] = outarray[counter,5]
									outarrayMaster[counterMaster,6] = outarray[counter,6]
									outarrayMaster[counterMaster,7] = outarray[counter,7]
									outarrayMaster[counterMaster,8] = outarray[counter,8]
									outarrayMaster[counterMaster,9] = outarray[counter,9]

									counter += 1
									counterMaster += 1


	
	np.savetxt("Input/Parameters/params_%i.csv"%real , outarray, delimiter=',',fmt=outputFormat)
	np.savetxt("Input/Parameters/params_%i.txt"%real, outarray,fmt=outputFormat,delimiter=' ')

np.savetxt("Input/Parameters/paramsMaster.csv" , outarrayMaster, delimiter=',',fmt=outputFormat)
np.savetxt("Input/Parameters/paramsMaster.txt", outarrayMaster,fmt=outputFormat,delimiter=' ')
print(counterMaster)
	









