import scipy.stats as sts
import numpy as np
import networkx as nx
import time

import matplotlib.pyplot as plt
import os

from norta import *
from functools import partial

#given a covariance matrix, a list of distributions, and a list with the parameters of each distributions generates n random samples of the random vector with those distributions and with tath covariance matrix.
def Cor_gen(cov,distrs, params,n):
	#print('Estoy por aca...')
	distrs=[getattr(sts, l) for l in distrs]
	t=time.time()
	Finv=[]
	def inv (d,par, x):
		return d.ppf(x, *par)

	for i,d in enumerate(distrs):
		par=params[i][:]    
		Finv.append(partial(inv,d,par))
	#print(cov)
	return NORTA(Finv, cov).gen(n=n)


'''

C=np.array([[1.,0.47155524, 0.4485085],
 [0.47155524, 1.,         0.67380935],
 [0.4485085,  0.67380935, 1.        ]])

C=np.eye(3)
print(C)
distrs=['lognorm','exponweib','burr']
params=[(2.374041293269499, -4.091759500327035e-25, 0.6057170920291481),(1.1935943026668159, 0.42895706224152175, -1.1535381722452478e-31, 0.9509420889710204),(165.91864571958797, 5826.322971606714, -4450.197844202376, 4254.579554419175)]
distrs=[getattr(sts, l) for l in distrs]


#print(help(partial))


Finv=[]
EX=[l.mean(*params[i]) for i,l in enumerate(distrs)]
print(EX)
SD=[l.var(*params[i])**(1/2) for i,l in enumerate(distrs)]
print(SD)

#scen=[sum(k) for k in dat]
#print(scen)
EX=[l.mean(*params[i]) for i,l in enumerate(distrs)]

SD=[l.var(*params[i])**(1/2) for i,l in enumerate(distrs)]

Finv=[]
def inv (d,par, x):
	return d.ppf(x, *par)

for i,d in enumerate(distrs):
	par=params[i][:]    
	Finv.append(partial(inv,d,par))



d=len(C)
CovX=np.zeros((d,d))
for i in range(d):
	for j in range(i,d):
		CovX[i,j]=C[i,j]*SD[i]*SD[j]
		CovX[j,i]=C[i,j]*SD[i]*SD[j]




norta_data = fit_NORTA_CovX(CovX=CovX, EX=EX, F_invs=Finv,output_flag=1)
NG = norta_data.gen(1000)
print(NG.mean(axis=0), EX)
print(np.corrcoef(NG, rowvar=0))

print(np.cov(NG, rowvar=False))
print(CovX)


'''

"""
Example of using NORTA
"""
'''
np.random.seed(0)
n_sample = 100
d_sample = 3
cov_sample = np.eye(d_sample) + np.random.rand(d_sample, d_sample)

sim_cov = cov_sample.transpose().dot(cov_sample)

print(cov_sample)
print(sim_cov)

data = np.random.exponential(size=(n_sample, d_sample)) + np.random.multivariate_normal(np.zeros(d_sample), sim_cov, size=n_sample)


n = len(data)
d = len(data[0])
norta_data = fit_NORTA(data, n, d-1, output_flag=1)
NG = norta_data.gen(1000)
print(NG.mean(axis=0), data.mean(axis=0))
print(np.corrcoef(NG, rowvar=0))
print(np.corrcoef(data, rowvar=0))
print(np.cov(NG, rowvar=False))
print(np.cov(data, rowvar=False))

'''




