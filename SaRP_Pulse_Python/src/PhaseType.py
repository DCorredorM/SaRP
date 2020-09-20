import math
import numpy as np
from operator import concat
from functools import reduce
import warnings
import matplotlib.pyplot as plt
from mpmath import *
#from scipy.linalg import 

class CustomError(Exception):
	pass


class PH:
	def __init__(self, T, alpha):
		self.T=matrix(T)
		self.alpha= matrix(alpha)
		if len(self.T)==0:
			self.a=matrix([[]])
		else:
			self.a=-self.T*ones(len(T),1)
	
	def cdf(self,x):
		#print(len(self.T))
		if len(self.T)>0:
			#print(self.T*x)
			eTT=expm(self.T*x)			
			one=ones(len(self.T),1)
			p=1-self.alpha.T*eTT*one
		else:		
			p=[1]
		return p[0]

	def pdf(self,x):
		eTT=expm(self.T*x,method="pade")
		p=self.alpha.T*eTT*self.a
		return p[0]
	def plot_pdf(self,max,color="lime"):
		xx = np.linspace(0, max, 50)
		y2 = list(map(self.pdf,xx))
		plt.plot(xx,y2, color=color,alpha=1)
		#plt.legend(['PH-fit'])
		#plt.show()

	def sum(self, v2):
		if len(self.T)==0:
			return v2
		else:
			gamma=matrix(self.alpha.rows+v2.alpha.rows,1)
			for i in range(self.alpha.rows):
				gamma[i]=self.alpha[i]
			for i in range(self.alpha.rows,self.alpha.rows+v2.alpha.rows):
				gamma[i]=v2.alpha[i-self.alpha.rows]*(1-sum(self.alpha))
			
			C=matrix(self.T.rows+v2.T.rows)
			for i in range(self.T.rows):
				for j in range(self.T.cols):
					C[i,j]=self.T[i,j]

			ab=self.a*v2.alpha.T
			for i in range(ab.rows):
				for j in range(ab.cols):
					C[i,j+self.T.cols]=ab[i,j]
			
			for i in range(v2.T.rows):
				for j in range(v2.T.cols):
					C[i+self.T.cols,j+self.T.rows]=v2.T[i,j]
			
			return(PH(C,gamma))
