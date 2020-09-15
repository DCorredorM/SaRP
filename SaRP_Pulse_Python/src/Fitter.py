from PhaseType import *
from numpy import *
from py4j.java_gateway import JavaGateway, GatewayParameters
import numpy as np
from py4j.java_collections import SetConverter, MapConverter, ListConverter,JavaArray
import scipy.stats as sts
from mpmath import *


def get_PH_HE (data,N):
	gateway = JavaGateway()
	fitter = gateway.entry_point.getEMHyperErlangFit()
	data=ListConverter().convert(data,gateway._gateway_client)
	emcont_fitter=fitter.EMHyperErlangFit(data)
	v=emcont_fitter.fit(N)
	#v=fitter.EMHyperErlangFit(data).fit(N)
	T=v.getMatrixArray()
	al=v.getVectorArray()
	loglike=emcont_fitter.getLogLikelihood()
	R=np.zeros((len(T),len(T)))
	R=array([[T[i][j]for j in range (len(T))] for i in range(len(T))])
	T=R
	alpha=array([al[i] for i in range(len(al))])

	T=matrix(T)
	alpha=matrix(alpha)
	return PH(T, alpha) ,loglike

def get_PH_MM (moments):
	gateway = JavaGateway()
	fitter = gateway.entry_point.getACPHFit()
	moments=ListConverter().convert(moments,gateway._gateway_client)
	v=fitter.MomentsACPHFit(moments).fit()
	T=v.getMatrixArray()
	al=v.getVectorArray()
	R=np.zeros((len(T),len(T)))
	R=array([[T[i][j]for j in range (len(T))] for i in range(len(T))])
	T=R
	alpha=array([al[i] for i in range(len(al))])
	return PH(T, alpha) 


#ph=get_PH_MM([2.0,6.0,25.0])
#print(ph.T)