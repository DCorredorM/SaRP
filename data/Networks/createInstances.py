'''
This file creates the insatnces for the given city
'''

import os,sys,networkx as nx,random as rnd, numpy as np,time
sys.path.append(os.path.abspath(f'../../SaRP_Pulse_Python/src'))
from pulse import *
from Fitter import *

global city,DG,pairs,tightness,alpha
cities='SiouxFalss','Chicago-Sketch','GoldCoast','Philadelphia','Sydney'

def creates_graphs():
	'''
	Creates a networkx graph for shortest path pourpuses.	
	'''
	global DG
	DG=nx.DiGraph()	
	##########################
	#Creates deterministic graph
	file1 = open(f'{city}/{city}_net.txt','r')
	for x in file1:
		x=x[:-1].split('\t')
		#print(x)
		if x[0]!='Tail':
			x=list(map(float,x))			
			DG.add_edge(int(x[0]),int(x[1]),Cost=x[2],Time=x[3],tmin=x[4])
	file1.close()

def createPairs():
	'''
	Creates possible (s,t) pairs
	'''
	global pairs	
	sp=dict(nx.shortest_path_length(DG,weight='Cost'))	
	pairs={i:max(sp[i].keys(),key=lambda x:sp[i][x]) for i in sp.keys()}	

def createPair(s,k=50):
	'''
	for a given node s returns a target node
	Args: 
		s(int): source node
		k(int): Optional for choosing the k-th farthest nodes from s randomly
	return:
		t(int): target node
	'''	
	sp=nx.shortest_path_length(DG,source=s,weight='Cost')
	longest=sorted(sp.keys(),key=lambda x:-sp[x])[:k]

	return rnd.choice(longest)

def createRandomInst(n,wb='w'):
	'''
	Creates n random instances for the SaRP
	Args:
		n(int): Number of insances to create
		wb(String): type of the file for printing the instances: 'r' read, 'w' erite...
	Return:
		None
	'''
	f=open(f'{city}/{city}_instances.txt',wb)
	f.write(f'{"#"*40}\n#Creatin {n} instances with tightness: {tightness}\n{"#"*40}\n')
	for i in range(n):
		s=rnd.choice(list(DG.nodes()))
		t=createPair(s)
		tMax=calcTMax(s,t)
		f.write(f'{s}\t{t}\t{tMax}\n')
	f.close()

def calcTMax1(s,t,timelimit=200):
	'''
	Computes an interesting TMax for the SaRP problem from s to t
	Args:
		s(int):Source node.
		t(int):Target node.
	Return:
		Tmax(double):Time budget for s,t pair.	
	'''	
	PG=pulse_graph(G=DG,T_max=None,source=s, target=t,tightness=1)
	pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()		
	
	T_max=PG.T_max
	T_max_star=T_max
	
	##############################################################################
	#To do: While the relaiability of the path from the CSP is greater than tightness, reduce the time constraint and find again this path.
	stop=False
	delta=0.8
	try:
		S_pat=nx.shortest_path(DG,s,t,'Time')
	except:
		stop=True

	time_creit=time.time()
	while not stop:									
		rho=eval_path(pulse_path,T_max)
		sp_rho=eval_path(S_pat,T_max)
		
		print("se tiene que , ",rho,sp_rho,T_max,S_pat)
		
		if rho!=0:
			T_max_star=T_max
		if rho<tightness and sp_rho>=inf_tightness:
			stop=True
			pri=True
		elif sp_rho>=inf_tightness:
			T_max=T_max*delta
			PG=pulse_graph(G=DG,T_max=T_max,source=s, target=t)
			pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()		
			T_max=PG.T_max
		
		elif sp_rho<=inf_tightness:
			T_max=T_max*(1/0.9)
			delta=delta +0.5*(1-delta)
			#T_max=T_max*delta
			PG=pulse_graph(G=DG,T_max=T_max,source=s, target=t)
			pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()		
			T_max=PG.T_max
		
		if time.time()-time_creit>=timelimit:
			stop=True
		
	
	return T_max_star

def calcTMax(s,t):
	'''
	Computes Time budget using the empirical distribution of the shortest time path
	Args:
		s(int):source node.
		t(int):target node.
	Return:
		TMax(double):Time budget for (s,t) pair.
	'''
	spT=nx.shortest_path(DG,s,t,'Time')
	spC=nx.shortest_path(DG,s,t,'Cost')
	
	TspT=alphaQuantile(list(zip(spT[:-1],spT[1:])),alpha)
	TspC=alphaQuantile(list(zip(spC[:-1],spC[1:])),alpha)	
	return TspT+(TspC-TspT)*tightness

def alphaQuantile(arcs,a):
	'''
	Computes the alpha-Quantile for the given arcs
	Args:
		arcs(list):Arcs to compute Quantile.
		a(double):Quantile to compute.
	Return:
		aQuantile(double):Alpha quantile for given arcs.
	'''
	data=readRealizations(arcs)
	tmin=sum(DG[i][j]['tmin'] for i,j in arcs)
	#print(list(data))
	
	#ph,loglike=get_PH_HE(list(data),10)
	#print('ph10 da:', ph.cdf(np.quantile(data,a)),f'\nT es: {tmin+np.quantile(data,a)}')

	return tmin+np.quantile(data,a)

def readRealizations(arcs):
	'''
	Reads the historical data for the given arcs
	Args:
		arcs(list): list of arcs to query
	Return:
		data(np.array: sum of the arc travel times
	'''
	f=open(f'{city}/Scenarios/scenTotal.txt')
	def readLine(l):
		l=l.replace('\n','').split('\t')
		i,j=tuple(map(int,l[:2]))
		data=list(map(float,l[-1].replace('[','').replace(']','').split(', ')))
		return i,j,np.array(data)
	
	i,j,data=readLine(f.readline())	
	data-=data

	for l in f:
		i,j,datai=readLine(l)
		if (i,j) in arcs:
			data+=datai
	
	return data

def eval_path(path1,pT_max):
	'''
	Computes the historical reliability of the given path 
	Args:
		path(list):Pth to evaluate
		pT_max(double): Time to compyte the reliability (P(t(p)<T))
	Return:
		reliability(double):Reliability of the given path (P(t(p)<T)).
	'''
	if len(path1)>0:
		pat=list(zip(path1[:-1],path1[1:]))
		data=readRealizations(pat)	
		
		
		tmin=sum(int(DG[i][j]['tmin']) for i,j in pat)
		ind=list(map(lambda x: x<=pT_max-tmin,data))
		

		prob= sum(ind)/len(data)
		
		if prob==0 and len(pat)!=0:
			plt.axvline(x=pT_max-tmin,ymin=0, ymax=1000,color="lime")
			plt.hist(data)
			print("El tmin calculado es ",tmin)
			#plt.show()
			plt.clf()
		
		if len(pat)==0:
			prob=0
		
		ph,loglike=get_PH_HE(list(data),5)
		print('El tmin es: ',tmin)
		print('ph10 da:', ph.cdf(pT_max-tmin))		
		#print('El valE es: ',ph.T)
		print('El valE es: ',np.mean(data))
		return prob		
	else:
		return 0
	

if __name__ == '__main__':
	########################################################
	'''
	Setting of some global parameters
	'''
	city=cities[1]
	creates_graphs()
	tightness=0.6
	#inf_tightness=0.85
	alpha=0.8
	########################################################
	rnd.seed(1)
	createRandomInst(n=40,wb='a')

	
	# p=[142,688]#, 698, 700, 812, 471, 470, 469, 468, 467, 457, 466, 465, 861, 888, 894, 348]
	# p=[i+1 for i in p]
	# tmax=calcTMax(s=p[0],t=p[-1])
	# print('Este es el T',tmax)
	# print('Con data es: ',eval_path(p,tmax))

	
	