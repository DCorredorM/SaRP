import networkx as nx
from PhaseType import *
from Fitter import *
import scipy.stats as sts
import networkx as nx
import time
import numpy as np
import os


#Charges the graph from a txt file and creates a graph from networkx
#For the pulse algorithm
def create_graph_det(instancia,headder=True):
	G=nx.DiGraph()
	f = open(instancia, "r")
	cont=0
	if headder:
		for x in f:
			if cont >0:
				#print(x)	
				x=list(map(float,x.split("\t")))
				G.add_edge(int(x[0]),int(x[1]),Cost=x[2],Time=x[3])
				#print(x)
			cont+=1
	else:
		for x in f:
			x=list(map(float,x.split("\t")))
			G.add_edge(int(x[0]),int(x[1]),Cost=x[2],Time=x[3])

	return(G)

#Creates a graph from the .txt file instancia, fits a PH with nPhases to each arc and saves a .gpickle file with the networkx graph
def is_diag(T):
	eye=np.ones(len(T))-np.eye(len(T))
	T=np.multiply(eye,T)
	if sum(T)==0:
		return True
	else:
		return False

#Reads a Grph from instance "Dimacs format" and generates a PH for the times of each arc. 
def create_graph_stoch(nPhases,instancia,headder=True):
	G=nx.DiGraph()
	f = open(instancia, "r")
	nlfile = len(f.readlines())
	print(nlfile)
	f.close()
	f = open(instancia, "r")
	#nlfile=15000
	#print(nlfile)
	cont=0
	tiempos=[]
	#print(f)
	for x in f:
		#print(x)
		var=sts.uniform.rvs()*10
		
		if headder and  cont ==0:
			pass
		else:				
			t=time.time()
			x=list(map(float,x.split("\t")))
			
			phi=sqrt(var+x[3]**2)
			mu=math.log(x[3]**2/phi)
			sigma=sqrt(math.log(phi**2/(x[3]**2)))
			data = np.random.lognormal(mu,sigma, 100)
			
			#plt.hist(data)
			#plt.show()

			tmin=min(data)
			data=data-tmin
			ph=get_PH_HE(data,nPhases)
			G.add_edge(int(x[0]),int(x[1]),Cost=x[2],alpha=ph.alpha.tolist(),T=ph.T.tolist(),tmin=tmin)
			tiempos.append(time.time()-t)


			#print('promemdio fit: ', mean(tiempos))
			if cont % 1000==0:
				print("Voy en: ",cont)
				print("Me he demorado en promedio: ",mean(tiempos))
				print('pronostico restante: ', mean(tiempos)*(nlfile-cont))
		cont+=1

	os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\SPulseD\santosInstances\Santos')
	file1 = open(instancia[:-4]+'_fit'+'.txt',"w") 
	for i,j in G.edges():
		text=str(i)+";"+str(j)+";"+str(int(G[i][j]["Cost"]))+";"+str(G[i][j]["tmin"])+";"+str([list(map(lambda x: round(float(x),10),k))[0]for k in G[i][j]["alpha"]])+";"+str([list(map(lambda x: round(float(x),10),k))for k in G[i][j]["T"]])+'\n'
		print(text)
		file1.write(text)
	file1.close()
	nx.write_gpickle(G, r"C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse\gpickle"+"\ "+instancia[:-4]+'.gpickle')
	return(G)

