#http://www.bgu.ac.il/~bargera/tntp/
import sys
sys.path.insert(1, r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse')
from PhaseType import *
from Fitter import *
import networkx as nx
import os
import matplotlib.pyplot as plt
import math
import numpy as np
from s_pulse import *
from PhaseType import *
from Fitter import *
from pulse import *
import scipy.stats as sts
import networkx as nx
import time
import numpy as np
from mpmath import *


def read(instance):	
	os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\3. Case_study\Chicago\TransportationNetworks-master\Chicago-Sketch')
	G=nx.DiGraph()
	f=open('ChicagoSketch_node.tntp')
	pos={}
	for x in f:
		if x[:4]!='node':
			x=list(map(int,x[:-3].split("\t")))
			G.add_node(x[0])
			pos[x[0]]=(x[1],x[2])
	f.close()

	f=open('ChicagoSketch_net.tntp')
	for x in f:
		x=x[1:-3].split("\t")
		#print(x)
		t=True
		try:
			type(int(x[0]))
				
		except Exception as e:
			t=False
		if t:
			G.add_edge(int(x[0]),int(x[1]),cap=int(x[2]),leng=float(x[3]),fftt=float(x[4])*60,B=float(x[5]),pow=float(x[6]))

			
	f.close()
	f=open('ChicagoSketch_flow.tntp')
	for x in f:
		
		x=x.split('\t')
		
		t=True
		try:
			type(int(x[0]))
				
		except Exception as e:
			t=False
		if t:
			G[int(x[0])][int(x[1])]["cost"]=int(float(x[3][:-3])*100)
			G[int(x[0])][int(x[1])]["flow"]=int(float(x[2][:-1]))	
	f.close()
	edg=list(G.edges()).copy()
	for i,j in edg:
		if G[i][j]["fftt"]==0:
				G[i][j]["fftt"]=60*(60*G[i][j]["leng"]/50)#(mph))
		#Link travel time = free flow time * ( 1 + B * (flow/capacity)^Power ).**2
		G[i][j]["time"]=G[i][j]["fftt"]*(1+G[i][j]["B"]*(G[i][j]["flow"]/G[i][j]["cap"])**G[i][j]["pow"])
		G[i][j]["SD"]=max(1/2*(G[i][j]["time"]-G[i][j]["fftt"])*(G[i][j]["flow"]/G[i][j]["cap"]),5)

		G[i][j]["cost"]=int(100*G[i][j]["cost"]*(G[i][j]["fftt"]/(G[i][j]["time"])) )

		#print(G[i][j])
		'''
		if G[i][j]["SD"]<=0.001:
			G.remove_edge(i,j)
	nod=list(G.nodes()).copy()	
	for i in nod:
		if G.degree(i)==0:
			G.remove_node(i)
'''
	
	return(G,pos)


def write_txt(G,file,nPhases):
	os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\SPulseD\networks')
	file1 = open(file+'_fit_'+str(nPhases)+'.txt',"w") 
	cont=0
	for i,j in G.edges():
		phi=math.sqrt(G[i][j]["SD"]**2+ G[i][j]["time"]**2)
		mu=math.log(G[i][j]["time"]**2/phi)
		sigma=math.sqrt(math.log(phi**2/(G[i][j]["time"]**2)))
		data = np.random.lognormal(mu,sigma, 100)	
		'''
		if cont % 100==0:
			plt.hist(data)
			plt.show()
		cont+=1	
		'''
		tmin=min(data)
		data=list(map(lambda x:max(0.0,x),data-G[i][j]["fftt"]))
		
		
		ph=get_PH_HE(data,nPhases)
		G[i][j]["alpha"]=ph.alpha.tolist()
		G[i][j]["T"]=ph.T.tolist()
		

		text=str(i)+";"+str(j)+";"+str(int(G[i][j]["cost"]))+";"+str(G[i][j]["fftt"])+";"+str([list(map(lambda x: round(float(x),10),k))[0]for k in G[i][j]["alpha"]])+";"+str([list(map(lambda x: round(float(x),10),k))for k in G[i][j]["T"]])+'\n'
		print(text)
		
		file1.write(text)
	file1.close()
	nx.write_gpickle(G,file+"_"+str(nPhases)+'.gpickle')
sfile="ChicagoSketch_3.gpickle"
print(os.getcwd())



G,pos=read("p")
os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\SPulseD\networks')
G=nx.readwrite.gpickle.read_gpickle(sfile)
print(os.getcwd())
#write_txt(G,'ChicagoSketch',3)


###########
#Creates deterministic graph and creates and solves an instance of the deterministic pulse graph
DG=nx.DiGraph()
for i,j in G.edges():
	DG.add_edge(i,j,Cost=G[i][j]["cost"],Time=G[i][j]["time"])

SG=DG.copy()
for i,j in G.edges():
	T=matrix(G[i][j]["T"])
	al=matrix(G[i][j]["alpha"])
	SG[i][j]["tRV"]=PH(T,al)
	SG[i][j]["tmin"]=G[i][j]["fftt"]


nx.draw_networkx(DG,pos,node_size=10,node_color='black',with_labels=False,arrows=False)
plt.show()

alpha=0.95
source=900
target=200
tightness=0.5



PG=pulse_graph(G=DG,T_max=0,source=source, target=target,tightness=tightness,pos=pos)
pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()

SPG=s_pulse_graph(G=SG,T_max=PG.T_max,alpha=alpha,source=source,target=target)
#print(SPG.run_pulse())

print("Pulse")
print("cost \t Time \t path")
print(str(pulse_cost),"\t",str(pulse_sol_time),"\t",str(pulse_path))
print("T_max", PG.T_max)
print("cota", PG.bound)
print("Factibilidad", PG.infeas)
print("Dominancia", PG.dom)

p=[899, 897, 442, 896, 891, 884, 856, 855, 854, 877, 874, 873, 870, 867, 864, 732, 414, 726, 720, 714, 390, 391, 392, 393, 583, 767, 756, 745, 199]
p=[i+1 for i in p]
print(len(p))

'''
The final path is: [899, 897, 442, 896, 891, 884, 856, 855, 854, 877, 874, 873, 870, 867, 864, 732, 414, 726, 720, 714, 390, 391, 392, 393, 583, 767, 756, 745, 199]
The total cost is: 4237.0
The time limit is: 12271.0sec
The time limit is: 204.51666666666668min
The probability of arriving in time is: 0.9764589854383421
The computational time is: 3663.4220287 s
'''



#print("la prob es: ", SPG.prob_path( PG.Fpath,PG.T_max))
PG.draw_graph(path2=p)



#def montecarlo(n_it,path):






















