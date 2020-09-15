from s_pulse_MC import *
from PhaseType import *
from Fitter import *
from pulse import *
import Distrs_Fitter as dfit
import scipy.stats as sts
import networkx as nx
import time
import numpy as np
import Data_handler as dh
import os
import time
import numpy as np
from mpmath import *
import scipy.linalg as sp




os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse\Santos')


def cargar_redes(n,k):
	for i in range(n, k+1):
		t=time.time()
		inst='network{}.txt'.format(i)
		print(inst)
		dh.create_graph_stoch(3,inst,False)
		print(inst, time.time()-t)

#cargar_redes(10, 12)

def correr_instancias(inst,tideness,alpha):
	f = open('File1.txt', "r")
	cont=0
	instancias={}
	for x in f:
		x=x[:-4]
		x=list(map(float,x.split(";")))
		if int(x[0]) in inst:			
			instancias[int(x[0])]={"nodes":x[1],"arcs":x[2],"t_min_cost":x[3],"min_time":x[4],"target":x[5],"T_max":(x[4]-x[5])*tideness +x[4]}
			print(x)
	f.close()
	
	for (i,inst) in instancias.items():
		sfile=' network{}.gpickle'.format(i)
		os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse\gpickle')
		SG=nx.readwrite.gpickle.read_gpickle(sfile)
		for a in SG.edges():			
			T=matrix(SG[a[0]][a[1]]["T"])
			al=matrix(SG[a[0]][a[1]]["alpha"])
			#print(T.rows,T.cols)

			SG.edges[a]["tRV"]=PH(T,al)
		SPG=s_pulse_graph(G=SG,T_max=inst["T_max"],alpha=alpha,source=0,target=int(inst["target"]))		
		os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse\Santos')
		print(os.curdir)
		file= 'network{}.txt'.format(i)
		G=dh.create_graph_det(file, headder=False)
		PG=pulse_graph(G=G,T_max=inst["T_max"],source=0,target=int(inst["target"]))

		t=time.time()
		pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()
		pulse_time=time.time()-t

		print("Pulse")
		print("Time\t cost \t Time \t path")
		print(str(pulse_time), "\t",str(pulse_cost),"\t",str(pulse_sol_time),"\t",str(pulse_path))
		print("cota", PG.bound)
		print("Factibilidad", PG.infeas)
		print("Dominancia", PG.dom)

		t=time.time()
		#s_pulse_path,s_pulse_cost,s_pulse_prob=SPG.run_pulse()
		pat=zip(pulse_path[:-1],pulse_path[1:])
		ph=SPG.G[pat[0][0]][pat[0][1]]["tRV"]
		for i,j in pat[1:]:
			ph=ph.sum(SPG.G[pat[0][0]][pat[0][1]]["tRV"])

		print(ph.T)
		print(ph.alpha)
		text="{"
		for i in range(ph.alph.rows):
			text+=str(float(alpha[i]))+", "

		print(text,"}")
		for i in range(ph.T.rows):
			text="{"
			for i in range(ph.T.cols):
				text+=str(float(T[i,j]))+", "
			print(text,"}")

		s_pulse_time=time.time()-t
 
		print("S_Pulse")
		print("Time\t cost \t Probability \t  path")
		print(str(s_pulse_time), "\t",str(s_pulse_cost),"\t",str(s_pulse_prob),"\t",str(s_pulse_path))
		print("tiempo expm", SPG.time_expm, " segundos")
		print("cota", SPG.bound)
		print("Factibilidad", SPG.feas)
		print("Dominancia", SPG.dom)


		if s_pulse_cost!= pulse_cost and s_pulse_cost!=np.inf:
			res = open('Results.txt', "a")
			text=str(i)+"\t"+str(tideness)+"\t"+str(alpha)+"\t"+str(pulse_cost)+"\t"+str(s_pulse_cost)+"\t"+str(s_pulse_prob)+"\t"+str(pulse_time)+"\t"+str(s_pulse_time)+"\t"+str(SPG.bound)+"\t"+str(SPG.infeas)+"\t"+str(SPG.dom)+"\t"+str(SPG.time_expm)
			res.write(text)
			res.write("\n")

#correr_instancias([11,12],0.3,0.99)


def plot_distrs(G,edges=[]):
	if edges==[]:
		edges=G.edges()

	for x in edges:
		xx = np.linspace(0, 1000, 100)
		y2 = list(map(G.edges[(x[0],x[1])]["tRV"].pdf,xx))
		plt.plot(xx,y2, color="lime")
		plt.legend(['PH-fit','Data'])


		#print(G.edges[(x[0],x[1])]["tmin"])
		#print(G.edges[(x[0],x[1])]["tRV"].T)
	plt.show()

def plot_sum(G,path):
	PH=G.edges[path[0]]["tRV"]
	for i in path[1:]:
		PH=PH.sum(G.edges[i]["tRV"])
	PH.plot_pdf()


def prueba():
	#correr_instancias(50, 50,0.9,0.95)
	os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse\gpickle')
	sfile=' network{}.gpickle'.format(12)
	SG=nx.readwrite.gpickle.read_gpickle(sfile)
	for a in SG.edges():			
			T=matrix(SG[a[0]][a[1]]["T"])
			al=matrix(SG[a[0]][a[1]]["alpha"])
			SG.edges[a]["tRV"]=PH(T,al)
	t=6500
	T_max=(9085-8458)*0.1+8458
	alpha=.95
	SPG=s_pulse_graph(G=SG,T_max=T_max,alpha=alpha,source=0,target=t)
	s_pulse_path,s_pulse_cost,s_pulse_prob=SPG.run_pulse()
	print(s_pulse_path)
	print(list(zip(s_pulse_path[:-1],s_pulse_path[1:])))
	path=list(zip(s_pulse_path[:-1],s_pulse_path[1:]))
	#os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse\Santos')


	#plot_distrs(SPG.G,path)
	for i in range(1,len(path[1:])):		
		
		print(SPG.G.edges[path[i]]["tRV"].T)
		print(SPG.G.edges[path[i]]["tRV"].alpha)
		#print("p")
		#SPG.G.edges[path[i]]["tRV"].plot_pdf()
		#print("s")
		#plot_sum(SPG.G,path[:i])

#prueba()



	#SPG.tRV.plot_pdf()
	#cargar_redes(30, 100)
def is_diag(T):
	eye=np.ones(len(T))-np.eye(len(T))
	T=np.multiply(eye,T)
	if sum(T)==0:
		return True
	else:
		return False


#prueba()
def prueba2():
	x=[100, 200, 123, 343, 456, 123, 234, 312, 341, 123, 122, 234, 124, 782,190,293,320,201,123,123,134,432,412,321,232,123,142,143,145,231,353,123,112]
	print(sum(x))
	datos=[]
	var=sts.uniform.rvs(loc=0.2,scale=1)
	data = np.random.lognormal(log(x[0])-1/2*var, var**(1/2), 100)
	tmin=min(data)
	data=data-tmin
	datos.append(data)
	nPhases=3
	PH=get_PH_HE(data,nPhases)
	#PH.plot_pdf()
	for i in range(1,len(x[1:])):
		var=0.5
		data = np.random.lognormal(log(x[0])-1/2*var, var**(1/2), 100)
		tmin=min(data)
		data=data-tmin
		datos.append(data)
		nPhases=3
		ph=get_PH_HE(data,nPhases)

		'''
		#ph.plot_pdf()
		xx = np.linspace(000, 4000, 5)
		y2 = list(map(PH.pdf,xx))
		plt.hist(x=data,density=True,bins='auto', color='#0504aa',alpha=0.7, rwidth=0.5)
		plt.plot(xx,y2, color="lime")
		plt.legend(['PH-fit'])
		
		print(ph.T)
		print(ph.alpha)
		'''
		PH=PH.sum(ph)


		'''
		#print(datos[0]+datos[1])
		#print(sum(datos[:i]))
		plt.hist(x=sum(datos[k] for k in range(i+1)),density=True,bins='auto', color='#0504aa',alpha=0.7, rwidth=0.5)
		#lt.show()
		xx = np.linspace(000, 4000, 100)
		y2 = list(map(PH.pdf,xx))
		plt.plot(xx,y2, color="lime")
		plt.legend(['PH-fit'])
		#print(G.edges[(x[0],x[1])]["tmin"])
		'''

	'''
	Java	
	text="{"
	for i in range(PH.alpha.rows):
		text+=", "+str(float(PH.alpha[i]))

	print(text,"}")
	for i in range(PH.T.rows):
		text="{"
		cont=0
		for j in range(PH.T.cols):
			if cont==0:
				text+=str(float(PH.T[i,j]))
			else:
				text+=", "+str(float(PH.T[i,j]))
			cont+=1
		print(text,"},")
	
	print(PH.T.rows)
	print(PH.alpha.rows)
	
	t=time.time()
	print(PH.cdf(7700))
	print("calculando cfd con ",str(PH.T.rows)," fases me tardo ",time.time() -t,"segundos")
	'''
	'''
	#Julia
	text="["
	for i in range(PH.alpha.rows):
		text+=" "+str(float(PH.alpha[i]))

	print(text,"]")
	for i in range(PH.T.rows):
		text=""
		cont=0
		for j in range(PH.T.cols):
			if cont==0:
				text+=str(float(PH.T[i,j]))
			else:
				text+=" "+str(float(PH.T[i,j]))
			cont+=1
		print(text,";")
	
	print(PH.T.rows)
	print(PH.alpha.rows)
	'''
	#x=list(map(float,x.split(";")))
	#for i  in PH.T.tolist():
	#	print(list(map(float,i)))

	T=np.array([list(map(float,i)) for i  in PH.T.tolist()])
	alpha=np.array([list(map(float,i)) for i  in PH.alpha.tolist()])
	print(T)
	

	print(sp.expm(T)*7700)
	print(expm(PH.T*7700))
	eTT=sp.expm(T)*7700			
	one=np.ones((len(T),1))
	p=[[1]]-np.matmul(np.matmul(alpha.transpose(),eTT),one)
	print(p)

	#print(PH.T.tolist())

	t=time.time()
	print(PH.cdf(7700))
	print("calculando cfd con ",str(PH.T.rows)," fases me tardo ",time.time() -t,"segundos")


def cargar_toy():
	os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse\Notebooks')
	inst="toy.txt"
	SG=dh.create_graph_stoch(3,inst,True)
	for a in SG.edges():			
			T=matrix(SG[a[0]][a[1]]["T"])
			al=matrix(SG[a[0]][a[1]]["alpha"])
			#print(T.rows,T.cols)

			SG.edges[a]["tRV"]=PH(T,al)
	SPG=s_pulse_graph(G=SG,T_max=120,alpha=0.95,source=1,target=7)
	print(SPG.run_pulse())


'''
inst=[12,23,24,28,29, 36,49]
correr_instancias(inst, 0.7, .9)

#data=distrs["lnorm2"][0]
data=sts.erlang.rvs(160,10,size=100)
#data=min(data)
data=data-min(data)
print(min(data))

#data
plt.hist(x=data,density=True,bins='auto', color='#0504aa',alpha=0.7, rwidth=0.5)
xx = np.linspace(0, max(data), 100)


v=get_PH_HE(data,3)
y2 = list(map(v.pdf,xx))
plt.plot(xx,y2, color="lime")


v=get_PH_HE(data,30)
y2 = list(map(v.pdf,xx))
plt.plot(xx,y2, color="purple")


plt.legend(['PH-fit 3 Phases','PH-fit 30 Phases','Data'])
plt.xlim(0, max(data))
plt.xlabel("Time")
plt.ylabel("Density")
plt.show()
'''


def load_and_fit_net():
	G,pos=dfit.read("fitted_net")
	dist_names = ['laplace','norm','lognorm', 'powerlognorm' ,'exponweib', 'pareto','johnsonsu','burr']
	G=dfit.fit_graph(G,'fitted_net',dist_names)
	return G,pos
def load_net():

	G,pos=dfit.read("fitted_net")

	os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\3. Case_study\Chicago')
	G=nx.readwrite.gpickle.read_gpickle("fitted_net.gpickle")
	return G,pos






G,pos=dfit.read("fitted_net")

DG=nx.DiGraph()
for i,j in G.edges():
	DG.add_edge(i,j,Cost=G[i][j]["cost"],Time=G[i][j]["time"])


#nx.draw_networkx(DG,pos,node_size=10,node_color='black',with_labels=False,arrows=False)
#plt.show()





G,pos =load_net()
G=nx.readwrite.gpickle.read_gpickle("fitted_net.gpickle")

for i,j in G.edges():
	G[i][j]["Cost"]=DG[i][j]["Cost"]



alpha=0.95
s=random.randint(1,int(len(G.nodes())/2))
t=random.randint(int(len(G.nodes())/2),len(G.nodes()))
tightness=0.3

PG=pulse_graph(G=DG,T_max=0,source=s, target=t,tightness=tightness,pos=pos)
pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()




print("Pulse")
print("cost \t Time \t path")
print(str(pulse_cost),"\t",str(pulse_sol_time),"\t",str(pulse_path))
print("T_max", PG.T_max)
print("cota", PG.bound)
print("Factibilidad", PG.infeas)
print("Dominancia", PG.dom)



SPG=s_pulse_graph(G,T_max=PG.T_max,alpha=.95,source=s,target=t,n_it=1000,tightness=tightness,pos=pos)
path=list(zip(pulse_path[:-1],pulse_path[1:]))
tmin=sum(SPG.G[i][j]["tmin"] for i,j in path)
print("Probabilidad de llegar a tiempo: ", SPG.monte_carlo(pulse_path,tmin,0,10000))

pulse_path,pulse_cost,pulse_sol_time=SPG.run_pulse()
print("S_Pulse")
print("cost \t Time \t path")
print(str(pulse_cost),"\t",str(pulse_sol_time),"\t",str(pulse_path))
print("T_max", SPG.T_max)
print("cota", SPG.bound)
print("Factibilidad", SPG.infeas)
print("Dominancia", SPG.dom)
path=list(zip(pulse_path[:-1],pulse_path[1:]))
tmin=sum(SPG.G[i][j]["tmin"] for i,j in path)
print("Probabilidad de llegar a tiempo: ", SPG.monte_carlo(pulse_path,tmin,0,10000))



PG.draw_graph(path2=pulse_path)

#p=[899, 897, 442, 896, 891, 884, 856, 855, 854, 877, 874, 873, 870, 867, 864, 732, 414, 726, 720, 714, 390, 391, 392, 393, 583, 767, 756, 745, 199]
#p=[i+1 for i in p]

#PG.draw_graph(path2=p)





