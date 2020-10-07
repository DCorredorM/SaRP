import os,time,random,sys
import numpy as np
import scipy.stats as st
import networkx as nx

sys.path.append(os.path.abspath('../../SaRP_Pulse_Python/src/'))

from pulse import *
import s_pulse_MC as MC

global city,timeLimit,CV,PG,SPG, alpha,n_it,refit

def creates_graphs():

	folder=os.path.abspath(f'../../data/Networks/{city}')
	
	
	DG=nx.DiGraph()
	#os.chdir(r'Networks/'+city)

	##########################
	#Creates deterministic graph
	file1 = open(f'{folder}/{city}_net.txt','r')
	for x in file1:
		x=x[:-1].split('\t')
		#print(x)
		if x[0]!='Tail':
			x=list(map(float,x))			
			DG.add_edge(int(x[0]),int(x[1]),Cost=x[2],Time=x[3],tmin=x[4])
	file1.close()

	##########################
	#Creates graph with distrs
	SG=DG.copy()
	file1 = open(f'{folder}/Independent/CV{CV}/DistrFit_cv{CV}.txt','r')
	
	for x in file1:
		x=x[:-1].split('\t')
		#print(x)
		if x[0]!='Tail':
			#print(x)
			params=tuple(map(float,x[-1][1:-1].split(', ')))
			#print(params)
			i,j=int(float(x[0])),int(float(x[1]))
			SG[i][j]["distr"]=x[2]
			SG[i][j]["params"]=params
	file1.close()
	
	global PG,SPG
	PG=pulse_graph(G=DG,T_max=0,source=0, target=0,tightness=0)
	SPG=MC.s_pulse_graph(SG,T_max=0,alpha=alpha,source=0,target=0,n_it=n_it,tightness=0)
	
	#return DG,SG

def solveCSP(s,t,T):
	'''
	Solves deterministic pulse for given instance
	Args:
		arg1(type):description.
	Return:
		Return(type):description.
	'''
	PG.source,PG.target,PG.T_max=s,t,T
	t_pulse=time.time()
	pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()
	t_pulse=time.time()-t_pulse

	prob=evalPath(path=pulse_path,nSim=n_it*10,pTMax=T)

	return t_pulse,pulse_cost,pulse_sol_time,prob,PG.bound, PG.infeas,PG.dom,pulse_path

def solveSaRPMC(s,t,T):
	'''
	Solves SaRP with montecarlo pulse for given instance
	Args:
		arg1(type):description.
	Return:
		Return(type):description.
	'''
	SPG.source,SPG.target,SPG.T_max=s,t,T
	t_pulse=time.time()	
	pulse_path,pulse_cost,pulse_prob=SPG.run_pulse(t_limit=1000)
	t_pulse=time.time()-t_pulse

	probPost=evalPath(path=pulse_path,nSim=n_it*10,pTMax=T)
	path=list(zip(pulse_path[:-1],pulse_path[1:]))
	Exp_time=sum(SPG.G[i][j]["Time"] for i,j in path)
	return t_pulse,pulse_cost,Exp_time,pulse_prob,probPost,PG.bound, PG.infeas,PG.dom,pulse_path

def createData(mu,fft,n):
	'''
	Given the expected value, and a CV, a dataset with n realizations of a lognormal variable with Expected value mu and with variance s=CV*mu
	Args:
		mu (float): Expected value
		
		fft (float): Free flow time 
		n (int): size of the dataset
	Returns
		data [list]: 
	'''
	sigma=CV*(mu-fft)
	phi=math.sqrt(sigma**2+ mu**2)
	muN=math.log(mu**2/phi)
	sigmaN=math.sqrt(math.log(phi**2/(mu**2)))
	
	data = np.random.lognormal(muN,sigmaN, n)		
	data=list(map(lambda x:max(0.0,x),data-fft))

	# print('Mean:\t',np.mean(data),mu-fft)
	# print('Stdev:\t',np.std(data),sigma)

	return np.array(data)

def evalPath(path,nSim,pTMax):
	'''
	Evaluates a given path with the real (assumed) distributions.
	Args:
		path(list):Path to evaluate
		nSim(int):Number of replications for montecarlo
		pTMax(float):Time budget to evaluate
	Return:
		prob(float):reliability of the path
	'''

	pat=list(zip(path[:-1],path[1:]))
	data=np.array([0.0 for i in range(nSim)])
	tmin=0
	for i,j in pat:
		data+=createData(mu=SPG.G[i][j]["Time"],fft=SPG.G[i][j]["tmin"],n=nSim)		
		tmin+=math.floor( SPG.G[i][j]["tmin"])

	ind=list(map(lambda x: x<=pTMax-tmin,data))
	prob= sum(ind)/nSim

	if len(pat)==0:
		prob=0
	return prob

def createFolder(name):	
	if os.name=='posix':
		os.popen(f'[ -d {name} ]  || mkdir {name}')
	else:
		name=name.replace('/','\\')
		os.popen(f'md {name}')

def runExperiment():
	'''
	Clas the java pulse for the config file located at f'../../data/Networks/{city}'
	Args:
		city(String): folder name with the data
	Return:
		None: prints out a file with the results:
			source+"\t"+target+"\t"+ running time (s)+"\t"+PrimalBound+"\t"+Free flow time+"\t"+ reliability+"\t"+Bound+"\t"+Infeasibility+"\t"+Dominance+"\t"+finalPath
	'''
	folder=os.path.abspath(f'../../data/Networks/{city}')
	d=os.popen(f'java -Xmx2048m -jar PH_Pulse_SaRP.jar {folder} {city}')
	d=d.read().replace('\n','').replace('pk menor que Double.MIN_VALUE','')#.split('\t')	
		

	return d

def runInstancesIndependent(nPhases,clearF=True):
	'''
	Runs all the insances for the given city
	Args:
		city(String): folder name with the data
	Return:
		None: prints out a file with the results:
			source+"\t"+target+"\t"+ running time (s)+"\t"+PrimalBound+"\t"+Free flow time+"\t"+ reliability(ex-ante)+"\t""+ reliability(ex-post)+"\t"+Bound+"\t"+Infeasibility+"\t"+Dominance+"\t"+finalPath
	'''
	#if clearF: clearResultFiles()
	
	folder=os.path.abspath(f'../../data/Networks/{city}')
	
	
	inst=open(f'{folder}/Independent/Instances/i.txt')
	#inst=open(f'{folder}/{city}_instances_{tightness}.txt','r')	
	results=f'../../results/{city}/Independent/alpha{alpha}_CV{CV}_tight{tightness}'
	createFolder(name=results)
		
	nn=0
	wb='w'
	for l in inst:
		if l[0]!='#':
			i=l.replace('\n','').split('\t')
			print(i)
			#Solve SaRP with ph
			s,t,tL,tU=int(i[0]),int(i[1]),float(i[2]),float(i[3])
			T=tL+(tU-tL)*(1-tightness)
			for np in nPhases:
				config=open(f'{folder}/{city}_config.txt','w')				
				text=f'PHFitFile:Independent/CV{CV}/PHFit{np}_cv{CV}.txt\nDataFile:Independent/CV{CV}/data_cv{CV}.txt\nNumber of Arcs:2950\nNumber of Nodes:933\nTime Constraint:{T}\nStart Node:{s}\nEnd Node:{t}\nNumber of Phases:{np}\nalpha:{alpha}\nTime Limit:{timeLimit}\nrefit:{refit}'
				config.write(text)				
				config.close()
				d=runExperiment()
				#Evals path
				info=d.split('\t')
				
				if info[-1]!='[]':
					path=list(map(lambda x: int(x)+1 ,info[-1].replace('[','').replace(']','').split(', ')))
					probPost=evalPath(path,nSim=n_it*10,pTMax=T)
				else:
					probPost=0
		
				t_pulse,pulse_cost,pulse_sol_time,probAnte,bound, infeas,dom,pulse_path=info[2:]
				d=f'{s}\t{t}\t{t_pulse}\t{pulse_cost}\t{pulse_sol_time}\t{probAnte}\t{probPost}\t{bound}\t{infeas}\t{dom}\t{pulse_path}\n'
				
				resultFile=open(f'{results}/PHFit{np}.txt',wb)				
				resultFile.write(d)
				print(f'PH{np}\t{d}')
				resultFile.close()
			
			s,t,T=int(i[0])+1,int(i[1])+1,float(i[2])
			#Solve CSP
			t_pulse,pulse_cost,pulse_sol_time,prob,bound, infeas,dom,pulse_path=solveCSP(s,t,T)
			d=f'{s}\t{t}\t{t_pulse}\t{pulse_cost}\t{pulse_sol_time}\t{prob}\t{prob}\t{bound}\t{infeas}\t{dom}\t{pulse_path}\n'
			resultFile=open(f'{results}/CSP.txt',wb)
			print(f'Pulse\t{d}')
			resultFile.write(d)
			
			#Solve SaRP with MC
			# t_pulse,pulse_cost,pulse_sol_time,probAnte,probPost,bound, infeas,dom,pulse_path=solveSaRPMC(s,t,T)
			# d=f'{s}\t{t}\t{t_pulse}\t{pulse_cost}\t{pulse_sol_time}\t{probAnte}\t{probPost}\t{bound}\t{infeas}\t{dom}\t{pulse_path}\n'
			# resultFile=open(f'{results}/MC_Pulse.txt',wb)
			# print(f'PulseMC\t{d}')
			# resultFile.write(d)
			wb='a'


			nn+=1
			if nn>=1000:
				break

def clearResultFiles():
	'''
	clears all result files from city path
	'''
	folder=os.path.abspath(f'../../data/Networks/{city}/Results/Independent/*.txt')
	f=open(folder[:-5]+'n.txt','w')
	f.write('i')
	f.close()

	if os.name=='posix':
		os.popen(f'rm {folder}')
	else:
		os.popen(f'del {folder}')


if __name__ == '__main__':
	########################################################
	'''
	Run scenario instances.

	Setting of some global parameters
	'''
	city='Chicago-Sketch'	
	timeLimit=5000	
	tightness=0.8
	CV=0.8
	alpha=0.9
	n_it=500	#Number of realizations for Montecarlo
	refit=3
	creates_graphs()
	########################################################		

	runInstancesIndependent(nPhases=[3,5])
	