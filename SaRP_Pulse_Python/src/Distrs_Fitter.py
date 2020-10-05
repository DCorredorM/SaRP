import scipy.stats as st
import numpy as np
import networkx as  nx
import os,math,sys,time

from Fitter import *

global CV #CV (float):	Coefficient of variation

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
		# print(x)
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

		G[i][j]["cost"]=int(100*G[i][j]["cost"]*(G[i][j]["fftt"]/(G[i][j]["time"]**2)) )

	
	return(G,pos)

def fit_graph(G,file,d_names):
	os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\SPulseD\networks')
	cont=0
	NG=nx.DiGraph()
	for i,j in G.edges():
		#Creates data with the given moments
		phi=math.sqrt(G[i][j]["SD"]**2+ G[i][j]["time"]**2)
		mu=math.log(G[i][j]["time"]**2/phi)
		sigma=math.sqrt(math.log(phi**2/(G[i][j]["time"]**2)))
		
		data = np.random.lognormal(mu,sigma, 100)		
		data=list(map(lambda x:max(0.0,x),data-G[i][j]["fftt"]))


		#Fits this data to one of the given distributions
		params = fit_to_all_distributions(data,d_names)
		# print()
		best_dist_chi, best_chi, params_chi, dist_results_chi = get_best_distribution_using_chisquared_test(data, params,d_names)

		#Saves the info in the new graph
		NG.add_edge(i,j, Cost=G[i][j]["cost"], tmin=G[i][j]["fftt"],distr=best_dist_chi,params=params_chi)

		#text=str(i)+";"+str(j)+";"+str(int(G[i][j]["cost"]))+";"+str(G[i][j]["fftt"])+";"+str([list(map(lambda x: round(float(x),10),k))[0]for k in G[i][j]["alpha"]])+";"+str([list(map(lambda x: round(float(x),10),k))for k in G[i][j]["T"]])+'\n'
		
	nx.write_gpickle(NG,file+'.gpickle')
	return NG

def rand_inst_to_txt(G,file,d_names):
	os.chdir(r'C:\Users\David Corredor M\Desktop\Thesis\2. Methodology\2.2 S-Pulse Algorithm\SPulseD\networks')
	cont=0
	for i,j in G.edges():
		#Creates data with the given moments
		phi=math.sqrt(G[i][j]["SD"]**2+ G[i][j]["time"]**2)
		mu=math.log(G[i][j]["time"]**2/phi)
		sigma=math.sqrt(math.log(phi**2/(G[i][j]["time"]**2)))
		
		data = np.random.lognormal(mu,sigma, 100)		
		data=list(map(lambda x:max(0.0,x),data-G[i][j]["fftt"]))
		#Fits this data to one of the given distributions
		params = fit_to_all_distributions(data,d_names)
		# print()
		best_dist_chi, best_chi, params_chi, dist_results_chi = get_best_distribution_using_chisquared_test(data, params,d_names)

		#Saves the info in the new graph
		NG.add_edge(i,j, Cost=G[i][j]["cost"], tmin=G[i][j]["fftt"],distr=best_dist_chi,params=params_chi)

		#text=str(i)+";"+str(j)+";"+str(int(G[i][j]["cost"]))+";"+str(G[i][j]["fftt"])+";"+str([list(map(lambda x: round(float(x),10),k))[0]for k in G[i][j]["alpha"]])+";"+str([list(map(lambda x: round(float(x),10),k))for k in G[i][j]["T"]])+'\n'
		
		
		
	
	nx.write_gpickle(NG,file+'.gpickle')

def fit_to_all_distributions(data,dist_names):    
    params = {}
    for dist_name in dist_names:
        try:
            dist = getattr(st, dist_name)
            param = dist.fit(data)

            params[dist_name] = param
        except Exception:
            print("Error occurred in fitting")
            params[dist_name] = "Error"

    return params 

def get_best_distribution_using_chisquared_test(data,params,dist_names):

    histo, bin_edges = np.histogram(data, bins='auto', normed=False)
    number_of_bins = len(bin_edges) - 1
    observed_values = histo
    
    dist_results = []

    for dist_name in dist_names:
        param = params[dist_name]
        if (param != "Error"):
            # Applying the SSE test
            arg = param[:-2]
            loc = param[-2]
            scale = param[-1]
            # print(getattr(st, dist_name))
            cdf = getattr(st, dist_name).cdf(bin_edges, loc=loc, scale=scale, *arg)
            expected_values = len(data) * np.diff(cdf)
            c , p = st.chisquare(observed_values, expected_values, ddof=number_of_bins-len(param))
            dist_results.append([dist_name, c, p])
    
    # select the best fitted distribution
    best_dist, best_c, best_p = None, sys.maxsize, 0

    for item in dist_results:
        name = item[0]
        c = item[1]
        p = item[2]
        if (not math.isnan(c)):
            if (c < best_c):
                best_c = c
                best_dist = name
                best_p = p

    # print the name of the best fit and its p value

    # print("Best fitting distribution: " + str(best_dist))
    # print("Best c value: " + str(best_c))
    # print("Best p value: " + str(best_p))
    # print("Parameters for the best fit: " + str(params[best_dist]))

    return best_dist, best_c, params[best_dist], dist_results

def loadNet(net='Chicago-Sketch'):
	'''
	Loads the {net}_net.txt file. Looks for the {net} folder in '../../data/Networks/'
	Args:
		net(String): Folder name
	Return:
		arcs (dict): keys:=(tail, head), values:=(Cost,E[Time],Fft,Sigma)
	'''		
	f=open(f'../../data/Networks/{net}/{net}_net.txt')
	arcs=dict()
	f.readline()
	for l in f:
		li=list(map(float,l.replace('\n','').split('\t')))
		arcs[int(li[0]),int(li[1])]=tuple(li[2:])
	return arcs

def createSceanrios(nScen,n,net='Chicago-Sketch'):
	'''
	Creates for the given network nScen scenarios with the following properties:
		-
		-
		-
	Args:
		nScen(int): Must be odd number in order to create (n-1)/2 scenarios under the expected value and (n-1)/2 above.
		n(int): Number of observations per scenario. Nothe that the agregated dataset has lenght of nScen*n for each arc.
		net(String): Folder name

	Returns:
		None: creates for each scenario a file with the data, and a file for the aggregated data. These files are saved in '../../data/Networks/{net}/Scenarios/'.

	'''
	files=[open(f'../../data/Networks/{net}/Scenarios/scen{k+1}.txt','w') for k in range(nScen)]+[open(f'../../data/Networks/{net}/Scenarios/scenTotal.txt','w')]
	arcs=loadNet(net=net)
	for i,j in arcs.keys():
		c,t,fft,sigma=arcs[i,j]
		delta=(t-fft)/(1.5*(nScen-1)/2)
		TData=[]
		for k in range(int(-(nScen-1)/2),int((nScen-1)/2+1)):
			data=createData(mu=t+k*delta,fft=fft,n=n)
			files[k+int((nScen-1)/2)].write(f'{i}\t{j}\t{data}\n')
			TData+=data
		files[-1].write(f'{i}\t{j}\t{TData}\n')		
	for f in files: f.close()

def createDataIndependent(n,net='Chicago-Sketch'):
	'''
	Creates the data files for the given network
	Args:
		n(int): Number of observations per . Nothe that the agregated dataset has lenght of nScen*n for each arc.
		net(String): Folder name

	Returns:
		None: creates for each scenario a file with the data, and a file for the aggregated data. These files are saved in '../../data/Networks/{net}/Scenarios/'.

	'''
	file=open(f'../../data/Networks/{net}/Independent/data_cv{CV}.txt','w') 
	arcs=loadNet(net=net)
	for i,j in arcs.keys():
		c,t,fft,sigma=arcs[i,j]
		
		data=createData(mu=t,fft=fft,n=n)
		file.write(f'{i}\t{j}\t{data}\n')
		
	file.close()

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

	return data

def fitPHScenarios(nScen,nPhases,net='Chicago-Sketch'):
	'''
	fits a PH random variable to each scenario created
	Args:
		nScen(int): Must be odd number in order to create (n-1)/2 scenarios under the expected value and (n-1)/2 above.		
		net(String): Folder name

	Returns:
		None: creates for each scenario a file with the data, and a file for the aggregated data. These files are saved in '../../data/Networks/{net}/Scenarios/'. 
	'''

	files=[open(f'../../data/Networks/{net}/Scenarios/scen{k+1}.txt','r') for k in range(nScen)]+[open(f'../../data/Networks/{net}/Scenarios/scenTotal.txt','r')]	
	filesOut=[open(f'../../data/Networks/{net}/Scenarios/PHFit{nPhases}_scen{k+1}.txt','w') for k in range(nScen)]+[open(f'../../data/Networks/{net}/Scenarios/PHFit{nPhases}_scenTotal.txt','w')]	
	arcs=loadNet(net=net)
	N=len(arcs.keys())
	n=0
	for i,j in arcs.keys():
		n+=1
		for fin,fout in zip(files,filesOut):
			l=fin.readline()			
			ii,jj,data=l.replace('\n','').split('\t')
			data=list(map(float,data.replace('[','').replace(']','').split(', ')))
			
			ph,loglike=get_PH_HE(data,nPhases)			
			fout.write(f'{i}\t{j}\t{arcs[i,j][0]}\t{arcs[i,j][2]}\t{[list(map(lambda x: round(float(x),10),k))[0]for k in ph.alpha.tolist()]}\t{[list(map(lambda x: round(float(x),10),k))for k in ph.T.tolist()]}\n')
		print(n/N)
	for fin,fout in zip(files,filesOut): (fin.close(),fout.close())

def fitIndependent(nPhases,net='Chicago-Sketch'):
	'''
	Fits a PH distribution and arbitrari distribution to each arc in the network.

	Args:	
		nPhases(list): list of number of phases to fit
		net(String): Folder name

	Returns:
		None: Creates a file with the aggregated data. These files are saved in '../../data/Networks/{net}/Independent/'. 
	'''
	file=open(f'../../data/Networks/{net}/Independent/data_cv{CV}.txt','r')
	fileOut=[open(f'../../data/Networks/{net}/Independent/PHFit{np}_cv{CV}.txt','w')for np in nPhases]+[open(f'../../data/Networks/{net}/Independent/DistrFit_cv{CV}.txt','w')]
	arcs=loadNet(net=net)
	N=len(arcs.keys())
	n=0

	fittingTimes={np:0 for np in nPhases} #Total fitting times
	fittingTimes['distr']=0

	for i,j in arcs.keys():
		n+=1		
		l=file.readline()			
		ii,jj,data=l.replace('\n','').split('\t')
		data=list(map(float,data.replace('[','').replace(']','').split(', ')))
		#fit PHtypes
		for ii,np in enumerate(nPhases):
			ph,loglike=get_PH_HE(data,np)
			tFit=time.time()			
			fileOut[ii].write(f'{i}\t{j}\t{arcs[i,j][0]}\t{arcs[i,j][2]}\t{[list(map(lambda x: round(float(x),10),k))[0]for k in ph.alpha.tolist()]}\t{[list(map(lambda x: round(float(x),10),k))for k in ph.T.tolist()]}\n')
			fittingTimes[np]+=time.time()-tFit
		
		#fit Distributions		
		tFit=time.time()			
		params = fit_to_all_distributions(data,d_names)		
		best_dist_chi, best_chi, params_chi, dist_results_chi = get_best_distribution_using_chisquared_test(data, params,d_names)
		fittingTimes['distr']+=time.time()-tFit
		fileOut[-1].write(f'{i}\t{j}\t{best_dist_chi}\t{params_chi}\n')		
		print(n/N)

	file.close()
	for f in fileOut: f.close()

	fTimes=open(f'../../data/Networks/Fitting_Times.txt','a')
	for np in nPhases:
		fTimes.write(f'PhaseType fit {np}\t{net}\t{fittingTimes[np]}\n')

	fTimes.write(f'Distrubution\t{net}\t{fittingTimes["distr"]}\n')
	fTimes.close()


if __name__ == '__main__':
	########################################################
	'''
	Setting of some global parameters
	'''
	CV=0.5 #CV (float):	Coefficient of variation
	d_names = ['lognorm','gamma','weibull_min']
	########################################################

	#createDataIndependent(n=500,net='Chicago-Sketch')
	fitIndependent(nPhases=[3,5],net='Chicago-Sketch')














