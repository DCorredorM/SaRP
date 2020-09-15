#http://www.bgu.ac.il/~bargera/tntp/
import sys
from os import path
sys.path.append(path.abspath(r'C:\Users\d.corredor\OneDrive - Universidad de los andes\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse'))
#sys.path.insert(1, r'C:\Users\David Corredor M\OneDrive - Universidad de Los Andes\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse')
#sys.path.insert(1, r'/Users/davidcorredor/OneDrive - Universidad de los Andes/Thesis/2. Methodology/2.2 S-Pulse Algorithm/S-PULSE implementation/Code_S-pulse')
from PhaseType import *
from Fitter import *
import networkx as nx
import os
import matplotlib.pyplot as plt
import math
import numpy as np
import s_pulse_MC as MC
import s_pulse  as PH_pulse
from PhaseType import *
from Fitter import *
from pulse import *
import scipy.stats as st
import networkx as nx
import time
import numpy as np
from mpmath import *
import random

##############################################################################
#Creates the basic net and writes a .txt (city_net.txt)
#Tail	Head	Cost	E[Time]	Fft	Sigma(Time)
#(i,j)->(Capacity=cap, Lenght=leng, FreeflowTime=fftt, B=B,Power=pow, Toll=toll,Flow*=flow, Cost=Cost, E[Time]=Time,Sigma(Time)=SD)
def create_net(city,write=False):
	os.chdir(r'TransportationNetworks-master/'+city)
	
	G=nx.DiGraph()
	toll_factor=100
	distance_factor=1.69*100
	#Fills the position of the nodes from the file ended by node.tntp
	nodes= [i for i in os.listdir() if "node" in i][0]
	f=open(nodes)
	pos={}
	for x in f:
		x=x.split("\t")
		if x[0]!='node':
			#print(x)
			x=list(map(float,x))
			G.add_node(int(x[0]))
			pos[x[0]]=(x[1],x[2])
	f.close()	
	net= [i for i in os.listdir() if "net" in i][0]
	f=open(net)
	for x in f:
		x=x[1:-3].split("\t")
		#print(x)
		t=True
		try:
			type(int(x[0]))
		except Exception as e:
			t=False
		if t:
			G.add_edge(int(x[0]),int(x[1]),cap=int(float(x[2])),leng=float(x[3]),fftt=float(x[4])*60,B=float(x[5]),pow=float(x[6]),toll=float(x[8]))
			G[int(x[0])][int(x[1])]["Cost"]=int(toll_factor*G[int(x[0])][int(x[1])]["toll"]+distance_factor*G[int(x[0])][int(x[1])]["leng"])
			if G[int(x[0])][int(x[1])]["B"]==0:
				G[int(x[0])][int(x[1])]["B"]=0.15
	f.close()

	#If there is information of the flow in the net uses that information to estimate the expected value of the time of each arc, if not creates a random instance of the time.
	try:
		flow= [i for i in os.listdir() if "flow" in i][0]
		f=open(flow)
		for x in f:
			x=x.split('\t')
			t=True
			try:
				type(int(x[0]))
			except Exception as e:
				t=False
			if t:				
				G[int(x[0])][int(x[1])]["flow"]=int(float(x[2][:-1]))	
		f.close()
		edg=list(G.edges()).copy()
		for i,j in edg:
			if G[i][j]["fftt"]==0:
					G[i][j]["fftt"]=60*(60*G[i][j]["leng"]/50)#(mph))
			#Link travel time = free flow time * ( 1 + B * (flow/capacity)^Power ).**2
			G[i][j]["Time"]=G[i][j]["fftt"]*(1+G[i][j]["B"]*(G[i][j]["flow"]/G[i][j]["cap"])**G[i][j]["pow"])+5
			G[i][j]["SD"]=max(1/5*(G[i][j]["Time"]-G[i][j]["fftt"])*(G[i][j]["flow"]/G[i][j]["cap"]),5)
			#G[i][j]["Cost"]=int(100*G[i][j]["Cost"]*(G[i][j]["fftt"]/(G[i][j]["time"]**2)) )
	except Exception as e:
		edg=list(G.edges()).copy()
		for i,j in edg:
			if G[i][j]["fftt"]==0:
					G[i][j]["fftt"]=60*(60*G[i][j]["leng"]/50)#(mph))
			#Link travel time = free flow time * ( 1 + B * (flow/capacity)^Power ).**2
			rand_FC_ratio=random.random()*0.5+0.5
			G[i][j]["Time"]=G[i][j]["fftt"]*(1+G[i][j]["B"]*(rand_FC_ratio)**G[i][j]["pow"])
			G[i][j]["SD"]=max(1/2*(G[i][j]["Time"]-G[i][j]["fftt"])*(rand_FC_ratio),5)
			
			#G[i][j]["Cost"]=int(100*G[i][j]["Cost"]*(G[i][j]["fftt"]/(G[i][j]["time"]**2)) )
		#raise e

	#nx.draw(G,pos,node_size=1,node_color="black",arrows=False,width=0.1)
	#plt.show()
	#Takes me back to the original directory
	

	if write:
		os.chdir('..')
		os.chdir('..')
		os.chdir(r'Networks/'+city)
		#Writes tis graph in a txt in the eachs net folder
		file1 = open(city+'_net'+'.txt',"w")
		file1.write('Tail\tHead\tCost\tE[Time]\tFft\tSigma(Time)\n')
		for i,j in G.edges():
			text=str(int(i))+"\t"+str(int(j))+"\t"+str(int(G[i][j]["Cost"]))+"\t"+str(int(G[i][j]["Time"]))+"\t"+str(G[i][j]["fftt"])+"\t"+str(int(G[i][j]["SD"]))+'\n'
			file1.write(text)
		file1.close()

		file1 = open(city+'_pos'+'.txt',"w")
		file1.write('Node\tX\tY\n')
		for i,(x,y) in pos.items():
			text=str(int(i))+"\t"+ str(x)+"\t"+str(y)+"\t"+'\n'
			print(text)
			file1.write(text)
		file1.close()

	os.chdir('..')
	os.chdir('..')
	return G, pos
def plot_net(city):
	G,pos=create_net(city)
	nx.draw_networkx_edges(G, pos=pos ,arrows=True)#,width=0.5
	nodes=nx.draw_networkx_nodes(G, pos=pos ,node_color='white')#,node_size=0)
	nodes.set_edgecolor('black')
	nx.draw_networkx_labels(G,pos=pos,font_size =8)#,font_color='white'
	plt.show()


##############################################################################
#Simulates times for each arc and fits distributions to each arc.
def fit_to_all_distributions(data,dist_names):    
    params = {}
    for dist_name in dist_names:
        try:
            dist = getattr(st, dist_name)
            param = dist.fit(data)

            params[dist_name] = param
        except Exception as e:
            print("Error occurred in fitting ",dist_name )
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
            #print(getattr(st, dist_name))
            

            cdf = getattr(st, dist_name).cdf(bin_edges, loc=loc, scale=scale, *arg)
            expected_values = len(data) * np.diff(cdf)
            Chi=sum((observed_values[i]-expected_values[i])**2/expected_values[i] for i in range(number_of_bins))
            p_val=1-st.chi2.cdf(Chi, number_of_bins-len(param))
            

            '''
            plt.bar(range(number_of_bins),expected_values)
            plt.bar(range(number_of_bins),observed_values)
            plt.title(dist_name)
            plt.legend(["Expected","Observed"])
            print("mi chi: ",Chi)
            print("su chi: ",c)
            
           	print("Chi critico: ", st.chi2.ppf(.95, number_of_bins-len(param)))
           	print("mi p-val: ", p_val)
           	print("Su p-val: ", p)
            plt.show()
            '''
            
            dist_results.append([dist_name, Chi, p_val])
    
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
                #print("best p: ",p)

    # print the name of the best fit and its p value
    '''
    print("Best fitting distribution: " + str(best_dist))
    print("Best c value: " + str(best_c))
    print("Best p value: " + str(best_p))
    print("Parameters for the best fit: " + str(params[best_dist]))
	'''

    return best_dist, best_c, params[best_dist], dist_results, best_p
def fit_and_print(G,city,Continue=False,plot=False):
	os.chdir(r'Networks/'+city)
	cont=0
	d_names = ['lognorm','gamma','weibull_min']#'weibull_min','burr']
	n_good_fits=0
	p_value_prom=0

	if Continue:
		for i,j in G.edges():
			G[i][j]["fitted"]=False
			#print(i,j)
		fited_distrs = open(city+'_fitts'+'.txt')
		cont = 0
		for x in fited_distrs:
			cont+=1
			
			x=x.split('\t')
			
			if x[0]!='Tail':
				print(int(float(x[0])),int(float(x[1])))
				G[int(float(x[0]))][int(float(x[1]))]["fitted"]=True
		
		fited_distrs.close()
		fited_distrs = open(city+'_fitts'+'.txt',"a")
		
		data_file=open(city+'_data'+'.txt',"a")


		t=time.time()
		tot=len(G.edges())
		s_time=1.5*cont
		print('Empiezo con ',cont/tot*100,'% ya fitteado y estimo ',(tot-cont)*1.5/60,'mins mas.')
		for i,j in G.edges():
			if not G[i][j]["fitted"]:	
				#Creates data with the given moments
				tt=time.time()
				mean=abs(G[i][j]["Time"]-G[i][j]["fftt"])+0.00001
			
				phi=math.sqrt(G[i][j]["SD"]**2+ mean**2)
				mu=math.log(mean**2/phi)
				sigma=math.sqrt(math.log(phi**2/(mean**2)))
				Good_fit=False
			while not Good_fit:
				data = np.random.lognormal(mu,sigma, 1000)	
				#data=list(map(lambda x:max(0.0,x),data-G[i][j]["fftt"]))
				mu_data=sum(data)/len(data)
				

				#Fits this data to one of the given distributions

				params = fit_to_all_distributions(data,d_names)
				best_dist_chi, best_chi, params_chi, dist_results_chi,dist_p_val = get_best_distribution_using_chisquared_test(data, params,d_names)

				
				dat=getattr(sts,best_dist_chi).rvs(size=1000,*params_chi)
				
				
				gen_mean=np.mean(dat) 
				gen_sd=math.sqrt(np.var(dat) )
				if abs(mu_data-gen_mean)/mu_data>1 and abs(mu_data-gen_mean)>100:	
					print("###########################################")
					print("distr ",best_dist_chi)
					print("Real mean ",mean)
					print("Data mean ",mu_data)
					print("Gener Mean ",np.mean(dat) )
					
					print("Rela sd ",G[i][j]["SD"])
					print("Data sd ",math.sqrt(np.var(data)))
					print("Gener sd ",math.sqrt(np.var(dat) ))

					print("p-val: ",dist_p_val)
					print("###########################################")
					xx = np.linspace(0, max(data), 100)
					
					arg = params_chi[:-2]
					c=params_chi[-3]
					loc = params_chi[-2]
					scale = params_chi[-1]


				else:
					Good_fit=True

			#Saves the info in the new graph
			#NG.add_edge(i,j, Cost=G[i][j]["Cost"], tmin=G[i][j]["fftt"],distr=best_dist_chi,params=params_chi)
			#text=str(i)+"\t"+str(j)+"\t"+best_dist_chi+"\t"+str([list(map(lambda x: round(float(x),10),k))[0]for k in G[i][j]["alpha"]])+"\t"+str([list(map(lambda x: round(float(x),10),k))for k in G[i][j]["T"]])+'\n'
			text1=str(int(i))+"\t"+str(int(j))+"\t"+best_dist_chi+"\t"+str(params_chi)+'\n'
			fited_distrs.write(text1)
			text2=str(int(i))+"\t"+str(int(j))+"\t"+str(list(data))+'\n'
			data_file.write(text2)
			
			s_time+=time.time()-tt
			if cont%100==0:
				print('voy ',cont, 'de ', tot,': ',cont/tot*100,'%')
				print('tiempo estimado: ',tot*(s_time/(cont+1))*(1-cont/tot)/60," mins")
				print("voy con tiempo de ",s_time)
				print("data mean ", sum(data)/(len(data)))
				print("real mean ",G[i][j]["Time"]-G[i][j]["fftt"])
				print("p-val: ",dist_p_val)

			cont+=1
		t_fit=time.time()-t
	else:
		fited_distrs = open(city+'_fitts'+'.txt',"w")
		fited_distrs.write('Tail\tHead\tDistribution\tParameters\n')
		data_file=open(city+'_data'+'.txt',"w")

		t=time.time()
		tot=len(G.edges())
		s_time=0
		for i,j in G.edges():
			#Creates data with the given moments
			tt=time.time()
			mean=abs(G[i][j]["Time"]-G[i][j]["fftt"])+0.00001
			
			phi=math.sqrt(G[i][j]["SD"]**2+ mean**2)
			mu=math.log(mean**2/phi)
			sigma=math.sqrt(math.log(phi**2/(mean**2)))
			
			#print("sd ",G[i][j]["SD"])
			#print("mean " ,mean)
			#print("s/e ",G[i][j]["SD"]**2/(mean**2))
			#mu=math.log(mean)-1/2*math.log(G[i][j]["SD"]**2/(mean**2) ) 
			#sigma=math.sqrt(math.log(G[i][j]["SD"]**2/(mean**2) ))
			Good_fit=False
			while not Good_fit:
				data = np.random.lognormal(mu,sigma, 1000)	
				#data=list(map(lambda x:max(0.0,x),data-G[i][j]["fftt"]))
				mu_data=sum(data)/len(data)
				

				#Fits this data to one of the given distributions

				params = fit_to_all_distributions(data,d_names)
				best_dist_chi, best_chi, params_chi, dist_results_chi,dist_p_val= get_best_distribution_using_chisquared_test(data, params,d_names)

				
				dat=getattr(sts,best_dist_chi).rvs(size=1000,*params_chi)

				if plot:
					if best_dist_chi=="lognorm":
						xx = np.linspace(0, max(data)*1.2, 100)
						s,loc,scale=params_chi[0:3]
						pdf=lambda x: getattr(sts,best_dist_chi).pdf(x=x, s=s, loc=loc, scale=scale)
						y2 = list(map(pdf,xx))
						plt.hist(data,bins=100, density=True)
						plt.plot(xx,y2, color="lime")
						plt.legend(['PH-fit'])
						plt.title("fitted "+best_dist_chi+" to "+str((i,j)))
						
						plt.savefig("./Distr_fit_plots/fit"+"_"+str(i)+"_"+str(j)+"_.png")
						plt.clf()

						#plt.show()
						#print("Arc: ",str((i,j)))
						#print("Distr: ",best_dist_chi)
						#print("Chi: ",best_chi)
						#print("P-val: ",dist_p_val)


					elif best_dist_chi=="gamma":
						xx = np.linspace(0, max(data)*1.2, 100)
						a,loc,scale=params_chi[0:3]
						pdf=lambda x: getattr(sts,best_dist_chi).pdf(x=x, a=a, loc=loc, scale=scale)
						y2 = list(map(pdf,xx))
						plt.hist(data,bins=100, density=True)
						plt.plot(xx,y2, color="lime")
						plt.legend(['PH-fit'])
						plt.title("fitted "+best_dist_chi+" to "+str((i,j)))
						plt.title("fitted "+best_dist_chi+" to "+str((i,j)))
						
						plt.savefig("./Distr_fit_plots/fit"+"_"+str(i)+"_"+str(j)+"_.png")
						plt.clf()
					elif best_dist_chi=="expon":
						xx = np.linspace(0, max(data)*1.2, 100)
						loc,scale=params_chi[0:2]
						pdf=lambda x: getattr(sts,best_dist_chi).pdf(x=x, loc=loc, scale=scale)
						y2 = list(map(pdf,xx))
						plt.hist(data,bins=100, density=True)
						plt.plot(xx,y2, color="lime")
						plt.legend(['PH-fit'])
						plt.title("fitted "+best_dist_chi+" to "+str((i,j)))
						plt.title("fitted "+best_dist_chi+" to "+str((i,j)))
						
						plt.savefig("./Distr_fit_plots/fit"+"_"+str(i)+"_"+str(j)+"_.png")
						plt.clf()

						#plt.show()
						print("Arc: ",str((i,j)))
						print("Distr: ",best_dist_chi)
						print("Chi: ",best_chi)
						print("P-val: ",dist_p_val)
						

					

				#plt.hist(data,bins=100)
				#plt.show()
				gen_mean=np.mean(dat) 
				gen_sd=math.sqrt(np.var(dat) )
				#if abs(mu_data-gen_mean)/mu_data>1 and abs(mu_data-gen_mean)>100:	
				if dist_p_val>=0.05:
					Good_fit=True
					n_good_fits+=1
					p_value_prom+=dist_p_val
				else:
					pass
					#print("Repeat")
					'''
					print("###########################################")
					print("distr ",best_dist_chi)
					print("Real mean ",mean)
					print("Data mean ",mu_data)
					print("Gener Mean ",np.mean(dat) )
					
					print("Rela sd ",G[i][j]["SD"])
					print("Data sd ",math.sqrt(np.var(data)))
					print("Gener sd ",math.sqrt(np.var(dat) ))
					print("P-val: ",dist_p_val)
					print("###########################################")
					'''
					
					

			
			

			#Saves the info in the new graph
			#NG.add_edge(i,j, Cost=G[i][j]["Cost"], tmin=G[i][j]["fftt"],distr=best_dist_chi,params=params_chi)
			#text=str(i)+"\t"+str(j)+"\t"+best_dist_chi+"\t"+str([list(map(lambda x: round(float(x),10),k))[0]for k in G[i][j]["alpha"]])+"\t"+str([list(map(lambda x: round(float(x),10),k))for k in G[i][j]["T"]])+'\n'
			
			text1=str(int(i))+"\t"+str(int(j))+"\t"+best_dist_chi+"\t"+str(params_chi)+'\n'
			fited_distrs.write(text1)
			text2=str(int(i))+"\t"+str(int(j))+"\t"+str(list(data))+'\n'
			data_file.write(text2)
			
			s_time+=time.time()-tt
			if cont%100==0:
				print('voy ',cont, 'de ', tot,': ',cont/tot*100,'%')
				print('tiempo estimado: ',tot*(s_time/(cont+1))*(1-cont/tot)/60," mins")
				print("voy con tiempo de ",s_time)
				print("data mean ", sum(data)/(len(data)))
				print("real mean ",G[i][j]["Time"]-G[i][j]["fftt"])

			cont+=1
		t_fit=time.time()-t

		#print("La prop de fits que pasan son: ",n_good_fits/len(G.edges))
		#print("El p-value promedio en la red es de  ",p_value_prom/len(G.edges))
	fited_distrs.close()
	data_file.close()
	os.chdir('..')
	fit_times = open('Fitting_Times'+'.txt',"a")
	fit_times.write('Distribution\t'+city+'\t'+str(t_fit)+'\t'+str(n_good_fits/len(G.edges)) + '\t' + str(p_value_prom/len(G.edges) )+'\n')
	fit_times.close()
	os.chdir('..')
	

##############################################################################
#Fits the city travel times to PH distributions and saves a txt file city_PHfit_nPhases.txt with the parameters of each PH
def fit_PH_and_print(G,city,nPhases,plot=False):
	os.chdir(r'Networks/'+city)
	data_file=open(city+'_data'+'.txt')
	file1 = open(city+'_PHfit_'+str(nPhases)+'.txt',"w") 
	cont=0
	tot=len(G.edges)
	s_time=0
	t=time.time()

	n_good_fits=0
	p_value_prom=0
	for x in data_file:
		tt=time.time()
		
		x=x.split('\t')
		#print(x)
		if len(x)>2:
			i,j=int(float(x[0])),int(x[1])
			s_Data=x[2].replace('[','')
			s_Data=s_Data.replace(']','')
			s_Data=s_Data.replace('\n','')
			
			lis=filter(lambda x: x!='',  s_Data.split(', '))
			#print(s_Data)
			data=list(map(float,lis))

			#print("len", len(data))
			
			
			phi=math.sqrt(G[i][j]["SD"]**2+ G[i][j]["Time"]**2)
			mu=math.log(G[i][j]["Time"]**2/phi)
			sigma=math.sqrt(math.log(phi**2/(G[i][j]["Time"]**2)))

			ph,log_like=get_PH_HE(data,nPhases)

			#################
			#Calcs GOF test for PH fit.
			histo, bin_edges = np.histogram(data, bins='auto', normed=False)
			number_of_bins = len(bin_edges) - 1
			observed_values = histo    

			nnn=len(data)
			#cdf = getattr(st, dist_name).cdf(bin_edges, loc=loc, scale=scale, *arg)
			cdf=np.array([float(ph.cdf(i)) for i in bin_edges])
			expected_values = len(data) * np.diff(cdf)
			
			#expected_values =  * np.diff(cdf)
			
			Chi=sum((observed_values[i]-expected_values[i])**2/expected_values[i] for i in range(number_of_bins))

			p_val=1-st.chi2.cdf(Chi, number_of_bins-2)
			#print("i,j ",(i,j))
			#print("Chi: ",Chi)
			#print("p_val: ",p_val)

			if p_val>0.05:
				n_good_fits+=1
			p_value_prom+=p_val
			
			if plot:
				
				plt.hist(data,bins=100,density=True)
				ph.plot_pdf(max=max(data)*1.4)
				plt.title("Phase-type Fit for "+str((i,j)))
				#plt.show()
				plt.savefig("./PH_fit_plots/PH_"+str(nPhases)+"/fit"+"_"+str(i)+"_"+str(j)+"_.png")
				plt.clf()
				


			G[i][j]["alpha"]=ph.alpha.tolist()
			G[i][j]["T"]=ph.T.tolist()

			text=str(i)+"\t"+str(j)+"\t"+str(int(G[i][j]["Cost"]))+"\t"+str(G[i][j]["fftt"])+"\t"+str([list(map(lambda x: round(float(x),10),k))[0]for k in G[i][j]["alpha"]])+"\t"+str([list(map(lambda x: round(float(x),10),k))for k in G[i][j]["T"]])+'\t'+str(log_like)+'\n'
			#print(text)
			file1.write(text)
			s_time+=time.time()-tt
			
			if cont%100==0:
				print('voy ',cont, 'de ', tot,': ',cont/tot*100,'%')
				print('tiempo estimado: ',tot*(s_time/(cont+1))*(1-cont/tot)/60," mins")
				print("voy con tiempo de ",s_time)
				print("data mean ", sum(data)/(len(data)))
				print("real mean ",G[i][j]["Time"]-G[i][j]["fftt"])
			cont+=1

	fit_time=time.time()-t
	os.chdir('..')
	fit_times = open('Fitting_Times'+'.txt',"a")
	fit_times.write('PhaseType fit '+str(nPhases)+ '\t'+city+'\t'+str(fit_time)+'\t'+str(n_good_fits/len(G.edges)) + '\t' + str(p_value_prom/len(G.edges) )+'\n')
	
	fit_times.close()
	os.chdir('..')
	file1.close()


##############################################################################
#Creates instances (source, targer, T_max) for the specified city
#Makes clusters based on the location
from operator import itemgetter

def split_x(pos):
	
	x=0	
	n=len(pos.values())
	print(pos)
	for i,j in pos.values():
		x+=i
	
	try:
		x=x/n
	except:
		x=0

	c1=[]
	c2=[]
	for n,(i,j) in pos.items():
		if i>=x:
			c1.append(n)
		else:
			c2.append(n)
	return c1,c2
def split_y(pos):
	y=0	
	n=len(pos.values())
	for i,j in pos.values():
		y+=j
	y=y/n

	c1=[]
	c2=[]
	for n,(i,j) in pos.items():
		if j>=y:
			c1.append(n)
		else:
			c2.append(n)
	return c1,c2

def create_random_instances(city,G,n_splits,pos,n_inst,tightness,inf_tightness,search_time,w_a,plot=False):
	
	part=list(split_x(pos))
	cont=1
	while cont<n_splits:
		
		parts=part.copy()
		part=[]
		for p in parts:
			post={k:pos[k] for k in p}
			if cont%2:
				part+=list(split_y(post))

			else:
				part+=list(split_x(post))
		cont+=1		
	
	part=list(filter(lambda x: x!=[],part))

	if plot:
		#nx.draw(G, pos=pos ,arrows=False,width=0.1,node_size=0.001)
		#plt.show()
		
		number_of_colors = 2**n_splits
		
		color = ["#"+''.join([random.choice('0123456789ABCDEF') for jj in range(6)])for ii in range(number_of_colors)]
		#color =['black']#'r','g','orange','b']
		nx.draw_networkx_edges(G, pos=pos ,arrows=False,width=0.1)
		
		for i,j in enumerate(part):		
			nx.draw_networkx_nodes(G, pos=pos, nodelist=j, node_size=1 ,node_color=color[i])
		plt.show()
		
	
	#d()
	os.chdir(r'Networks/'+city)
	file1 = open(city+'_Instances'+'.txt',w_a)
	#file1.write('Source\tTarget\tT_max\tTightness\n')
	num_it=0
	done=False
	while not done:
		zones=list(range(len(part)))
		s_zone=random.choice(zones)

		zones.pop(s_zone)
		t_zone=random.choice(zones)
		s=random.choice(part[s_zone])
		t=random.choice(part[t_zone])

		#s=2
		#t=13

		tightness=tightness

		
		PG=pulse_graph(G=G,T_max=0,source=s, target=t,tightness=1,pos=pos)
		pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()		
		
		T_max=PG.T_max
		T_max_star=T_max
		#if plot: PG.draw_graph()
		##############################################################################
		#To do: While the relaiability of the path from the CSP is greater than tightness, reduce the time constraint and find again this path.
		stop=False
		delta=0.8
		try:
			S_pat=nx.shortest_path(G,s,t,'Time')
		except:
			stop=True

		

		time_creit=time.time()
		pri=False
		if len(S_pat)<5:
			stop=True
			print('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n',S_pat)
		while not stop:

			
			os.chdir('..')
			os.chdir('..')
			
			rho=eval_path(pulse_path,city,1000,T_max)
			sp_rho=eval_path(S_pat,city,1000,T_max)
			print("se tiene que , ",rho,sp_rho,T_max,S_pat)
			if rho!=0:
				T_max_star=T_max

			os.chdir(r'Networks/'+city)
			
			if rho<tightness and sp_rho>=inf_tightness:
				stop=True
				pri=True
			elif sp_rho>=inf_tightness:
				T_max=T_max*delta
				PG=pulse_graph(G=G,T_max=T_max,source=s, target=t,pos=pos)
				pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()		
				T_max=PG.T_max
			
			elif sp_rho<=inf_tightness:
				T_max=T_max*(1/0.9)
				delta=delta +0.5*(1-delta)
				#T_max=T_max*delta
				PG=pulse_graph(G=G,T_max=T_max,source=s, target=t,pos=pos)
				pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()		
				T_max=PG.T_max
			
			if time.time()-time_creit>=search_time:
				stop=True
		
		
		if num_it==n_inst:
			done=True
			
		if pri:	
			num_it+=1
			print("#############")
			print("El tiempo maximo es de ",T_max_star)
			print("El min time de la red es de ",PG.G.nodes[s]["s_time"])
			print("Estamos a  ",(T_max_star-PG.G.nodes[s]["s_time"]))
			
			print("#############")
			text=str(s)+'\t'+str(t)+'\t'+str(T_max_star)+'\t'+str(tightness)+'\n'
			print(text)
			file1.write(text)
	file1.close()
	os.chdir('..')
	os.chdir('..')


##############################################################################
#Run instances for the given city

def creates_graphs(city):
	DG=nx.DiGraph()
	#os.chdir(r'Networks/'+city)

	##########################
	#Creates deterministic graph
	file1 = open(city+'_net.txt','r')
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
	file1 = open(city+'_fitts.txt','r')
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

	##########################
	#Creates graph with PH 
	SG_PH=DG.copy()
	file1 = open(city+'_PHfit_3.txt','r')
	for x in file1:
		x=x[:-1].split('\t')
		#print(x)
		if x[0]!='Tail':		
			def s_to_list(s):
				s=s.replace('[','')
				s=s.replace(']','')
				return(list(map(float,s.split(', '))))			
			t=s_to_list(x[4])			
			T=x[5].split('], [')		
			T=list(map(s_to_list,T))
			

			i,j=int(float(x[0])),int(float(x[1]))
			SG_PH[i][j]["tRV"]=PH(T,t)
			
	file1.close()

	##########################
	#Charges nodes positions
	pos={}
	file1 = open(city+'_pos.txt','r')
	for x in file1:
		x=x[:-2].split('\t')
		#print(x)
		if x[0]!='Node':
			#print(x)
			x=list(map(float,x))
			pos[int(x[0])]=(x[1],x[2])
			
	file1.close()


	os.chdir('..')
	os.chdir('..')
	return DG,SG,SG_PH,pos


def run_instances(city,alpha,n_it,new_cont='w'):
	os.chdir(r'Networks/'+city)
	
	DG,SG,SG_PH,pos=creates_graphs(city)
	os.chdir(r'Networks/'+city)
	file1 = open(city+'_Instances'+'.txt','r')
	os.chdir(r'Results')
	file2 = open(city+'_pulse'+'.txt',new_cont)
	file3 = open(city+'_MC_Not_corr'+'.txt',new_cont)
	file4= open(city+'_PH_3'+'.txt',new_cont)
	if new_cont=='w':
		file2.write('Time(s)\tCost(PB)\tTime\tP(T<T_max)\tBound\tInfeasibility\tDominance\tPath\n')
		file3.write('Time(s)\tCost(PB)\tTime\tP(T<T_max)\tBound\tInfeasibility\tDominance\tPath\n')		
	
	for x in file1:
		x=x[:-1].split('\t')
		#print(x)
		if x[0]!='Source':
			x=list(map(float,x))
			source, target, T_max=int(x[0]),int(x[1]),x[2]

			PG=pulse_graph(G=DG,T_max=T_max,source=source, target=target,tightness=0,pos=pos)
			SPG=MC.s_pulse_graph(SG,T_max=T_max,alpha=alpha,source=source,target=target,n_it=n_it,tightness=0,pos=pos)
			SPG_PH=PH_pulse.s_pulse_graph(SG_PH,T_max=T_max,alpha=alpha,source=source,target=target,pos=pos)

			###################################################
			#Runs CSP with pulse algo
			print("Running CSP")
			t_pulse=time.time()
			pulse_path,pulse_cost,pulse_sol_time=PG.run_pulse()
			t_pulse=time.time()-t_pulse
			
			#The probability of arriving on time is estimated using montecarlo
			#print(T_max)
			path=list(zip(pulse_path[:-1],pulse_path[1:]))
			tmin=sum(SPG.G[i][j]["tmin"] for i,j in path)
			#pulse_prob= SPG.monte_carlo(pulse_path,tmin,0,10000)
			os.chdir('..')
			os.chdir('..')
			os.chdir('..')
			pulse_prob=eval_path(path1=pulse_path,city=city,n_sim=10000,pT_max=T_max)
			os.chdir(r'Networks/'+city)
			os.chdir(r'Results')
			#print("Probabilidad de llegar a tiempo: ", pulse_prob)
			text= str(source) +'\t'+str(target)+'\t'+ str(t_pulse) +'\t'+str(pulse_cost)+"\t"+str(pulse_sol_time)+'\t'+str(pulse_prob)+'\t'+str(PG.bound)+'\t'+str(PG.infeas)+'\t'+str(PG.dom)+"\t"+str(pulse_path)+'\n'
			print(text)
			file2.write(text)

			###################################################
			#Runs CCSP with MC pulse algo
			
			print("Running CSP MC")
			t_s_pulse=time.time()
			pulse_path,pulse_cost,pulse_sol_time=SPG.run_pulse(t_limit=1000)
			t_s_pulse=time.time()-t_s_pulse

			path=list(zip(pulse_path[:-1],pulse_path[1:]))
			tmin=sum(SPG.G[i][j]["tmin"] for i,j in path)
			Exp_time=sum(SPG.G[i][j]["Time"] for i,j in path)
			os.chdir('..')
			os.chdir('..')
			os.chdir('..')
			pulse_prob=eval_path(path1=pulse_path,city=city,n_sim=10000,pT_max=T_max)
			os.chdir(r'Networks/'+city)
			os.chdir(r'Results')

			text= str(source) +'\t'+str(target)+'\t'+ str(t_s_pulse) +'\t'+str(pulse_cost)+"\t"+str(Exp_time)+'\t'+str(pulse_prob)+'\t' +str(SPG.Resource)+'\t'+str(SPG.bound)+'\t'+str(SPG.infeas)+'\t'+str(SPG.dom)+"\t"+str(pulse_path)+'\n'
			print(text)
			print("la prob del mc es ",SPG.Resource)
			print("Tiempo en Simulacion ",SPG.time_expm)
			file3.write(text)
			

			'''			
			###################################################
			#Runs CCSP with PH pulse algo
			print("Running CSP PH")
			t_s_pulse=time.time()
			pulse_path,pulse_cost,pulse_sol_time=SPG_PH.run_pulse(t_limit=1000)
			t_s_pulse=time.time()-t_s_pulse

			path=list(zip(pulse_path[:-1],pulse_path[1:]))
			tmin=sum(SPG.G[i][j]["tmin"] for i,j in path)
			Exp_time=sum(SPG.G[i][j]["Time"] for i,j in path)
			os.chdir('..')
			os.chdir('..')
			os.chdir('..')
			pulse_prob=eval_path(path1=pulse_path,city=city,n_sim=10000,T_max=T_max)
			os.chdir(r'Networks/'+city)
			os.chdir(r'Results')

			text= str(t_s_pulse) +'\t'+str(pulse_cost)+"\t"+str(Exp_time)+'\t'+str(pulse_prob)+'\t'+str(SPG.bound)+'\t'+str(SPG.infeas)+'\t'+str(SPG.dom)+"\t"+str(pulse_path)+'\n'
			print(text)
			print("Tiempo en Simulacion ",SPG.time_expm)
			file4.write(text)
			'''
			
	file1.close()
	file2.close()
	file3.close()
	file4.close()

	os.chdir('..')
	os.chdir('..')
	os.chdir('..')
	
#############################################################################
#Plays god: Simulates and evaluates paths
def eval_path(path1,city,n_sim,pT_max):
	
	G,pos=create_net(city)	
	
	pat=list(zip(path1[:-1],path1[1:]))
	data=np.array([0.0 for i in range(n_sim)])
	tmin=0
	'''print('Hola, voy a evaluar ',pat)
				print('con tmax ',pT_max)	'''
	for i,j in pat:
		mean=abs(G[i][j]["Time"]-G[i][j]["fftt"])+0.00001			
		phi=math.sqrt(G[i][j]["SD"]**2+ mean**2)
		mu=math.log(mean**2/phi)
		sigma=math.sqrt(math.log(phi**2/(mean**2)))
		data += np.random.lognormal(mu,sigma, n_sim)
		
		tmin+=math.floor( G[i][j]["fftt"])
	
	
	
	ind=list(map(lambda x: x<=pT_max-tmin,data))
	prob= sum(ind)/n_sim
	
	if prob==0 and len(pat)!=0:
		plt.axvline(x=pT_max-tmin,ymin=0, ymax=1000,color="lime")
		plt.hist(data)
		print("El tmin calculado es ",tmin)
		#plt.show()
		plt.clf()

	if len(pat)==0:
		prob=0
	return prob
	#print("sd ",G[i][j]["SD"])
	#print("mean " ,mean)
	#print("s/e ",G[i][j]["SD"]**2/(mean**2))
	#mu=math.log(mean)-1/2*math.log(G[i][j]["SD"]**2/(mean**2) ) 
	#sigma=math.sqrt(math.log(G[i][j]["SD"]**2/(mean**2) ))	


##############################################################################
#Reads results from txx and makes plots and computes statistics

def sum_up_fits(city):
	DG,SG,pos=creates_graphs(city)
	d_names = ['laplace','norm','lognorm', 'powerlognorm' ,'exponweib', 'pareto','johnsonsu','burr']
	frecs=[0 for i in d_names]
	for i,j in SG.edges():
		for k,d in enumerate(d_names):
			try:
				if SG[i][j]["distr"]==d:
					frecs[k]+=1
			except Exception as e:
				print('error con ',(i,j))
			

	plt.bar(d_names,frecs)
	plt.title("Fitted distributions in "+city)
	plt.show()

def read_results(city,nPhases):
	G,pos=create_net(city)
	os.chdir(r'Networks/'+city+'/Results')
	
	file1 = open(city+'_pulse'+'.txt','r')
	file2 = open(city+'_MC_Not_corr'+'.txt','r')
	file3 = open(city+'_PHfit_'+str(nPhases)+'.txt','r')
	sols_p=[]
	info=[]
	for x in file1:
		sol={}
		x=x.split('\t')
		if x[0]=='Time(s)':
			info=x
		else:
			for i,j in enumerate(info[:-1]):
				sol[j]=float(x[i])
			x[-1]=x[-1].replace('\n','')
			sol[info[-1]]=list(map(float,x[-1][1:-1].split(', ')))
			sols_p.append(sol)
	
	sols_MC_NC=[]
	for x in file2:
		sol={}
		x=x.split('\t')
		if x[0]=='Time(s)':
			info=x
		else:
			for i,j in enumerate(info[:-1]):
				sol[j]=float(x[i])
			x[-1]=x[-1].replace('\n','')
			
			#print("path: ",x[-1][1:-2].split(', '))
			try:
				sol[info[-1]]=list(map(float,x[-1][1:-1].split(', ')))
			except:
				sol[info[-1]]=[]
			sols_MC_NC.append(sol)

	sols_PH=[]
	for x in file3:
		sol={}
		x=x.split('\t')
		if x[0]=='Time(s)':
			info=x
		else:
			for i,j in enumerate(info[:-1]):
				sol[j]=float(x[i])
			x[-1]=x[-1].replace('\n','')
			
			#print("path: ",x[-1][1:-2].split(', '))
			try:
				sol[info[-1]]=list(map(float,x[-1][1:-1].split(', ')))
			except:
				sol[info[-1]]=[]
			sols_PH.append(sol)
	file1.close()
	file2.close()
	os.chdir('Plots')
	PG=pulse_graph(G=G,T_max=0,source=1, target=1,tightness=0,pos=pos)
	inst=1
	
	print("Guarde ")
	print (len(sols_p),len(sols_MC_NC),len(sols_PH))
	for p,MC,PHsol in list(zip(sols_p,sols_MC_NC,sols_PH)):
		PG.Fpath=p[info[-1]]
		##fig=PG.draw_graph(path2=MC[info[-1]],path3=[[info[-1]]])
		print(PHsol)
		fig=PG.draw_graph(path2=MC[info[-1]],path3=list(map(lambda x:x+1,PHsol[info[-1]])))
		#plt.legend(['CSP (Expected value)','CCSP (Montecarlo Independent)'])
		ax = fig.add_subplot(1,1,1)
		x=sum(p[0] for p in pos.values())/len(pos.values())
		y=sum(p[1] for p in pos.values())/len(pos.values())
		ax.plot([x],[y],color="lime",label='CSP (Expected value)')
		ax.plot([x],[y],color="blue",label='CCSP (Montecarlo)')
		ax.plot([x],[y],color="red",label='CCSP (Phase-Type)')
		plt.axis('off')
		fig.set_facecolor('w')
		plt.legend()
		fig.tight_layout()
		
		plt.savefig(city+'_'+str(inst)+'.png')

		#plt.show()
		inst+=1

	os.chdir('..')


	os.chdir('..')	
	os.chdir('..')
	os.chdir('..')


def eval_paths_ph(city,phases):
	os.chdir(r'Networks/'+city)
	file1 = open(city+'_Instances'+'.txt','r')


	os.chdir(r'Results')
	file = open(city+'_PHfit_'+str(phases)+'.txt','r')
	file2 = open(city+'_PHfit_'+str(phases)+'eval'+'.txt','w')

	instances=[]
	for x in file1:
		x=x[:-1].split('\t')
		
		if x[0]!='Source':
			x=list(map(float,x))
			#source target, T_max
			instances.append([int(x[0]),int(x[1]),x[2]])
			#print('esto',[int(x[0]),int(x[1]),x[2]])
	
	info=[]
	n=0
	for x in file:
		x=x[:-1].split('\t')
		x=list(filter(lambda y: y!='',x))
		if x!=[]:
			
			path=x[-1]
			x=list(map(float,x[:-1]))

			path=path.replace("[",'')
			path=path.replace("]",'')
			
			if path=='':
				prob=0
				path=[]
			else:
				path=list(map(lambda x:int(x)+1,path.split(', ')))
				

				os.chdir('..')
				os.chdir('..')
				os.chdir('..')
				print(instances[n])
				print(path)
				prob=eval_path(path,city,10000,instances[n][2])
				#print(prob)
				os.chdir(r'Networks/'+city)
				os.chdir(r'Results')
			
			if prob==0 and path!=[]:
				print(n,"T de ",instances[n][2])
				print('el t min que viene es de ',x[4])
				#print(n,prob,path)
			n+=1
			#print(x[:5]+[prob]+x[5:]+[path])
			info.append(x[:5]+[prob]+x[5:]+[path])
	

	for lin in info:
		text=str(lin[:-1]).replace(', ','\t')
		text=text.replace('[','')
		text=text.replace(']','')
		path=str(lin[-1])
		
		file2.write(text+ '\t' +path+'\n')
	os.chdir('..')	
	os.chdir('..')
	os.chdir('..')


def eval_insts(city):
	
	os.chdir(r'Networks/'+city)	
	file1 = open(city+'_Instances'+'.txt','r')
	instances=[]
	for x in file1:
		x=x[:-1].split('\t')
		
		if x[0]!='Source':
			x=list(map(float,x))
			#source target, T_max
			instances.append([int(x[0]),int(x[1]),x[2]])

	file1.close()

	DG,SG,SG_PH,pos=creates_graphs(city)
	info=[]
	for inst in instances:
		sp_time=nx.shortest_path(DG,inst[0],inst[1],"Time")
		sp_cost=nx.shortest_path(DG,inst[0],inst[1],"Cost")
		
		LET=sum(DG[i][j]["Time"] for i,j in zip(sp_time[:-1],sp_time[1:]))
		time_sc=sum(DG[i][j]["Time"] for i,j in zip(sp_cost[:-1],sp_cost[1:]))
		info.append(inst+[LET,time_sc])

	
	os.chdir(r'Networks/'+city)	
	file2 = open(city+'_Instances_mod'+'.txt','w')
	for lin in info:
		text=str(lin).replace(', ','\t')
		text=text.replace('[','')
		text=text.replace(']','')
			
		file2.write(text+'\n')
	os.chdir('..')
	os.chdir('..')
		

		

	




##############################################################################
#Use the functions to create files, run insances and save results.


############################################################################
#creates random instances for all the cities


'''

n_inst=10
clusters=4
tightness=0.93
inf_tightness=0.95
search_time=10
cities=["SiouxFalls",'Chicago-Sketch','GoldCoast','Philadelphia','Sydney']
cities=[cities[0]]

print(cities)
for city in cities:
	G,pos=create_net(city,True)
	#print(pos)
	print(city)
	print(len(G.nodes))
	print(len(G.edges))
	print('emoece con ',city)
	create_random_instances(city,G,clusters,pos,n_inst,tightness,inf_tightness,search_time,"a",False)
'''

############################################################################
#Fits the distributions to all cities
#cities=['Chicago-Sketch','GoldCoast','Philadelphia','Sydney']

'''


cities=["SiouxFalls",'Chicago-Sketch','GoldCoast','Philadelphia','Sydney']
cities=cities[2:4]



for city in cities:
	G,pos=create_net(city)
	print('emoece con ',city)
	fit_and_print(G,city,Continue=False)


'''

############################################################################
#Fits the Phase-type distributions to all cities


'''
cities=["SiouxFalls",'Chicago-Sketch','GoldCoast','Philadelphia','Sydney']
cities=cities[1:2]
nPhases=10

for city in cities:
	G,pos=create_net(city)
	print(len(G.nodes))
	print(len(G.edges))
	print('emoece con ',city)
	fit_PH_and_print(G,city,nPhases,False)


'''

##########################################################################
#runs the instances for all the cities 

'''

cities=["SiouxFalls",'Chicago-Sketch','GoldCoast','Philadelphia','Sydney']
cities=cities[1:2]
print(cities)
alpha=0.90
n_it=1000

for city in cities:
	print('empece con: ',city)
	run_instances(city,alpha,n_it)

'''
##########################################################################
#Evaluates PH-paths
#Evaluates instances

'''
city='Chicago-Sketch'
phases= 5
G,pos=create_net(city,True)
nx.draw(G,pos,node_color="black",node_size=25,arrows=False)
plt.show()

eval_paths_ph(city,phases)
eval_insts(city)
'''

##########################################################################

##########################################################################


if False:
	#Plot results from all cities
	#g=w

	city='Chicago-Sketch'
	plot_net(city)
	os.chdir(r'Networks/'+city)
		
	DG,SG,SG_PH,pos=creates_graphs(city)
	G=nx.DiGraph()
	for i,j in SG_PH.edges():
		G.add_edge(i,j,Time=SG_PH[i][j]["Time"],Cost=SG_PH[i][j]["Cost"],tmin=SG_PH[i][j]["tmin"],T=SG_PH[i][j]["tRV"].T.tolist(),alpha=SG_PH[i][j]["tRV"].alpha.tolist())

	nx.write_gpickle(G, r"C:\Users\David Corredor M\OneDrive - Universidad de Los Andes\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse\Notebooks\siouxfalls.gpickle")

	#print(pos)


	cities=["SiouxFalls",'Chicago-Sketch','GoldCoast','Philadelphia','Sydney']
	cities=cities[1:2]
	for city in cities:
		read_results(city,5)
		sum_up_fits(city)

