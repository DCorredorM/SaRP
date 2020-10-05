import os,time,random
import numpy as np
import scipy.stats as st
import networkx as nx



global city,timeLimit

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

def createFolder(name):
	os.popen(f'[ -d {name} ]  || mkdir {name}')

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
	d=os.popen(f'java -Xmx1024m -jar PH_Pulse_SaRP.jar {folder} {city}')
	d=d.read().replace('\n','').replace('pk menor que Double.MIN_VALUE','')#.split('\t')	
	return d

def runInstancesIndependent(nPhases,nScen,alpha,clearF=True):
	'''
	Runs all the insances for the given city
	Args:
		city(String): folder name with the data
	Return:
		None: prints out a file with the results:
			source+"\t"+target+"\t"+ running time (s)+"\t"+PrimalBound+"\t"+Free flow time+"\t"+ reliability+"\t"+Bound+"\t"+Infeasibility+"\t"+Dominance+"\t"+finalPath
	'''
	#if clearF: clearResultFiles()
	
	folder=os.path.abspath(f'../../data/Networks/{city}')
	inst=open(f'{folder}/{city}_instances_{tightness}.txt','r')	
	
	createFolder(name=f'{folder}/Results/Independent/PHFit{nPhases}_{tightness}')
	
	scen=list(range(1,nScen+1))#+['Total']
	nn=0
	for l in inst:
		if l[0]!='#':
			i=l.replace('\n','').split('\t')
			for k in scen:
				config=open(f'{folder}/{city}_config.txt','w')
				print(i)
				text=f'PHFitFile:Independent/PHFit{nPhases}_scen{k}.txt\nDataFile:Independent/scen{k}.txt\nNumber of Arcs:2950\nNumber of Nodes:933\nTime Constraint:{float(i[2])}\nStart Node:{i[0]}\nEnd Node:{i[1]}\nNumber of Phases:{nPhases}\nalpha:{alpha}\nTime Limit:{timeLimit}'
				config.write(text)				
				config.close()
				d=runExperiment()
				resultFile=open(f'{folder}/Results/Independent/PHFit{nPhases}_{tightness}/scen{k}.txt','a')
				print(d)
				resultFile.write(d)
				resultFile.close()

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
	tightness=0.4
	nPhases=3
	########################################################		

	#runInstancesIndependent(nPhases=10,nScen=5,alpha=0.8)
	runInstancesIndependent(nPhases=5,nScen=1,alpha=0.8)
	#runInstancesIndependent(nPhases=3,nScen=5,alpha=0.8)

