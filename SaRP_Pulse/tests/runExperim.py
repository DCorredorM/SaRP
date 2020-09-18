import os

global city,timeLimit

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
	d=os.popen(f'java -jar PH_Pulse_SaRP.jar {folder} {city}')
	d=d.read().replace('\n','').split('\t')	


def runInstances(nPhases,nScen,alpha,clearF=True):
	'''
	Runs all the insances for the given city
	Args:
		city(String): folder name with the data
	Return:
		None: prints out a file with the results:
			source+"\t"+target+"\t"+ running time (s)+"\t"+PrimalBound+"\t"+Free flow time+"\t"+ reliability+"\t"+Bound+"\t"+Infeasibility+"\t"+Dominance+"\t"+finalPath
	'''
	if clearF: clearResultFiles()
	
	folder=os.path.abspath(f'../../data/Networks/{city}')
	inst=open(f'{folder}/{city}_instances.txt','r')

	gamma=0.9

	scen=list(range(1,nScen+1))+['Total']
	nn=0
	for l in inst:
		if l[0]!='#':
			i=l.replace('\n','').split('\t')
			for k in scen[::-1]:
				config=open(f'{folder}/{city}_config.txt','w')
				print(i)
				text=f'DataFile:Scenarios/PHFit{nPhases}_scen{k}.txt\nNumber of Arcs:2950\nNumber of Nodes:933\nTime Constraint:{float(i[2])*gamma}\nStart Node:{i[0]}\nEnd Node:{i[1]}\nNumber of Phases:{nPhases}\nalpha:{alpha}\nTime Limit:{timeLimit}'
				config.write(text)
				config.close()
				runExperiment()
			nn+=1
			if nn>=1000:
				break

def clearResultFiles():
	'''
	clears all result files from city path
	'''
	folder=os.path.abspath(f'../../data/Networks/{city}/Results/Scenarios/*.txt')
	f=open(folder[:-5]+'n.txt','w')
	f.write('i')
	f.close()

	if os.name=='posix':
		os.popen(f'rm {folder}')
	else:
		os.popen(f'del {folder}')







if __name__ == '__main__':
	city='Chicago-Sketch'
	timeLimit=5000
	#runInstances(3,5,0.8)
	#print(os.name)
	
	runInstances(nPhases=10,nScen=5,alpha=0.8)
	runInstances(nPhases=5,nScen=5,alpha=0.8)
	runInstances(nPhases=3,nScen=5,alpha=0.8)

