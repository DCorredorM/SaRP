import os

global city,timeLimit

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
	d=os.popen(f'java -Xmx1024m -jar PH_Pulse_SaRP.jar {folder} {city}')
	d=d.read().replace('\n','').replace('pk menor que Double.MIN_VALUE','')#.split('\t')	
	return d

def runInstancesScenarios(nPhases,nScen,alpha,clearF=True):
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
	
	createFolder(name=f'{folder}/Results/Scenarios/PHFit{nPhases}_{tightness}')
	
	scen=list(range(1,nScen+1))#+['Total']
	nn=0
	for l in inst:
		if l[0]!='#':
			i=l.replace('\n','').split('\t')
			for k in scen:
				config=open(f'{folder}/{city}_config.txt','w')
				print(i)
				text=f'PHFitFile:Scenarios/PHFit{nPhases}_scen{k}.txt\nDataFile:Scenarios/scen{k}.txt\nNumber of Arcs:2950\nNumber of Nodes:933\nTime Constraint:{float(i[2])}\nStart Node:{i[0]}\nEnd Node:{i[1]}\nNumber of Phases:{nPhases}\nalpha:{alpha}\nTime Limit:{timeLimit}'
				config.write(text)				
				config.close()
				d=runExperiment()
				resultFile=open(f'{folder}/Results/Scenarios/PHFit{nPhases}_{tightness}/scen{k}.txt','a')
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
	folder=os.path.abspath(f'../../data/Networks/{city}/Results/Scenarios/*.txt')
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
	########################################################		

	#runInstancesScenarios(nPhases=10,nScen=5,alpha=0.8)
	#runInstancesScenarios(nPhases=5,nScen=5,alpha=0.8)
	for i in [0.2,0.6,0.8]:
		tightness=i
		runInstancesScenarios(nPhases=3,nScen=5,alpha=0.8)

