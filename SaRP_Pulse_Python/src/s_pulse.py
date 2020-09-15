import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from matplotlib import cbook, docstring, rcParams
import pylab
import matplotlib.animation
import math
import time as tiemm
import PhaseType
from PhaseType	import *
from importlib import reload
from mpmath import *
import time

reload(PhaseType)
#print("hola1")

#Stochastic pulse class
class s_pulse_graph():
	"""docstring for pulse_graph"""
	def __init__(self, G,T_max,alpha,source,target,pos=[]):
		#super(s_pulse_graph, self).__init__()
		self.G = G
		self.Primal_Bound=float("inf")
		self.Fpath=[]
		self.Resource=0
		self.minimum_time=0
		self.tRV=PH(matrix([[]]),matrix([[]]))
		self.T_max=T_max
		self.alpha=alpha
		self.source=source
		self.target=target
		self.anim=[]
		self.pos=pos
		self.time_expm=0
		self.bound=0
		self.infeas=0
		self.dom=0
		self.time_limit=float('inf')
		self.pulse_time=0
		self.feasibel=False

	#Preprocess the graph labeling every node with the shortest cost to the end node (target)
	def preprocess(self):
		p=nx.shortest_path_length(self.G,weight="Cost",target=self.target)
		attrs={i:{"labels":[],"s_cost":p[i]} for i in p.keys()}
		attrs.update({i:{"labels":[],"s_cost":float("inf")} for i in self.G.nodes() if i not in p.keys()})
		nx.set_node_attributes(self.G, attrs)

		p=nx.shortest_path_length(self.G,weight="tmin",target=self.target)
		attrs={i:{"s_time":p[i]} for i in p.keys()}
		attrs.update({i:{"s_time":float("inf")} for i in self.G.nodes() if i not in p.keys()})
		nx.set_node_attributes(self.G, attrs)
		
		self.minimum_time=attrs[self.source]["s_time"]

		for i in self.G.nodes:			
			self.G.nodes[i]["labels"]=[]		

	#Check Dominance
	def C_Dominance(self,vk,c,prob):
		bool=True	#We assume that the current path passes the dominance check (i.e is not dominated)
		for (i,(cc,probb)) in enumerate(self.G.nodes[vk]["labels"]):
			if c>cc and prob<probb:				
				bool=False
				self.dom+=1
		return bool

	#Check Feasibility
	def C_Feasibility(self,vk,prob):
		bool=True
		if prob<self.alpha:			
			bool=False
			self.infeas+=1
		return bool
	#Check Bounds
	def C_Bounds(self,vk,c):
		bool=True
		if c+self.G.nodes[vk]["s_cost"]>self.Primal_Bound:		
			bool=False
			self.bound+=1
		return bool
	#Update the labels of a given node vk
	def update_labels(self,vk,c,prob):	
		self.G.nodes[vk]["labels"].append((c,prob))
	
	def sort(self,sons):
		if self.feasibel:
			return(sorted(sons,key=lambda x: self.G.nodes[x]["s_cost"] ))
		else:
			return(sorted(sons,key=lambda x: self.G.nodes[x]["s_time"] ))

	def update_primal_bound(self,c,prob,P,tRV):
		if prob>=self.alpha and c<=self.Primal_Bound:
			self.Primal_Bound=c
			self.Fpath=P
			self.Resource=prob
			self.tRV=tRV
			self.feasibel=True
			#print("Primal_Bound: ",c, prob)
	
	def pulse(self,vk,c,tRV,tmin,ExptT,P):		
		cota=(ExptT-tmin)/(self.T_max-tmin)	#Makkovitz inequality
		if (self.C_Bounds(vk,c) and vk not in P):
			
			if cota<1-self.alpha:
				prob=1-cota
				#print("Se usa la cota de ", prob )
			else:				
				t=time.time()
				#print("T_max ",self.T_max)
				#print("t_min ",tmin)
				#print("t_min_in ",self.G.nodes[vk]["s_time"])
				#print("Constr ",self.T_max-tmin-self.G.nodes[vk]["s_time"])		
				#print("NO se usa la cota de ", 1-cota )
				prob=tRV.cdf(self.T_max-tmin-self.G.nodes[vk]["s_time"])	
				self.time_expm+=time.time()-t
			
			#print("Entro a ",vk, "con ", c,prob ,tmin,P)
			self.update_labels(vk,c,prob)		
			if vk==self.target:
				self.update_primal_bound(c,prob,P+[vk],tRV)
			elif self.C_Feasibility(vk,prob) and self.C_Dominance(vk,c,prob) :
				PP=P.copy()
				PP.append(vk)
				for i in self.sort(self.G.successors(vk)):			
					cc=c+self.G.edges[vk,i]["Cost"]
					
					ntRV=tRV.sum(self.G.edges[vk,i]["tRV"])
					ntmin=tmin+self.G.edges[vk,i]["tmin"]
					nExptT=ExptT+self.G.edges[vk,i]["Time"]
					self.pulse(i,cc,ntRV,ntmin,nExptT,PP)
	
	def run_pulse(self,t_limit=1000):
		self.time_limit=t_limit
		self.Fpath=[]
		self.Primal_Bound=float("inf")
		self.Resource=0
		#print(self.Fpath)
		self.preprocess()
		#print("s",self.source)
		if self.G.nodes[self.source]["s_cost"]!=float("inf"):
			self.pulse_time=time.time()
			self.pulse(vk=self.source,c=0,tRV=self.tRV,tmin=0,ExptT=0,P=self.Fpath)
			self.pulse_time=time.time()-self.pulse_time
		else:
			print("The instance is infeasible")

		return self.Fpath, self.Primal_Bound,self.Resource

	def prob_path(self,path, T_max):
		RV=self.G[path[0]][path[1]]["tRV"]
		arcs=list(zip(path[1:-1],path[2:]))

		tmin=0

		for i,j in arcs:
			RV=RV.sum(self.G[i][j]["tRV"])
			tmin+=self.G[i][j]["tmin"]
		

		return(RV.cdf(T_max-tmin))
	
class graphDrawing(s_pulse_graph):
	"""docstring for graphDrawing"""
	def __init__(self):		
		pass		
		

	#Draws the graph with a given position and a given path
	def draw_graph(self,path=[],bgcolor="white",edge_color="black",arc_color="gray",path_color="red"):
		if self.pos==[]:
			self.pos = nx.random_layout(self.G)
		if path==[]:
			path=self.Fpath

		fig= plt.figure(figsize=(8,6))
		#edge_labels={e:(int(self.G.edges[e]["Cost"])) for e in self.G.edges}
		BGnodes=set(self.G.nodes()) - set(path)
		nx.draw_networkx_edges(self.G, pos=self.pos, edge_color=arc_color)
		null_nodes = nx.draw_networkx_nodes(self.G, pos=self.pos, nodelist=BGnodes, node_color=bgcolor)#node_size=1000
		null_nodes.set_edgecolor(edge_color)
		nx.draw_networkx_labels(self.G, pos=self.pos, labels=dict(zip(BGnodes,BGnodes)),  font_color="black",font_size=7)
		try:
			query_nodes=nx.draw_networkx_nodes(self.G,pos=self.pos,nodelist=path,node_color=path_color)#,node_size=1000
			query_nodes.set_edgecolor(path_color)
		except:
			pass
		nx.draw_networkx_labels(self.G,pos=self.pos,labels=dict(zip(path,path)),font_color="black",font_size=7)
		edgelist = [path[k:k+2] for k in range(len(path) - 1)]
		nx.draw_networkx_edges(self.G, pos=self.pos, edgelist=edgelist, width=4,edge_color=path_color)

	    #Edge labels
		#nx.draw_networkx_edge_labels(self.G, self.pos, edge_labels = edge_labels ,font_size=7)

		plt.axis('off')
		#plt.figure(num=None, figsize=(6, 3), dpi=80*3, facecolor='w', edgecolor='k')
		plt.show()
	
	def pulse_anim(self,vk,c,tRV,tmin,P):
		#print(tmin)
		prob=tRV.cdf(self.T_max-tmin)
		self.update_labels(vk,c,prob)
		if vk==self.target:
			self.update_primal_bound(c,prob,P+[vk],tRV)
			if prob>=self.alpha and c<=self.Primal_Bound:
				self.anim.append(["Primal_Bound",P+[vk],c,round(prob,4),self.Fpath,self.Primal_Bound,self.Resource])

		elif not self.C_Dominance(vk,c,prob):
			self.anim.append(["Dom",P+[vk],c,round(prob,4),self.Fpath,self.Primal_Bound,self.Resource])
		elif not self.C_Feasibility(vk,prob):
			self.anim.append(["Feas",P+[vk],c,round(prob,4),self.Fpath,self.Primal_Bound,self.Resource])
		elif not self.C_Bounds(vk,c):
			self.anim.append(["Bound",P+[vk],c,round(prob,4),self.Fpath,self.Primal_Bound,self.Resource])
		elif vk in P:
			self.anim.append(["Loop",P+[vk],c,round(prob,4),self.Fpath,self.Primal_Bound,self.Resource])

		if (self.C_Dominance(vk,c,prob) and self.C_Feasibility(vk,prob) and self.C_Bounds(vk,c) and vk not in P and vk!=self.target):
			PP=P.copy()
			PP.append(vk)
			for i in self.sort(self.G.successors(vk)):			
				cc=c+self.G.edges[vk,i]["Cost"]
				ntRV=tRV.sum(self.G.edges[vk,i]["tRV"])
				probi=ntRV.cdf(self.T_max)
				ntmin=tmin+self.G.edges[vk,i]["tmin"]
				#self.pulse(i,cc,ntRV,ntmin,PP)
				if i not in P:
					self.anim.append([PP+[i],cc,round(probi,4),self.Fpath,self.Primal_Bound,self.Resource])
					self.pulse_anim(i,cc,ntRV,ntmin,PP)

	def onClick(self,event):
		global pause
		pause ^= True
		if pause:
			ani.event_source.stop()
		elif not pause:
			ani.event_source.start()
	def update(self,num):
		edge_labels= {e:(self.G.edges[e]['Cost']) for e in self.G.edges}
		ax.clear()
		axCurrent_p.clear()
		axPrimal_B.clear()
		
		axCurrent_p.set_facecolor('white')
		axPrimal_B.set_facecolor('white')

		fontsize=11
		prune=False
		primB=False
		opt=False
		
		y=(1/len(self.pos.values()))*sum(yy for xx,yy in self.pos.values())
		x=0.05
		dy=100000/2

		#ax.text(-1.2, 1,'Stop/Start \n button',fontsize=fontsize-5,fontweight="bold",color='white',bbox={'facecolor': 'green', 'alpha': 1, 'pad': 5})
		if not isinstance(self.anim[num][0], str):
			path = self.anim[num][0]
			c=self.anim[num][1]
			t=self.anim[num][2]
			FP=self.anim[num][3]
			PB=self.anim[num][4]
			RC=self.anim[num][5]
		elif self.anim[num][0]=="Primal_Bound":
			
			primB=True
			axPrimal_B.set_xticks([])
			axPrimal_B.set_yticks([])
			axPrimal_B.set_title("Best Bound information",fontweight="bold")
			axPrimal_B.text(x, y+4*dy, 'New integer solution found!',fontsize=fontsize-2,fontweight="bold")
			path=self.anim[num][1]
			c=self.anim[num][2]
			t=self.anim[num][3]
			FP=self.anim[num][4]
			PB=self.anim[num][5]
			RC=self.anim[num][6]
			cp=""
			for i in FP:
				cp+=str(i)+"-"
			axPrimal_B.text(x, y+2*dy, cp[:-1],fontsize=fontsize)
			axPrimal_B.text(x, y+dy, "Primal Bound:",fontsize=fontsize)
			axPrimal_B.text(x, y, str(PB),fontsize=fontsize)
			axPrimal_B.text(x, y-1*dy, "Probability:",fontsize=fontsize)
			axPrimal_B.text(x, y-2*dy, str(round(t,4)),fontsize=fontsize)
		elif self.anim[num][0]=="Opt":
			opt=True
			axPrimal_B.set_xticks([])
			axPrimal_B.set_yticks([])       
			axPrimal_B.set_title("Best Bound information",fontweight="bold")
			axPrimal_B.text(x, y+4*dy, 'Optimal solution found!',fontsize=fontsize-2,fontweight="bold")
			path=self.anim[num][1]
			PB=self.anim[num][2]
			RC=self.anim[num][3]
			cp=""
			for i in path:
				cp+=str(i)+"-"
			axPrimal_B.text(x, y+2*dy, cp[:-1],fontsize=fontsize)
			axPrimal_B.text(x, y+1*dy, "Primal Bound:",fontsize=fontsize)
			axPrimal_B.text(x, y-0*dy, str(PB),fontsize=fontsize)
			axPrimal_B.text(x, y-1*dy, "Probability:",fontsize=fontsize)
			axPrimal_B.text(x, y-2*dy, str(round(RC,4)),fontsize=fontsize)
		else:
			
			if self.anim[num][0]=="Dom":
				prune=True
				axCurrent_p.text(x, y+4*dy, r'Pulse is discarted by dominance:' ,fontsize=fontsize-3)
				axCurrent_p.text(x, y+2*dy, r'$c(\mathcal{P}´)$     $P[t(\mathcal{P}´)>T]$' ,fontsize=fontsize-1)

				path=self.anim[num][1]
				c=self.anim[num][2]
				t=self.anim[num][3]
				FP=self.anim[num][4]
				PB=self.anim[num][5]
				RC=self.anim[num][6]
				axCurrent_p.text(x+0.1, y+1*dy, str(c)+'        '+str(t) ,fontsize=fontsize-1)
				tit="Labels at node "+str(path[-1])
				cost=""
				tiemp=""
				for (i,j) in self.G.nodes[path[-1]]["labels"]:
					cost+=str(int(i))+"    "
					tiemp+=str(round(j,2))+"    "
				axCurrent_p.text(x, y-2*dy,tit+'\n\n'+'Cost: '+cost+'\n'+'Resource: '+tiemp,fontsize=fontsize-5,bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5})
				#plt.pause(3)


			elif self.anim[num][0]=="Feas":
				prune=True
				path=self.anim[num][1]
				c=self.anim[num][2]
				t=self.anim[num][3]
				FP=self.anim[num][4]
				PB=self.anim[num][5]
				RC=self.anim[num][6]
				prune=True
				axCurrent_p.text(x, y+4*dy, 'Pulse is discarted by Feasibility',fontsize=fontsize-3,color='black'),
				axCurrent_p.text(x+0.1, y+1*dy,r'$P[t (\mathcal{P}´)<T] <\alpha $',fontsize=fontsize,bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5})
				axCurrent_p.text(x+0.1, y-1*dy,str(t) +'<'+ str(self.alpha) ,fontsize=fontsize,bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5})
				#plt.pause(3)
			elif self.anim[num][0]=="Bound":
				
				axCurrent_p.text(x, y+4*dy, 'Pulse is discarted by Bound',fontsize=fontsize-3,color='black')

				path=self.anim[num][1]
				c=self.anim[num][2]
				t=self.anim[num][3]
				FP=self.anim[num][4]
				PB=self.anim[num][5]
				RC=self.anim[num][6]

				axCurrent_p.text(x+0, y+1*dy,r'$c(\mathcal{P}´) \plus c \underbar (v_k) \geq \bar{c}$',fontsize=fontsize+2,bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5})
				axCurrent_p.text(x+0, y-1*dy,str(c)+"+"+str(self.G.nodes[path[-1]]["s_cost"])+r'$\geq$'+str(PB) ,fontsize=fontsize+2,bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5})

				prune=True
				#plt.pause(3)
			elif self.anim[num][0]=="Loop":
				
				axCurrent_p.text(x, y-.8, 'Pulse is discarted by loop',fontsize=fontsize,color='black')
				path=self.anim[num][1]
				c=self.anim[num][2]
				t=self.anim[num][3]
				FP=self.anim[num][4]
				PB=self.anim[num][5]
				RC=self.anim[num][6]
				prune=True
		# Background nodes
		BGnodes=set(self.G.nodes()) - set(path)
		nx.draw_networkx_edges(self.G, pos=self.pos, ax=ax, edge_color="royalblue")
		null_nodes = nx.draw_networkx_nodes(self.G, pos=self.pos, nodelist=BGnodes, node_color="white",  ax=ax)
		null_nodes.set_edgecolor("orange")
		nx.draw_networkx_labels(self.G, pos=self.pos, labels=dict(zip(BGnodes,BGnodes)),  font_color="black", ax=ax)
		try:
			query_nodes = nx.draw_networkx_nodes(self.G, pos=self.pos, nodelist=path, node_color="red", ax=ax)
			query_nodes.set_edgecolor("red")
		except:
			pass
		nx.draw_networkx_labels(self.G, pos=self.pos, labels=dict(zip(path,path)),  font_color="black", ax=ax)
		edgelist = [path[k:k+2] for k in range(len(path) - 1)]
		nx.draw_networkx_edges(self.G, pos=self.pos, edgelist=edgelist, width=4,edge_color="red", ax=ax)

	    #Edge labels
		#nx.draw_networkx_edge_labels(self.G, self.pos, edge_labels = edge_labels ,ax=ax)
	    # Scale plot ax
		ax.set_title("Graph", fontweight="bold")
		ax.set_xticks([])
		ax.set_yticks([])
		#Current path info}
		axCurrent_p.set_xticks([])
		axCurrent_p.set_yticks([])
		axCurrent_p.set_title("Current path information ",fontweight="bold")
		if prune and not opt:
			if opt:
				print("entre aca 1")
			axCurrent_p.set_facecolor('red')
		elif not prune and not opt:
			if opt:
				print("entre aca 2")
			axCurrent_p.set_facecolor('white')
			axCurrent_p.text(.05, y+3*dy, 'Current Path ' + r'$( \mathcal{P} )$:',fontsize=fontsize)
			cp=""
			for i in path:
				cp+=str(i)+"-"
			axCurrent_p.text(0.05, y+2*dy, cp[:-1],fontsize=fontsize)
			axCurrent_p.text(0.05, y+dy, "Current Cost:",fontsize=fontsize)
			axCurrent_p.text(0.05, y-dy, str(c),fontsize=fontsize)
			axCurrent_p.text(.05, y-2*dy, "Probability:",fontsize=fontsize)
			axCurrent_p.text(.05, y-3*dy, str(t),fontsize=fontsize)
			#########################################################################################################################
		#Primal bound info
		if primB and not opt:
			axPrimal_B.set_facecolor('lime')
		elif not primB and not opt:
			axPrimal_B.set_xticks([])
			axPrimal_B.set_yticks([])
			axPrimal_B.set_title("Best Bound information",fontweight="bold")

			axPrimal_B.text(.05, y+4*dy, 'Best Path ' + r'$( \mathcal{P}* )$:',fontsize=fontsize)
			cp=""
			for i in FP:
				cp+=str(i)+"-"
			axPrimal_B.text(0.05, y+2*dy, cp[:-1],fontsize=fontsize)
			axPrimal_B.text(0.05, y+dy, "Primal Bound:",fontsize=fontsize)
			axPrimal_B.text(0.05, y, str(PB),fontsize=fontsize)
			axPrimal_B.text(.05, y-dy, "Probability:",fontsize=fontsize)
			axPrimal_B.text(.05, y-2*dy, str(round(RC,4)),fontsize=fontsize)
			global pause
		if prune:
			ax.text(1.6,-.30,"Press here to continue",color="white",fontweight="bold",bbox={'facecolor': 'green', 'alpha': 1, 'pad': 5})
			if sp:
				pause^=True
				ani.event_source.stop()

		if primB:
			ax.text(1.6,-.30,"Press here to continue",color="white",fontweight="bold",bbox={'facecolor': 'green', 'alpha': 1, 'pad': 5})
			if sp:
				pause^=True
				ani.event_source.stop()
				plt.pause(5)
		if opt:
			axPrimal_B.set_facecolor('lime')
			plt.pause(10)
		return fig

	def animation(self,speed,stop_prune=False,pos=[]):
		if pos==[]:
			pos = nx.random_layout(self.G)
		c=0
		tRV=PH(matrix(np.array([[]])),matrix(np.array([[]])))
		global sp,pause
		sp=stop_prune
		pause=False
		self.Primal_Bound=float("inf")
		self.Resource=0
		self.Fpath=[]
	    
		P=[]
		pause = False
		self.anim=[]
		self.preprocess()
		tmin=0
		self.pulse_anim(self.source,c,tRV,tmin,P)
		self.anim.append(["Opt",self.Fpath,self.Primal_Bound,self.Resource])

	    # Build plot
		global fig, ax,axCurrent_p,axPrimal_B,ani
		fig= plt.figure(figsize=(12*1.5,6))
		fig.canvas.mpl_connect('button_press_event', self.onClick)
		grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
		ax = fig.add_subplot(grid[:, :3])
		axCurrent_p = fig.add_subplot(grid[:2, 3], xticklabels=[],yticklabels=[], sharey=ax)
		axPrimal_B = fig.add_subplot(grid[2:, 3], xticklabels=[],yticklabels=[], sharey=ax)
		edge_labels= {e:(self.G.edges[e]['Cost']) for e in self.G.edges}
		ani = matplotlib.animation.FuncAnimation(fig, self.update, frames=len(self.anim), interval=speed, repeat=True,blit=False)
		return(ani)
	def draw_path(self, path,color="yellow",edge_color="black",arc_color="gray",):
		edgelist = list(zip(path,path[1:]))
		edge_labels={e:(self.G.edges[e]["Cost"]) for e in self.G.edges if e in edgelist}
		pos={i:(path.index(i),1) for i in path}

		fig= plt.figure(figsize=(12,6))
		
		query_nodes=nx.draw_networkx_nodes(self.G,pos=pos,nodelist=path,node_color=color,node_size=1000)
		query_nodes.set_edgecolor(edge_color)
		nx.draw_networkx_labels(self.G,pos=pos,labels=dict(zip(path,path)),font_color="black")
		
		#edgelist = [path[k:k+2] for k in range(len(path) - 1)]
		
		nx.draw_networkx_edges(self.G, pos=pos, edgelist=edgelist, width=4,edge_color=arc_color)
		#Edge labels
		nx.draw_networkx_edge_labels(self.G, pos=pos,edgelist=edgelist, edge_labels = edge_labels )
		plt.axis('off')
		plt.figure(num=None, figsize=(6, 3), dpi=80*3, facecolor='w', edgecolor='k')
		plt.show()

