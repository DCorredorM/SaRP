
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from matplotlib import cbook, docstring, rcParams

import pylab
import matplotlib.animation
from scipy.linalg import expm, sinm, cosm

import math
import time as tiemm


		

class pulse_graph():
	"""docstring for pulse_graph"""
	def __init__(self, G,T_max,source,target,tightness=0.5,pos=[]):
		super(pulse_graph, self).__init__()
		self.G = G
		self.Primal_Bound=float("inf")
		self.Fpath=[]
		self.Resource=0
		self.T_max=T_max
		self.source=source
		self.target=target
		self.anim=[]
		self.pos=pos
		self.bound=0
		self.infeas=0
		self.dom=0
		self.tightness=tightness

		
	def calc_cost(self,p):
		edges=zip(p,p[1:])
		cost=0
		time=0
		for (i,j) in edges:
			cost+=self.G[i][j]["Cost"]
			time+=self.G[i][j]["Time"]
		return (cost, time)
	#Preprocess the graph labeling every node with the shortest cost to the end node (target)
	def preprocess(self):

		p=nx.shortest_path_length(self.G,weight="Cost",target=self.target)
		attrs={i:{"labels":[],"s_cost":p[i]} for i in p.keys()}
		attrs.update({i:{"labels":[],"s_cost":float("inf")} for i in self.G.nodes() if i not in p.keys()})
		nx.set_node_attributes(self.G, attrs)

		t=nx.shortest_path_length(self.G,weight="Time",target=self.target)
		attrs={i:{"s_time":t[i]} for i in p.keys()}
		attrs.update({i:{"s_time":float("inf")} for i in self.G.nodes() if i not in p.keys()})
		

		nx.set_node_attributes(self.G, attrs)
		#self.minimum_time=attrs[self.source]["s_time"]

		if self.T_max==None:
			try:
				self.T_max=t[self.source]*(1+self.tightness)
			except:
				print("Infeasible") 
		
		for i in self.G.nodes:
			self.G.nodes[i]["labels"]=[]

	#Check Dominance
	def C_Dominance(self,vk,c,t):
		bool=True	#We assume that the current path passes the dominance check (i.e is not dominated)
		for (i,(cc,tt)) in enumerate(self.G.nodes[vk]["labels"]):
			if c>cc and t>tt:
				bool=False
				self.dom+=1
		return bool
	
	#Check Feasibility
	def C_Feasibility(self,vk,t):
		bool=True #We assume that the current path passes the Feasibility check (i.e is feasible)
		if t+self.G.nodes[vk]["s_time"]>self.T_max:
			bool=False
			self.infeas+=1
		return bool
	#Check Bounds
	def C_Bounds(self,vk,c):
		bool=True #We assume that the current path passes the primal_Bound check (i.e in the bes casenario the path is better than the PB)
		if c+self.G.nodes[vk]["s_cost"]>self.Primal_Bound:
			#print("bound")
			bool=False
			self.bound+=1
		return bool
	#Check path completion
	def path_completion(self,vk,c,t,P):
		bool=True #We assume that the current path passes the path_completion check (i.e is not possible to complete the path)
		if (c + self.G.nodes[vk]["s_cost"]<self.Primal_Bound) and (t+self.G.nodes[vk]["s_time"]<=self.T_max):
			#self.update_primal_bound(c + self.G.nodes[vk]["s_cost"],t+self.G.nodes[vk]["s_time"],P)
			bool=True
			#print("its working")
		return(bool)
	#Update the labels of a given node vk
	def update_labels(self,vk,c,t):	
		self.G.nodes[vk]["labels"].append((c,t))
	def sort(self,sons):
		return(sorted(sons,key=lambda x: self.G.nodes[x]["s_cost"] ))

	def update_primal_bound(self,c,t,P):
		if t<=self.T_max and c<=self.Primal_Bound:
			self.Primal_Bound=c
			self.Fpath=P
			self.Resource=t
			#print("Nuevo PB, costo: ",self.Primal_Bound,"tiempo: ",self.Resource)
	
	def pulse(self,vk,c,t,P):
		self.update_labels(vk,c,t)
		if vk==self.target:
			self.update_primal_bound(c,t,P+[vk])
			#print("LLegue a ",vk, "Con tiempo de ",t, " Y costo de ",c)
		if (self.C_Dominance(vk,c,t) and self.C_Feasibility(vk,t) and self.C_Bounds(vk,c) and self.path_completion(vk,c,t,P) and vk not in P):
			PP=P.copy()
			PP.append(vk)
			for i in self.sort(self.G.successors(vk)):			
				cc=c+self.G.edges[vk,i]["Cost"]
				tt=t+self.G.edges[vk,i]["Time"]
				self.pulse(i,cc,tt,PP)
	
	def run_pulse(self):
		self.Fpath=[]
		self.Primal_Bound=float("inf")
		self.Resource=0
		#print(self.Fpath)
		self.preprocess()
		if self.G.nodes[self.source]["s_cost"]!=np.Infinity:
			self.pulse(self.source,0,0,self.Fpath)
		else:
			print("The instance is not feasible")
		
		return self.Fpath, self.Primal_Bound,self.Resource

	
	#Draws the graph with a given position and a given path
	def draw_graph(self,path=[],bgcolor="white",edge_color="black",arc_color="gray",path_color="red"):
		if self.pos==[]:
			self.pos = nx.random_layout(self.G)
		if path==[]:
			path=self.Fpath

		fig= plt.figure(figsize=(12,6))
		edge_labels={e:(int(self.G.edges[e]["Cost"]),int(self.G.edges[e]["Time"])) for e in self.G.edges}
		BGnodes=set(self.G.nodes()) - set(path)
		nx.draw_networkx_edges(self.G, pos=self.pos, edge_color=arc_color)
		null_nodes = nx.draw_networkx_nodes(self.G, pos=self.pos, nodelist=BGnodes, node_color=bgcolor)#node_size=1000
		null_nodes.set_edgecolor(edge_color)
		nx.draw_networkx_labels(self.G, pos=self.pos, labels=dict(zip(BGnodes,BGnodes)),  font_color="black")
		try:
			query_nodes=nx.draw_networkx_nodes(self.G,pos=self.pos,nodelist=path,node_color=path_color)#,node_size=1000
			query_nodes.set_edgecolor(path_color)
		except:
			pass
		nx.draw_networkx_labels(self.G,pos=self.pos,labels=dict(zip(path,path)),font_color="black")
		edgelist = [path[k:k+2] for k in range(len(path) - 1)]
		nx.draw_networkx_edges(self.G, pos=self.pos, edgelist=edgelist, width=4,edge_color=path_color)

	    #Edge labels
		nx.draw_networkx_edge_labels(self.G, self.pos, edge_labels = edge_labels )
		plt.axis('off')
		#plt.figure(num=None, figsize=(6, 3), dpi=80*3, facecolor='w', edgecolor='k')
		#plt.show()	
		#return fig
		




'''
insts=["USA-road-BAY.txt","network10.txt"]

G=create_graph(insts[1], headder=False)
#PG=pulse_graph(G=G,T_max=3500,source=1,target=3301)
# 6236
#(4-5)*tiddes +tmin
PG=pulse_graph(G=G,T_max=9000,source=0,target=6236)

#pos=nx.spring_layout(G)
sol=PG.run_pulse()
print(sol)
#PG.draw_graph(path=sol[0])

'''




