#http://www.bgu.ac.il/~bargera/tntp/
#https://data.cityofchicago.org/Transportation/Chicago-Traffic-Tracker-Historical-Congestion-Esti/sxs8-h27x/data
import sys
from os import path

sys.path.append(path.abspath(r'C:\Users\d.corredor\OneDrive - Universidad de los andes\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse'))
#sys.path.insert(1, r'C:\Users\David Corredor M\OneDrive - Universidad de Los Andes\Thesis\2. Methodology\2.2 S-Pulse Algorithm\S-PULSE implementation\Code_S-pulse')
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

import pandas as pd
import csv
import datetime
from shapely.geometry import Point
#print(help(pd.read_csv))
#data =pd.read_csv(r'times.csv')
#print(data.head(10) )
#print(10)
def df_empty(columns, dtypes, index=None):
    assert len(columns)==len(dtypes)
    df = pd.DataFrame(index=index)
    for c,d in zip(columns, dtypes):
        df[c] = pd.Series(dtype=d)
    return df

def import_k_f(k):
	dtypes=['datetime64[ns, UTC]',int,float,str,str,str,str,float,str,str,int,int,int,int,int,str,float,float,float,float,'O','O']
	dtypesf=[lambda x:datetime.datetime.strptime(x, '%m/%d/%Y %I:%M:%S %p'),int,float,str,str,str,str,float,str,str,int,int,int,int,int,str,float,float,float,float,eval,eval]#,lambda x:Point(float(x))
	print(len(dtypes))
	print(dtypes)
	with open('times.csv') as csv_file:
	    csv_reader = csv.reader(csv_file, delimiter=',')
	    line_count = 0
	    for row in csv_reader:
	        if line_count == 0:
	            #print(f'Column names are {", ".join(row)}')
	            print(len(row))
	            print(row)
	            line_count += 1
	            col_names=row
	            df=df_empty(row,dtypes)
	            print(df.head())

	        else:
	            #print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
	            #print(row)
	            jj=0
	            for (k,r,f) in zip(col_names,row,dtypesf):
	            	print(k)
	            	print(k,f(r))
	            	jj+=1
	            	print(jj)
	            df=df.append({k:f(r) for (k,r,f) in zip(col_names,row,dtypesf)},ignore_index=True)
	            line_count += 1
	        if line_count==k:
	        	print(df.head(k))
	        	break
	    print(f'Processed {line_count} lines.')

import_k_f(10)

