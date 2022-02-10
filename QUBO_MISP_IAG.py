from dwave_qbsolv import QBSolv
from dimod import BinaryQuadraticModel
from dwave.system import EmbeddingComposite, DWaveSampler
from dwave.system import LeapHybridSampler
from itertools import combinations
import pandas as pd
import networkx as nx
import time
import math

# Define the sampler that will be used to run the problem
#sampler = EmbeddingComposite(DWaveSampler())
#sampler = LeapHybridSampler()
sampler = QBSolv()
 
lines = pd.read_csv('lines52.csv')
nodes = pd.read_csv('nodes52_3products.csv')
# =============================================================================
# lines = pd.read_csv('lines52.csv')
# nodes = pd.read_csv('nodes52_3products.csv')
# =============================================================================
options = pd.read_csv('options_3products.csv')
forbidden = pd.read_csv('forbiddenPairs_3products.csv')

allOptions = []
for option in options['Option']:
    if option != 'discontinuity':
        allOptions += [int(option)]

iaG = nx.Graph()
for u in nodes['Node']:
    comb = combinations(allOptions, 2)
    for (p,q) in list(comb):
        iaG.add_edge((u,p),(u,q))

for line in lines['Line']:
    [source_line] = lines.loc[lines['Line']==line,'Source']
    [target_line] = lines.loc[lines['Line']==line,'Target']
    for pair in forbidden['Pair']:
        [source_pair] = forbidden.loc[forbidden['Pair']==pair,'Source']
        [target_pair] = forbidden.loc[forbidden['Pair']==pair,'Target']
        iaG.add_edge((source_line,source_pair),(target_line,target_pair))
        iaG.add_edge((source_line,target_pair),(target_line,source_pair))

print('# of nodes of IAG = ',len(iaG.nodes()))
print('# of edges of IAG = ',len(iaG.edges()))

bigM = 0.0
for u in nodes['Node']:
    ubRev_u = 0.0
    for p in allOptions:
        [product_p] = options.loc[options['Option']==str(p),'Product']
        [price_p] = options.loc[options['Option']==str(p),'Price']
        [value_u_p] = nodes.loc[nodes['Node']==u,'Value%s'%product_p]
        if value_u_p < float(price_p):
            iaG.remove_node((u,p))
        if value_u_p >= float(price_p):
            if ubRev_u < float(price_p):
                ubRev_u = float(price_p)
    bigM += ubRev_u
    
bigM += 1.0

Q1 = {}
Q2 = {}

for (u,p) in iaG.nodes():
    [price_p] = options.loc[options['Option']==str(p),'Price']
    Q1['X[%s,%s]'%(u,p)] = -float(price_p)
    
for ((u,p),(v,q)) in iaG.edges():
    Q2['X[%s,%s]'%(u,p),'X[%s,%s]'%(v,q)] = bigM
    
    
bqm = BinaryQuadraticModel(Q1, Q2, 0.0,'BINARY')
sampleset = sampler.sample(bqm,num_reads=100)
samplesetPandas = sampleset.to_pandas_dataframe(sample_column=True)   
print('near optimal =',-min(sampleset.data_vectors['energy']))
#print(sampleset.lowest(atol=.5))
print('check if same as above =',-samplesetPandas.loc[[0],'energy'])
isFeasible = 1
for ((u,p),(v,q)) in iaG.edges():
    xVal_u_p = samplesetPandas.loc[[0],'sample'][0]['X[%s,%s]'%(u,p)]
    xVal_v_q = samplesetPandas.loc[[0],'sample'][0]['X[%s,%s]'%(v,q)]
    if xVal_u_p * xVal_v_q > 0:
        isFeasible = 0
print('Feasibility = ',isFeasible) # isFeasibility = 0 means the solution is infeasible

            







