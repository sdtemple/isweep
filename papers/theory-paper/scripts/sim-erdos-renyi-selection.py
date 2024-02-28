import sys
from isweep import *
import networkx as nx

fileout = sys.argv[1] # output file prefix
sample_size = sys.argv[2] # diploid sample size
Ne_file = sys.argv[3]  # file with Ne dictionary values
ibd_length = sys.argv[4] # IBD length cutoff
num_sim = sys.argv[5] # number of simulations (these are for histogram of each replicate)
num_rep = sys.argv[6] # number of replicates (each is an experiment)

num_sim = int(float(num_sim))
num_rep = int(float(num_rep))
sample_size = int(float(sample_size))
sample_size2 = 2 * sample_size # erdos renyi specific
ibd_length = float(ibd_length)

s = float(sys.argv[7])
p = float(sys.argv[8])

print('number of replicates')
print(str(num_rep))
print('\n')
print('number of simulations')
print(str(num_sim))
print('\n')
print('sample size')
print(str(sample_size))
print('\n')
print('ibd length cutoff')
print(str(ibd_length))
print('\n')
print('simulation iterator')

def is_tree(the_subgraph):
    '''Check if subgraph is a tree'''
    num_nodes = the_subgraph.number_of_nodes()
    num_edges = the_subgraph.number_of_edges()
    if num_nodes == (num_edges + 1):
        return True
    return False

def is_complete(the_subgraph):
    '''Check if subgraph is complete'''
    num_nodes = the_subgraph.number_of_nodes()
    num_edges = the_subgraph.number_of_edges()
    complete_edges = num_edges * (num_edges - 1) / 2
    if num_edges >= complete_edges:
        return True
    return False

Ne = read_Ne(Ne_file)
prop = probability_ibd_isweep(s,p,Ne,long_ibd=ibd_length)

fileout=fileout+'-'+'selcoef'+str(s)+'-freq'+str(p)+'-'+Ne_file+'-size'+str(sample_size)+'-cM'+str(ibd_length)+'-sim'+str(num_sim)+'-rep'+str(num_rep)+'.txt'
with open(fileout,'w') as f:
	f.write('selection coefficient: '+str(s)+'\n')
	f.write('allele frequency: '+str(p)+'\n')
	f.write('description of comma-separated inputs\n')
	f.write('num_tracts,largest_component,tree_order2,tree_order3,tree_order4,tree_order5,complete3,complete4,complete5,iscomplete\n')
	f.write('\n')
	for sim in range(num_sim):
		if sim % 10 == 0:
			print(sim)
		for rep in range(num_rep-1):
			the_graph = nx.erdos_renyi_graph(sample_size2, prop)
			num_tracts = the_graph.number_of_edges()
			f.write(str(num_tracts))
			f.write(',')
			tree2 = 0
			tree3 = 0
			tree4 = 0
			tree5 = 0
			complete3 = 0
			complete4 = 0
			complete5 = 0
			hascycle = 0
			complete3plus = 0
			mx = 0
			for component in nx.connected_components(the_graph):
				the_subgraph = the_graph.subgraph(component)
				degree_component = len(component)
				if degree_component > mx:
					mx = degree_component
				# look into components
				if is_tree(the_subgraph):
					if degree_component == 2:
						tree2 += 1
					if degree_component == 3:
						tree3 += 1
					if degree_component == 4:
						tree4 += 1
					if degree_component == 5:
						tree5 += 1
				# look into complete graphs
				if is_complete(the_subgraph):
					if degree_component >= 3:
						complete3plus += 1
					if degree_component == 3:
						complete3 += 1
					if degree_component == 4:
						complete4 += 1
					if degree_component == 5:
						complete5 += 1
			f.write(str(mx))
			f.write(',')
			f.write(str(tree2))
			f.write(',')
			f.write(str(tree3))
			f.write(',')
			f.write(str(tree4))
			f.write(',')
			f.write(str(tree5))
			f.write(',')
			f.write(str(complete3))
			f.write(',')
			f.write(str(complete4))
			f.write(',')
			f.write(str(complete5))
			f.write(',')
			f.write(str(complete3plus))
			f.write('\t')
		the_graph = nx.erdos_renyi_graph(sample_size2,prop)
		num_tracts = the_graph.number_of_edges()
		f.write(str(num_tracts))
		f.write(',')
		tree2 = 0
		tree3 = 0
		tree4 = 0
		tree5 = 0
		complete3 = 0
		complete4 = 0
		complete5 = 0
		hascycle = 0
		complete3plus = 0
		mx = 0
		for component in nx.connected_components(the_graph):
			the_subgraph = the_graph.subgraph(component)
			degree_component = len(component)
			if degree_component > mx:
				mx = degree_component
			# look into components
			if is_tree(the_subgraph):
				if degree_component == 2:
					tree2 += 1
				if degree_component == 3:
					tree3 += 1
				if degree_component == 4:
					tree4 += 1
				if degree_component == 5:
					tree5 += 1
			# look into complete graphs
			if is_complete(the_subgraph):
				if degree_component >= 3:
					complete3plus += 1
				if degree_component == 3:
					complete3 += 1
				if degree_component == 4:
					complete4 += 1
				if degree_component == 5:
					complete5 += 1
		f.write(str(mx))
		f.write(',')
		f.write(str(tree2))
		f.write(',')
		f.write(str(tree3))
		f.write(',')
		f.write(str(tree4))
		f.write(',')
		f.write(str(tree5))
		f.write(',')
		f.write(str(complete3))
		f.write(',')
		f.write(str(complete4))
		f.write(',')
		f.write(str(complete5))
		f.write(',')
		f.write(str(complete3plus))
		f.write('\n')
