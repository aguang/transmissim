#! /usr/python
# AUGUST GUANG
# NOVEMBER-DECEMBER 2013
# networkx drawing added JULY 2014

# tests the sensitivity of fablast under the following criteria:
# a clustering is a COMPLETE CLUSTER if:
#    1) all genes in the cluster are from the same tree (i.e. evolved from the same ancestor)
#    2) multiple trees are not included in the cluster
#    3) all genes in the tree are in the cluster
# a clustering is a MULTIPLE CLUSTER if:
#    1) multiple trees are included in the cluster
# a clustering is an INCOMPLETE CLUSTER if:
#    1) all genes in the cluster are from the same tree
#    2) not all genes in the tree are in the cluster

# also draws networkx graph with nodes colored by gene since there's only 10 of them (temporary)

# additional assumptions:
#    -genes evolve independently from each other

import sys, getopt, string
import itertools
import glob
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pprint
import networkx as nx

# creates networkx graph from mcl file
def nxGraph(mclFile):
    G = nx.Graph()
    with open(mclFile, 'r') as fin:
        for line in fin:
            nodes = line.split()
            edges = itertools.combinations(nodes, 2)
            G.add_nodes_from(nodes)
            G.add_edges_from(edges)
    return G

def components(graph, min_nodes=0):
	"""
	Returns an iterator of subgraphs for each connected component of `graph`
	with >= `min_nodes`.
	"""
	for subgraph in nx.connected_component_subgraphs(graph):
		if len(subgraph) >= min_nodes:
			yield subgraph

def print_edges(graph, output=sys.stdout):
	"""
	Prints a tab-separated list of edges in `graph`. For debug purposes.
	"""
	for e in graph.edges_iter(data=True):
		print >>output, "%s\t%s" % (e[0], e[1])

def assign_node_colors(graph, c=10):
	"""
	Assigns colors to nodes based on which gene the node is from. For temporary use with homology simulation fablast data, limited to c=10 colors.
	"""
	colors = ['b','g','r','c','m','y','k','w','0.5','#ffa500']
	node_color = []
	for n in graph.nodes_iter(data=False):
            [tmp, geneID] = n.split('-')
            node_color.append(colors[(int(geneID)+1)%c])
            # +1 in int(geneID) mod is to make congruent with blastp scheme
	return node_color

def process(mclFile, numSpecies, freq):
    """
    Processes mcl cluster file to
    1) aggregate statistics: check if all gene IDs are all equivalent and whether all gene IDs are there
    2) generate networkx graphs with (if applicable) nodes colored by gene for visual understanding of clustering
    note:
    0 = complete cluster
    1 = multiple cluster
    2 = incomplete cluster
    """
    clust0 = 0
    clust1 = 0
    clust2 = 0
    G = nxGraph(mclFile)

    with open(mclFile, 'r') as fin:
        for line in fin:
            genes = line.split()
            geneDict = {}
            count = 0
            for gene in genes:
                count = count + 1
                [speciesID, geneID] = gene.split('-')
                if geneID in geneDict:
                    geneDict[geneID].append(speciesID)
                else:
                    geneDict[geneID] = [speciesID]
            # aggregate statistics
            if len(geneDict) > 1: clust1 = clust1 + 1
            elif count == numSpecies:
                clust0 = clust0 + 1
                freq.append(count)
            else:
                clust2 = clust2 + 1
                freq.append(count)
    return clust0, clust1, clust2, G

# computes percentage of complete clusters as well as ratio of
# complete clusters to multiple/incomplete clusters
def stats(clust0, clust1, clust2, freq):
    total = clust0 + clust1 + clust2
    cc = float(clust0)/total*100
    print "Percantage of complete clusters is: ", cc, "%."
    mc = float(clust1)/total*100
    print "Percentage of multiple clusters is: ", mc, "%."
    ic = float(clust2)/total*100
    print "Percentage of incomplete clusters is: ", ic, "%."
    if (clust1 + clust2) != 0:
        ratio = float(clust0)/(clust1 + clust2)
        print "Ratio of complete clusters to multiple/incomplete clusters is: ", ratio
    else:
        print "All clusters are complete! Yay!"

    # generate histogram
    n, bins, patches = plt.hist(freq, histtype='bar')
    zipped = zip(n, bins)
    print "Frequency of clusters at each level is: ", pprint.pprint(zipped)

    plt.xlabel('Number of Species in Cluster')
    plt.ylabel('Number of Clusters')
    plt.title('Histogram of the Number of Clusters with a Specific Number of Species in them')
    plt.savefig('clustHist.png', bbox_inches='tight')
    plt.clf()
    plt.close()

# main
def main(argv):
    mclFile = ''
    numSpecies = 0
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:n:", ["mclFile=", "numSpecies="])
    except getopt.error, msg:
        print 'homologize-stats.py -n numSpecies -m mclFile'
        sys.exit(2)
    # process arguments
    for opt, arg in opts:
        if opt in ("-n", "--numSpecies"):
            numSpecies = arg
        elif opt in ("-m", "--mcl"):
            mclFile = arg

    freq = []
    clust0, clust1, clust2, G = process(mclFile, numSpecies, freq)
    stats(clust0, clust1, clust2, freq)

#    print_edges(G)
    i = 0
    for subgraph in components(G):
        nc = assign_node_colors(subgraph)
        nx.draw(subgraph, node_color=nc)
        plt.axis("off")
        plt.savefig("cluster%d.png" % i)
        plt.clf()
        i = i+1


if __name__ == "__main__":
    main(sys.argv[1:])
