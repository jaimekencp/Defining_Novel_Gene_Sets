import pandas as pd
import networkx as nx
from create_network_map_scripts.pyvis_for_networkx import draw_graph3
from collections import Counter

""" 
This script creates the full gene network matrix and outputs it as an html file
WARNING: Only computers with a 16 gb of RAM will be able to view the gene network map
"""

# Dictionary gives the me the list of traits each gene is connected to
gene_significant_traits = pd.read_pickle("data_files/trait_gene_dictionaries/gene_significant_traits.pickle")
# Gives the gene - gene connections
connected_genes = pd.read_pickle("data_files/trait_gene_dictionaries/connected_genes.pickle")

def gene_network(connected_genes, gene_significant_traits):

    for gene, gene_list in connected_genes.items():
        G.add_node(gene, size=len(gene_significant_traits[gene]))
        for connected_gene in gene_list:
            C = Counter(gene_list)
            G.add_edge(gene, connected_gene, width=C[connected_gene])

    draw_graph3(G, output_filename='data_files/network_data/network_map/gene_network_map.html', notebook=False)

    print("# of edges: {}".format(G.number_of_edges()))
    print("# of nodes: {}".format(G.number_of_nodes()))


def combined_dicts(d1, d2):
    ds = [d1, d2]
    d = {}
    for k in d1.keys():
        d[k] = tuple(d[k] for d in ds)

    return d


# This section calculates the degree centrality of the nodes in order
# to find which genes act as central point clusters

#To analyse the network without actually having to output it, I'll be using the
#networkX library
gene_network_matrix = pd.read_csv("data_files/network_matrices/gene_network_matrix.csv", index_col=0)
#Here I create a networkx object that contains my gene network map
gene_network_map = nx.from_pandas_adjacency(gene_network_matrix)

#gives me the number of nodes/edges/average degree
print(nx.info(gene_network_map))
print(nx.density(gene_network_map))

# Get gene centrality score dict
degree_centrality = nx.degree_centrality(gene_network_map)

# create node size dictionary
node_size_dict = {}
for gene, trait_list in gene_significant_traits.items():
    node_size_dict[gene] = len(trait_list)

# combine both dictionaries
degree_centrality_node_size_dict = combined_dicts(degree_centrality, node_size_dict)

#output into a dataframe which will be used in R to plot the relationship between node size and degree centrality
degree_centrality_pd = pd.DataFrame.from_dict(degree_centrality_node_size_dict, orient='index', columns=["Degree_Centrality", "Node_Size"])
degree_centrality_pd.to_csv("data_files/network_data/network_analysis/gene_node_degree_centrality.csv")
