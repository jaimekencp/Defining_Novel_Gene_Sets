import pandas as pd
import mygene

"""
In this script I'm changing the row and column ENSEMBL gene IDS for the human gene ID names.
I also only include genes that have a minimum of 7000 connections. This is because
the full gene network map is too big and leads to memory errors. As such, I create a smaller
matrix that can be inputed into GEPHI as a network matrix and the network can be viewed there
"""

#read in the gene network matrix
gene_network_matrix = pd.read_csv("data_files/network_matrices/gene_network_matrix.csv", index_col=0)

#Get the gene ID for column and rows
ensembl_gene_IDs = [gene for gene in gene_network_matrix.columns]

mg = mygene.MyGeneInfo()
geneSyms= mg.querymany(ensembl_gene_IDs, scopes='ensembl.gene', fields='symbol', species='human')

gene_ID_dic = {}

for gene in ensembl_gene_IDs:
    for query_search in geneSyms:
        if "symbol" in query_search:
            if gene == query_search["query"]:
                gene_ID_dic[gene] = query_search["symbol"]
        elif "notfound" in query_search:
            if gene == query_search["query"]:
                gene_ID_dic[gene] = gene

#change the ENSEMBL row ID and column ID's for their respective gene IDs
gene_network_matrix = gene_network_matrix.rename(columns=gene_ID_dic)
gene_network_matrix = gene_network_matrix.rename(index=gene_ID_dic)

# Now I only include genes that have more than 7000 connections
gene_network_matrix['connection_sum'] = gene_network_matrix.sum(axis=1)

#get average weighted gene node degree
average_weighted_degree = sum(gene_network_matrix['connection_sum']) / len(gene_network_matrix.index)

#for visualization purposes we are going to subset the genes that have over 7000 connections
gene_network_matrix = gene_network_matrix[gene_network_matrix["connection_sum"] > 7000]
gene_network_matrix_columns = list(gene_network_matrix.columns)
gene_network_matrix_rows = list(gene_network_matrix.index)
columns_to_removed = list(set(gene_network_matrix_columns) - set(gene_network_matrix_rows))
gene_network_matrix = gene_network_matrix.drop(columns_to_removed, axis=1)

#Matrix will be used for Gephi
gene_network_matrix.to_csv("data_files/network_data/network_matrices/minimised_gene_network_matrix.csv")



