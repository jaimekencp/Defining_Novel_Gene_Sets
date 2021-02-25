import pandas as pd
import glob
import itertools
import csv

#Read in T-SNE dbscan
gene_tsne_dbscan_df = pd.read_csv("T_SNE/T_SNE_DBSCAN_result/gene_network_matrix_tsne_dbscan.csv", index_col=0)

#transform genes into their respective gene id name
#here I make a dictionary with the ensembl id as key and their gene symbol as value
#for ensembl id where gene symbol was not found, we just keep their ensembl gene id
gene_id_df = pd.read_csv("FUMA_gene_sets/geneIDs.txt", sep="\t")

ensembl_id = gene_id_df["ensg"]
gene_symbol = gene_id_df["symbol"]
gene_dict = dict(zip(ensembl_id, gene_symbol))

# Now I'm going to change the row and column names in my gene network matrix to gene id symbol:
gene_tsne_dbscan_df = gene_tsne_dbscan_df.rename(columns=gene_dict, index=gene_dict)

# I'm going to read in the gene_set data files from each of my 23 gene clusters and stored in a list
path = r"C:/Users/jaime/PycharmProjects/Defining_Novel_Gene_Sets/FUMA_gene_sets/GS_"
file_names = glob.glob(path + "*.txt")

gene_set_data_files = []

# Here I read in the files and multiply the adjP column by 23 to correct for multiple gene clusters
# I then drop rows/genes whose P value is > 0.05
for file in file_names:
    gene_set_df = pd.read_csv(file, sep="\t")
    print(gene_set_df)
    gene_set_df["corrected_P"] = gene_set_df['adjP'] * 23
    gene_set_df = gene_set_df.drop(gene_set_df[gene_set_df.corrected_P > 0.05].index)
    gene_set_data_files.append(gene_set_df)

print(gene_set_data_files)
# Here I save the new gene set datafiles
for i, gene_set_df in enumerate(gene_set_data_files):
    gene_set_df.to_csv("FUMA_gene_sets/GS_" + str(0+i) + ".csv")

# Now I'm only interested in the genes column of the gene set datasets
gene_sets = []
for gene_set_df in gene_set_data_files:
    gene_set_df = gene_set_df["genes"]
    gene_sets.append(gene_set_df)

print(gene_sets)
# Here I'm inputing all genes from the gene sets into 23 separate lists
for i, gene_set in enumerate(gene_sets):
    locals()["gene_sets_in_cluster_" + str(i)] = gene_set
    locals()["genes_for_each_cluster_" + str(i)] = []
    for gene_rows in locals()["gene_sets_in_cluster_" + str(i)]:
        gene_rows = gene_rows.split(":")
        locals()["genes_for_each_cluster_" + str(i)].append(gene_rows)
        locals()["gene_sets_of_" + str(i)] = list(itertools.chain(*locals()["genes_for_each_cluster_" + str(i)]))


# Now I want to store all the genes from each of my clusters in a dictionary where
# keys are clusters and values are the genes within the clsuters
gene_ids = gene_tsne_dbscan_df.index
cluster_ids = gene_tsne_dbscan_df["DBSCAN_cluster"]

cluster_id_dict = dict(zip(gene_ids, cluster_ids))

gene_cluster_dict = {}
for gene, cluster_id in cluster_id_dict.items():
    gene_cluster_dict.setdefault(cluster_id, []).append(gene)

# In the next few lines I basically see the difference in genes

with open('FUMA_gene_sets/non_overlapping_genes.csv', 'w') as f1:
    writer = csv.writer(f1, lineterminator='\n')
    for i in range(0, 23):
        try:
            set_difference = list(set(gene_cluster_dict[i]) - (set(locals()["gene_sets_of_" + str(i)])))
            writer.writerows([str(i), set_difference])
        except:
            KeyError
            cluster_15 = list(gene_cluster_dict[15])
            print(len(cluster_15))
            writer.writerows([str(i), cluster_15])

