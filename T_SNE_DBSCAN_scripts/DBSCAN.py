import pandas as pd
from sklearn.cluster import DBSCAN

""" Here I perform DBSCAN on the tsne gene network matrix"""

# perform DBSCAN - try to just get the number of clusters as an output
gene_network_matrix_tsne = pd.read_csv(r"T-SNE/t-sne_results/tsne_results_for_gene_network_matrix/GNM_tsne_results.csv", index_col=0)
GNM_tsne_cols = gene_network_matrix_tsne[["tsne_1", "tsne_2"]]
model = DBSCAN(eps=6, min_samples=50).fit(GNM_tsne_cols)
colors = model.labels_
gene_network_matrix_tsne["DBSCAN_cluster"] = colors

# make a new columns which specifies the colour of each gene based on the cluster
# it belongs to. The keys are the clusters and values are ggplot color codes

cluster_color_dict = {-1: '#000000', 0: '#FFCC00', 1: '#FF6600', 2: '#FF9966',
                      3: '#CC3300', 4: '#99CC00', 5: '#CCFF00', 6: '#333300',
                      7: '#666600', 8: '#FFFF00', 9: '#CC9933', 10: '#996600',
                      11: "#990000", 12: "#FF0000", 13: "#FF6699", 14: "#FF0099",
                      15: "#CC6699", 16: "#660066", 17: "#00FFCC", 18: "#00FFFF",
                      19: "#009999", 20: "#3366CC", 21: "#0000FF", 22: "#999999"}

gene_network_matrix_tsne['color_codes'] = gene_network_matrix_tsne['DBSCAN_cluster'].map(cluster_color_dict)

gene_network_matrix_tsne.to_csv("data_files/gene_network_matrix_tsne_dbscan.csv")
