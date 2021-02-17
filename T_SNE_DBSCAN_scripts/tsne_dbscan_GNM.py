import pandas as pd
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt


""" To build my alternate hypothesis distribution here I perform tsne on the gene network matrix
100 times followed by DBSCAN to output the number of clusters for each of the 100 tsne plots"""

def read_data(df):

    df.rename(columns={list(df)[0]: 'GENES'}, inplace=True)
    df = df.set_index('GENES')

    row_name_index = df.index
    row_name_index = row_name_index.to_flat_index()

    return df, row_name_index


def perform_tsne(df, row_name_index):
    tsne = TSNE(n_components=2, verbose=1, perplexity=50, n_iter=10000)
    tsne_results = tsne.fit_transform(df)

    tsne_1 = tsne_results[:, 0]
    tsne_2 = tsne_results[:, 1]

    tsne_results_df = pd.DataFrame(tsne_results, columns=['tsne_x', 'tsne_y'], index=row_name_index)

    return tsne_results_df, tsne_1, tsne_2


def plot_tsne(tsne_df, index):
    fig, ax = plt.subplots(figsize=(80, 80))

    ax.scatter(tsne_df['tsne_x'], tsne_df['tsne_y'])

    plt.title('Gene Clusters')
    plt.xlabel('T-sne component 1')
    plt.ylabel('T-sne component 2')
    fig.savefig(
        r"T_SNE/T_SNE_DBSCAN_result/tsne_plots/gene_clusters_" + str(0+index+1) + ".png")
    plt.close(fig)

def perform_DBSCAN(gene_network_matrix_tsne):

    GNM_tsne_cols = gene_network_matrix_tsne[["tsne_1", "tsne_2"]]
    model = DBSCAN(eps=6, min_samples=50).fit(GNM_tsne_cols)
    colors = model.labels_

    number_of_clusters = len(set(colors))

    return number_of_clusters, colors


gnm_df = pd.read_csv("data_files/network_data/network_matrices/gene_network_matrix.csv")

for i in range(100):
    gene_network_matrix, gene_ids = read_data(gnm_df)
    gene_network_matrix_tsne_df, tsne_component_1, tsne_component_2 = perform_tsne(gene_network_matrix, gene_ids)
    gene_network_matrix["tsne_1"] = tsne_component_1
    gene_network_matrix["tsne_2"] = tsne_component_2
    plot_tsne(gene_network_matrix, i)
    cluster_number, colors = perform_DBSCAN(gene_network_matrix)
    gene_network_matrix["DBSCAN cluster"] = colors
    gene_network_matrix.to_csv("T_SNE/T_SNE_DBSCAN_result/HA_distribution/gnm_tsne_dbscan_" + str(0+i+1) + ".csv")
