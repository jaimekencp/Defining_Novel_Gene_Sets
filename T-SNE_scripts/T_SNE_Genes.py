import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

""" This script performs t-SNE on the gene network matrix"""

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


def plot_tsne(tsne_df):
    fig, ax = plt.subplots(figsize=(80, 80))

    ax.scatter(tsne_df['tsne_1'], tsne_df['tsne_2'])

    plt.title('Gene Clusters')
    plt.xlabel('T-sne component 1')
    plt.ylabel('T-sne component 2')
    plt.show()
    fig.savefig("T-SNE/t-sne_results/tsne_results_for_gene_network_matrix/gene_clusters_50.png")
    plt.close(fig)

#read in gene network matrix
gene_network_matrix = pd.read_csv(r"data_files/network_data/network_matrices/gene_network_matrix.csv")
gene_network_matrix, gene_names = read_data(gene_network_matrix)
tsne_df, tsne_component_1, tsne_component_2 = perform_tsne(gene_network_matrix, gene_names)
plot_tsne(tsne_df)

#read in gene network matrix again under different variable name to put the t-sne coordinates
gene_network_matrix_tsne_coordinates = pd.read_csv(r"data_files/network_data/network_matrices/gene_network_matrix.csv")
gene_network_matrix_tsne_coordinates, gene_names = read_data(gene_network_matrix_tsne_coordinates)
gene_network_matrix_tsne_coordinates['tsne_1'] = tsne_component_1
gene_network_matrix_tsne_coordinates['tsne_2'] = tsne_component_2
gene_network_matrix_tsne_coordinates.\
    to_csv(r"T-SNE/t-sne_results/tsne_results_for_gene_network_matrix/GNM_tsne_results.csv")


gnm_tsne = pd.read_csv("T-SNE/t-sne_results/tsne_results_for_gene_network_matrix/GNM_tsne_results.csv")
plot_tsne(gnm_tsne)
