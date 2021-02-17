import pandas as pd
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import glob

""" Here I perform t-SNE on each of the 100 permuted gene network matrices,
followed by DBSCAN to output the number of clusters for each of the 100 tsne plots"""

path_second = r"C:/Users/jaime/PycharmProjects/master_project/randomized_GNM"
all_files_second = glob.glob(path_second + "/*.csv")

randomized_matrices = []

for matrix in all_files_second:
    df = pd.read_csv(matrix)
    randomized_matrices.append(df)

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

    plt.title('Gene Clusters - Randomized')
    plt.xlabel('T-sne component 1')
    plt.ylabel('T-sne component 2')
    fig.savefig(
        r"T_SNE/randomized_TSNE_plots/gene_clusters_randomized_" + str(0+index+1) + ".png")
    plt.close(fig)

def perform_DBSCAN(gene_network_matrix_tsne):

    GNM_tsne_cols = gene_network_matrix_tsne[["tsne_1", "tsne_2"]]
    model = DBSCAN(eps=6, min_samples=50).fit(GNM_tsne_cols)
    colors = model.labels_

    number_of_clusters = len(set(colors))

    return number_of_clusters, colors

for i, randomized_gnm in enumerate(randomized_matrices):
    randomized_gnm_df, gene_ids = read_data(randomized_gnm)
    randomized_gnm_df_tsne, tsne_component_1, tsne_component_2 = perform_tsne(randomized_gnm_df, gene_ids)
    plot_tsne(randomized_gnm_df_tsne, i)
    randomized_gnm_df["tsne_1"] = tsne_component_1
    randomized_gnm_df["tsne_2"] = tsne_component_2
    cluster_number, colors = perform_DBSCAN(randomized_gnm_df)
    randomized_gnm_df["DBSCAN cluster"] = colors
    randomized_gnm_df.to_csv("T_SNE/randomized_T_SNE_DBSCAN_results/gene_network_matrix_tsne_dbscan_" + str(0+i+1) + ".csv")




