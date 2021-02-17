import pandas as pd
from collections import Counter

""" This script performs the trait enrichment analysis of each cluster"""

def create_cluster_dict(df):
    for i, cluster in enumerate(set(df["DBSCAN_cluster"])):
        if cluster == -1:
            pass
        for gene in df.index:
            if df["DBSCAN_cluster"][gene] == i:
                cluster_dict.setdefault("cluster_%s" %i, []).append(gene)

    return cluster_dict

def produce_cluster_counter_df(cluster):
    cluster_trait_counter = count_traits_per_cluster(cluster)
    cluster_count_df = get_cluster_trait_table(cluster_trait_counter)

    return cluster_count_df


def count_traits_per_cluster(cluster):

    traits = []
    for gene in cluster:
        for trait in gene_significant_traits[gene]:
            traits.append(trait)

    cluster_trait_counter = Counter(traits)

    return cluster_trait_counter


def get_cluster_trait_table(cluster_counter):

    trait_name = []
    trait_count = []
    for i, trait in enumerate(cluster_counter):
        trait_name.append(trait)
        trait_count.append(cluster_counter[trait])

    cluster_count_dict = {'Trait': trait_name, 'Trait_Count': trait_count}
    cluster_count_df = pd.DataFrame(cluster_count_dict)

    return cluster_count_df


gene_significant_traits = pd.read_pickle("data_files/trait_gene_dictionaries/gene_significant_traits.pickle")

GNM_tsne_DBSCAN_df = pd.read_csv("data_files/T_SNE/T_SNE_DBSCAN_result/gene_network_matrix_tsne_dbscan.csv", index_col=0)

cluster_dictionary = create_cluster_dict(GNM_tsne_DBSCAN_df)

for cluster, genes in cluster_dictionary.items():
    print(cluster)
    df_cluster = produce_cluster_counter_df(genes)
    print(df_cluster)
    df_cluster.to_csv("cluster_analysis/analysis_of_" + str(cluster) + ".csv")


