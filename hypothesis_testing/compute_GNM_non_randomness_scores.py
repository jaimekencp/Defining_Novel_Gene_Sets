import glob
import pandas as pd
import networkx as nx

"""
Here I am computing the non randomness scores for each of the 100
randomized network matrices as well as the original network matrix.
I then output the null distribution of the non randomness scores 
into a csv file. I do the same for my original network matrix 
"""

def read_data_files(file_path):
    path = file_path
    all_files = glob.glob(path + "/*.csv")

    randomized_GNMs = []
    for matrix in all_files:
        df = pd.read_csv(matrix, index_col=0)
        randomized_GNMs.append(df)

    return randomized_GNMs


def get_random_distribution_array(matrix):
    non_random_scores = []
    for random_GNM in matrix:
        randomized_G = nx.from_pandas_adjacency(random_GNM)
        randomized_random_score_G = nx.non_randomness(randomized_G)
        non_random_scores.append(randomized_random_score_G)

    non_random_score = []
    for i, (score_1, score_2) in enumerate(non_random_scores):
        non_random_score.append([i + 1, score_1])

    return non_random_score


def get_non_random_score(df):
    G = nx.from_pandas_adjacency(df)
    random_score_G = nx.non_randomness(G)
    (score_1, score_2) = random_score_G
    GNM_random_score = [[1, score_1]]

    return GNM_random_score


path = "data_files/randomized_GNM"

randomized_gn_matrices = read_data_files(path)
non_random_score = get_random_distribution_array(randomized_gn_matrices)

null_distribution_df = pd.DataFrame(non_random_score, columns=["randomized_GNM", "R(G)"])
null_distribution_df.to_csv("hypothesis_testing/non_randomness_distribution.csv")

# The two data files are going to be plotted in R and the hypothesis testing is performed there

gnm = pd.read_csv("data_files/network_data/network_matrices/gene_network_matrix.csv", index_col=0)
GNM_random_score = get_non_random_score(gnm)

gnm_randomness_score_df = pd.DataFrame(GNM_random_score, columns=["GNM", "R(G)"])
gnm_randomness_score_df.to_csv("hypothesis_testing/non_randomness_gnm_score.csv")




