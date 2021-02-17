import glob
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import levene

""" Here I perform a two sample t test on the 
alternative hypothesis distribution and null distribution"""

def read_files(file_path):

    path = file_path
    all_files = glob.glob(path + "/*.csv")

    tsne_results = []

    for matrix in all_files:
        print(matrix)
        df = pd.read_csv(matrix, index_col=0)
        tsne_results.append(df)

    return tsne_results

def array_distribution(results_list, name_of_distribution):

    """name_of_distribution: string"""

    rows = []
    for i, matrix in enumerate(results_list):
        dbscan_clusters = set(matrix["DBSCAN cluster"])
        if -1 in dbscan_clusters:
            dbscan_clusters.remove(-1)
            number_of_clusters = len(dbscan_clusters)
            rows.append([i + 1, number_of_clusters])

    distribution_df = pd.DataFrame(rows, columns=["randomized_GNM", "number_of_clusters"])
    distribution_df.to_csv("T_SNE/tsne_dbscan_distributions/" + name_of_distribution + ".csv")
    distribution = list(distribution_df["number_of_clusters"])

    return distribution

def perform_t_test(distribution_1, distriubtion_2):

    stat, P = levene(distribution_1, distriubtion_2)
    print(stat, P)

    if P > 0.05:
        test = stats.ttest_ind(a=distribution_1, b=distriubtion_2, equal_var=True)
        t_score = test[0]
        p_value = test[1]
        print(t_score, p_value)
        with open("equal_variance.txt", "w") as text_file:
            print("T score: {}".format(t_score), file=text_file)
            print("P value: {}".format(p_value), file=text_file)
    else:
        test = stats.ttest_ind(a=distribution_1, b=distriubtion_2, equal_var=False)
        t_score = test[0]
        p_value = test[1]
        print(t_score, p_value)
        with open("unequal_variance.txt", "w") as text_file:
            print("T score: {}".format(t_score), file=text_file)
            print("P value: {}".format(p_value), file=text_file)

def plot_histogram(x, y):

    y1 = x
    y2 = y
    colors = ['r', 'b']

    # plots the histogram
    fig, ax1 = plt.subplots()
    ax1.hist([y1, y2], color=colors)
    ax1.set_xlim(15, 35)
    ax1.set_ylabel("Count")
    ax1.set_xlabel("Number of Clusters")
    plt.tight_layout()
    plt.legend(["Permuted GNM", "Non Permuted GNM"])
    plt.show()

def plot_boxplot(x, y):

    sns.boxplot(data=[x, y])
    box_plot = sns.boxplot(
        data=[x, y],
        palette=[sns.xkcd_rgb["red"], sns.xkcd_rgb["blue"]],
        showmeans=True,
    )
    box_plot.set(ylim=(20, 35))

path_1 = r"C:/Users/jaime/PycharmProjects/Defining_Novel_Gene_Sets/T_SNE/T_SNE_DBSCAN_result/HA_distribution"
path_2 = r"C:/Users/jaime/PycharmProjects/Defining_Novel_Gene_Sets/T_SNE/randomized_T_SNE_DBSCAN_results"

tsne_results_non_randomized = read_files(path_1)
tsne_results_randomized = read_files(path_2)

tsne_results_non_randomized_distribution_1 = array_distribution(tsne_results_non_randomized, "tsne_dbscan_non_random_distribution")
tsne_results_randomized_distribution_2 = array_distribution(tsne_results_randomized, "tsne_dbscan_random_distribution")

perform_t_test(tsne_results_non_randomized_distribution_1, tsne_results_randomized_distribution_2)
perform_t_test_2(tsne_results_non_randomized_distribution_1, tsne_results_randomized_distribution_2)

permuted_distribution = pd.read_csv("data_files/T_SNE/tsne_dbscan_distributions/tsne_dbscan_random_distribution.csv")
permuted_distribution = permuted_distribution["number_of_clusters"].tolist()

non_permuted_distribution = pd.read_csv(
    "data_files/T_SNE/tsne_dbscan_distributions/tsne_dbscan_non_random_distribution.csv")
non_permuted_distribution = non_permuted_distribution["number_of_clusters"].tolist()

plot_histogram(permuted_distribution, non_permuted_distribution)
plot_boxplot(permuted_distribution, non_permuted_distribution)

perform_t_test(non_permuted_distribution, permuted_distribution)

