import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE


def read_data(trait_gene_matrix, categorical_df):
    trait_gene_matrix.rename(columns={list(trait_gene_matrix)[0]: 'TRAITS'}, inplace=True)
    trait_gene_matrix = trait_gene_matrix.set_index('TRAITS')

    categorical_df = categorical_df[['Domain', 'uniqTrait']]
    categorical_df = categorical_df.set_index("Domain")
    categorical_df = categorical_df.drop_duplicates(keep='first')
    categorical_dict = dict(zip(categorical_df.uniqTrait, categorical_df.index))

    trait_gene_matrix["Domain"] = pd.Series(categorical_dict)
    category = trait_gene_matrix["Domain"]
    category = category.values.flatten()
    trait_gene_matrix = trait_gene_matrix.drop('Domain', axis=1)

    row_name_index = trait_gene_matrix.index
    row_name_index = row_name_index.to_flat_index()

    """ In this part of the function, we are going to insert the domain 
    column again, but this time we are converting the domain name to a
    number that corresponds to the trait domain"""

    trait_gene_matrix["Domain"] = pd.Series(categorical_dict)
    categorical_domains = sorted(list(set(trait_gene_matrix["Domain"])))
    categorical_domains_dict = {}
    for i, domain in enumerate(categorical_domains):
        categorical_domains_dict[domain] = i

    trait_gene_matrix = trait_gene_matrix.replace({"Domain": categorical_domains_dict})
    category_number = trait_gene_matrix["Domain"]
    category_number = category_number.values.flatten()
    category_number = [str(i) for i in category_number]
    trait_gene_matrix = trait_gene_matrix.drop('Domain', axis=1)

    return trait_gene_matrix, row_name_index, category, category_number, categorical_domains_dict


def perform_tsne(df, row_name_index, trait_category, trait_category_number):
    tsne = TSNE(n_components=2, verbose=1, perplexity=50, n_iter=5000)
    tsne_results = tsne.fit_transform(df)

    tsne_1 = tsne_results[:, 0]
    tsne_2 = tsne_results[:, 1]

    tsne_results_df = pd.DataFrame(tsne_results, columns=['tsne_x', 'tsne_y'], index=row_name_index).assign(
        category=trait_category).groupby('category')
    tsne_results_df_numbered = pd.DataFrame(tsne_results, columns=['tsne_x', 'tsne_y'], index=row_name_index).assign(
        category=trait_category_number).groupby('category')

    return tsne_results_df, tsne_results_df_numbered, tsne_1, tsne_2


def plot_tsne_results(tsne_results_df, tsne_results_df_numbered):
    fig, ax = plt.subplots(figsize=(80, 80))

    for name, points in tsne_results_df:
        ax.scatter(points.tsne_x, points.tsne_y, label=name)
        for i, trait_name in enumerate(points.index):
            plt.annotate(trait_name, (points.tsne_x[i], points.tsne_y[i]))

    ax.legend()
    plt.title('Trait Cluster')
    plt.xlabel('T-SNE component 1')
    plt.ylabel('T-SNE component 2')
    plt.show()
    fig.savefig(r"T-SNE/t-sne_results/tsne_results_for_trait_network_matrix/trait_cluster_labeled.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(80, 80))

    for name, points in tsne_results_df_numbered:
        ax.scatter(points.tsne_x, points.tsne_y, label=name)
        for i, trait_name in enumerate(points.index):
            plt.annotate(name, (points.tsne_x[i], points.tsne_y[i]))

    ax.legend()
    plt.title('Trait Cluster')
    plt.xlabel('T-SNE component 1')
    plt.ylabel('T-SNE component 2')

    plt.show()
    fig.savefig(r"T-SNE/t-sne_results/tsne_results_for_trait_network_matrix/trait_cluster_point_labeled.png")
    plt.close(fig)


trait_network_matrix = pd.read_csv(r"data_files\network_data\network_matrices\trait_network_matrix")
categorical_data = pd.read_csv(r"Data_Input_Files\gwasATLAS_v20191115.txt", sep="\t")
trait_network_matrix, gene_name, domain_category, domain_number, categorical_dictionary = \
    read_data(trait_network_matrix, categorical_data)

pd.DataFrame(list(categorical_dictionary.items()), columns=['Trait', 'Index']).\
    to_csv("T-SNE/t-sne_results/tsne_results_for_trait_network_matrix/legend.csv")

tsne_result, tsne_results_numbered, tsne_component_1, tsne_component_2 = \
    perform_tsne(trait_network_matrix, row_name_index=gene_name, trait_category=domain_category,
                 trait_category_number=domain_number)
plot_tsne_results(tsne_results_df=tsne_result, tsne_results_df_numbered=tsne_results_numbered)

# I read in the trait network matrix under different variable name to add t-sne coordinates
trait_network_matrix_tsne_coordinates = pd.read_csv(r"data_files\network_data\network_matrices\trait_network_matrix")
trait_network_matrix_tsne_coordinates.rename(columns={list(trait_network_matrix_tsne_coordinates)[0]: 'TRAITS'},
                                             inplace=True)
trait_network_matrix_tsne_coordinates = trait_network_matrix_tsne_coordinates.set_index('TRAITS')
trait_network_matrix_tsne_coordinates['tsne_1'] = tsne_component_1
trait_network_matrix_tsne_coordinates['tsne_2'] = tsne_component_2
trait_network_matrix_tsne_coordinates.to_csv(
    "T-SNE/t-sne_results/tsne_results_for_trait_network_matrix/trait_network_matrix_tsne_results.csv")
