import glob
import pandas as pd


""" Here I'm creating the 100 randomized gene-gene network matrices.
 I read in the 100 randomized gene-trait matrices outputed by the matlab
 'randomization_script.m' script"""

def check_matrix_similarity(randomized_matrices):
    for i, matrix_i in enumerate(randomized_matrices):
        for j, matrix_j in enumerate(randomized_matrices):
            if i == j:
                continue
            elif matrix_i.equals(matrix_j) == True:
                text_file = open("matrices_are_equal.txt", "w")
                text_file.write("matrix " + str(i + 1) + " equals matrix " + str(j + 1))
                text_file.close()
            else:
                text_file = open("data_files/matrices_are_not_equal.txt", "w")
                text_file.write("matrices are not equal")
                text_file.close()



def make_randomized_gene_network_matrices(randomized_matrices):
    for i, matrix in enumerate(randomized_matrices):

        significant_traits = {}
        for gene in matrix.columns:
            significant_traits_list = matrix.loc[matrix[gene] == 1].index
            for trait in significant_traits_list:
                significant_traits.setdefault(gene, []).append(trait)

        genes_to_remove = []

        for gene, trait_list in significant_traits.items():
            if len(trait_list) < 2:
                genes_to_remove.append(gene)

        for gene in genes_to_remove:
            del significant_traits[gene]

        connected_genes = {}

        significant_traits_copy = significant_traits.copy()

        for gene, trait_list in significant_traits.items():
            for gene_copy, trait_list_copy in significant_traits_copy.items():
                for trait in trait_list:
                    if trait in trait_list_copy:
                        connected_genes.setdefault(gene, []).append(gene_copy)

        # If the gene value is the same as it's key, then the gene is removed (this removes diagonal)
        connected_genes = {k: [vi for vi in v if k != vi] for k, v in connected_genes.items()}
        edges = [(a, b) for a, bs in connected_genes.items() for b in bs]
        network_matrix = pd.DataFrame(edges)
        adj_matrix = pd.crosstab(network_matrix[0], network_matrix[1])

        adj_matrix.to_csv(r"data_files/randomized_GNM/randomized_GNM_0" + str(i+1) + ".csv")


path = 'data_files/randomized_gene_trait_matrices/'
all_files = glob.glob(path + "*.csv")

randomized_matrices = []

for matrix in all_files:
    print(matrix)
    df = pd.read_csv(matrix, index_col=0)
    randomized_matrices.append(df)

check_matrix_similarity(randomized_matrices=randomized_matrices)
make_randomized_gene_network_matrices(randomized_matrices=randomized_matrices)

