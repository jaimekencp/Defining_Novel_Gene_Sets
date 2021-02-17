import pandas as pd
import pickle

Trait_GC = pd.read_csv(r"Data_Input_Files/Trait_GC.csv")
Trait_GC = Trait_GC.pivot_table(columns='uniqTrait2', index='uniqTrait1', values='rg')
Trait_GC = Trait_GC.fillna(0)

# For any traits that exhibit a GC score > 0.80, they are stored as dictionary keys for
# the trait that they are being compared to
correlated_traits = Trait_GC.apply(lambda x: x.index[x >= 0.80].tolist(), axis=1).to_dict()
correlated_traits_df = pd.DataFrame(list(correlated_traits.items()), columns=['Key-Trait', 'Correlated_Traits'])

# Traits that do not contain traits as dictionary keys mean that they are
# the representative trait for that group of highly correlated traits
# These representative traits are stored in the non_correlated_traits list
non_correlated_traits = []
for trait, trait_correlated_list in correlated_traits.items():
    if not trait_correlated_list:
        non_correlated_traits.append(trait)


df = pd.read_csv(r"Data_Input_Files/17518_genes_401_traits.csv")
df.rename(columns={list(df)[0]: 'TRAITS'}, inplace=True)
df = df.set_index('TRAITS')
# Traits that are not in the non_correlated_traits list are removed,
# That way we are only working with non-correlated traits
df = df[[c for c in df.columns if c in non_correlated_traits]]

# A gene is said to be significantly associated to a trait if their measured
# p value for that trait is smaller than 5E-08
P_VALUE_CUTOFF = 5E-08

# In these next few lines of code I create a dictionary with traits(keys)
# and gene list (key value) storing all genes that are significantly
# associated to the trait of interest
trait_significant_genes = {}
for trait in df.columns:
    significant_genes_list = df.loc[df[trait] < P_VALUE_CUTOFF].index
    for gene in significant_genes_list:
        trait_significant_genes.setdefault(trait, []).append(gene)

#Make a bipartite dataframe - This file will be used as data input for the randomization matlab scripts
trait_significant_genes_list = [(key, x) for key, val in trait_significant_genes.items() for x in val]
trait_significant_genes_DF = pd.DataFrame(trait_significant_genes_list, columns=['Trait', 'Significant_Genes'])
trait_significant_genes_DF.to_csv("data_files/bipartite_trait_gene.csv")


# In the next few lines, I build a dataframe that tells me the number of
# associated genes per trait. It shows that some traits have a considerable
# amount of genes associated to them.
trait_number_of_associated_genes = {}
for trait, gene in trait_significant_genes.items():
    trait_number_of_associated_genes[trait] = len(gene)

with open('data_files/trait_gene_dictionaries/trait_number_of_associated_genes.pickle', 'wb') as handle:
    pickle.dump(trait_number_of_associated_genes, handle, protocol=pickle.HIGHEST_PROTOCOL)


trait_number_of_associated_genes_pd = pd.DataFrame(trait_number_of_associated_genes.items(), columns=['Trait', 'Number of associated genes'])
trait_number_of_associated_genes_pd.to_csv(r"data_files/trait_number_of_associated_genes.csv")


# In the next couple of lines I build my trait-gene matrix
trait_significant_genes_copy = trait_significant_genes.copy()

connected_traits = {}
for trait, gene_list in trait_significant_genes.items():
    for trait_copy, gene_list_copy in trait_significant_genes_copy.items():
        for gene in gene_list:
            if gene in gene_list_copy:
                connected_traits.setdefault(trait, []).append(trait_copy)

# Gets rid of diagonal
connected_traits = {k: [vi for vi in v if k != vi] for k, v in connected_traits.items()}
# Gets rid of traits that do not share any significantly associated genes
connected_traits = {k: v for k, v in connected_traits.items() if v}

# I want to store this dictionary in a pickle as it is going to be used later
# on for making the network maps
with open('data_files/trait_gene_dictionaries/connected_traits.pickle', 'wb') as handle:
    pickle.dump(connected_traits, handle, protocol=pickle.HIGHEST_PROTOCOL)

# These next few lines of code build the matrix
edges = [(a, b) for a, bs in connected_traits.items() for b in bs]
network_matrix = pd.DataFrame(edges)
adj_matrix = pd.crosstab(network_matrix[0], network_matrix[1])
adj_matrix.to_csv(r"data_files/network_data/network_matrices/trait_network_matrix.csv")


