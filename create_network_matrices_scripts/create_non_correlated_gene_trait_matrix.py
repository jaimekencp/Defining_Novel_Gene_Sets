import pandas as pd
import pickle

Trait_GC = pd.read_csv(r"Data_Input_Files/Trait_GC.csv")
Trait_GC = Trait_GC.pivot_table(columns='uniqTrait2', index='uniqTrait1', values='rg')
Trait_GC = Trait_GC.fillna(0)

# For any traits that exhibit a GC score > 0.80, they are stored as dictionary keys for
# the trait that they are being compared to
correlated_traits = Trait_GC.apply(lambda x: x.index[x >= 0.80].tolist(), axis=1).to_dict()

# Traits that do not contain traits as dictionary keys mean that they are
# the representative trait for that group of highly correlated traits
# These representative traits are stored in the non_correlated_traits list
non_correlated_traits = []
for trait, trait_correlated_list in correlated_traits.items():
    if not trait_correlated_list:
        non_correlated_traits.append(trait)


df = pd.read_csv(r"Data_Input_Files/17518_genes_401_traits.csv")
df.rename(columns={list(df)[0]: 'GENES'}, inplace=True)
df = df.set_index('GENES')

# Traits that are not in the non_correlated_traits list are removed,
# That way we are only working with non-correlated traits
df = df[[c for c in df.columns if c in non_correlated_traits]]

# A gene is said to be significantly associated to a trait if their measured
# p value for that trait is smaller than 5E-08
P_VALUE_CUTOFF = 5E-08

# A dictionary where Genes are keys and with a trait list as dict key that stores all
# traits that the gene of interest is significantly associated to
significant_traits = df.apply(lambda x: x.index[x < P_VALUE_CUTOFF].tolist(), axis=1).to_dict()
# Genes that are non associated to any trait are removed from the dictionary
significant_traits = {k: v for k, v in significant_traits.items() if v}

# Genes that are only associated to one trait are removed, this is because
# we want to analyse genes across traits. If we were to keep genes that
# are only associated to one trait, then it is very possible that we
# would get gene clusters that are driven by a single trait
genes_to_remove = []
for gene, trait_list in significant_traits.items():
    if len(trait_list) < 2:
        genes_to_remove.append(gene)

for gene in genes_to_remove:
    del significant_traits[gene]


# I stored my gene(key)-trait_list(key value) in a pickle
with open('data_files/trait_gene_dictionaries/gene_significant_traits.pickle', 'wb') as handle:
    pickle.dump(significant_traits, handle, protocol=pickle.HIGHEST_PROTOCOL)

# The next several lines is where I build my gene - trait matrix
connected_genes = {}

significant_traits_copy = significant_traits.copy()

for gene, trait_list in significant_traits.items():
    for gene_copy, trait_list_copy in significant_traits_copy.items():
        for trait in trait_list:
            if trait in trait_list_copy:
                connected_genes.setdefault(gene, []).append(gene_copy)


# Stored connected_genes dict as it is going to be used to build the network map
with open('data_files/trait_gene_dictionaries/connected_genes.pickle', 'wb') as handle:
    pickle.dump(connected_genes, handle, protocol=pickle.HIGHEST_PROTOCOL)


# If the gene value is the same as it's key, then the gene is removed (this removes diagonal)
connected_genes = {k: [vi for vi in v if k != vi] for k, v in connected_genes.items()}
edges = [(a, b) for a, bs in connected_genes.items() for b in bs]
network_matrix = pd.DataFrame(edges)
adj_matrix = pd.crosstab(network_matrix[0], network_matrix[1])


# To csv file
adj_matrix.to_csv(r"data_files/network_data/network_matrices/gene_network_matrix.csv")