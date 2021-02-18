# defining_novel_gene_sets - workflow

Here I describe which scripts to run in order.

Preparing the data input files to make the Gene and Trait network matrices:

- The first script to run are found in the R_code directory under the name "create_necessary_dataframes.Rmd". 
This script creates the trait by gene matrix and the trait correlation dataframe.

Creating the gene and trait network matrices:

- Next the scripts "create_non_correlated_gene_trait_matrix.py" creates the gene network matrix.
The "create_non_correlated_trait_gene_matrix.py" creates the trait network matrix.

Creating gene and trait network map:
- to create the gene network map, run the "gene_network_map.py" script
- to create the trait network map, run the "trait_network_map.py" script

Make trait t-SNE network:

- Within the "T-SNE_scripts" directory, the "T_SNE_Traits.py" creates the trait t-SNE network map


Hypothesis testing for Gene network map:
- test whether the topology of the gene network map deviates from random, we first need to
randomize our gene network matrix. This can be done using the  matlab code scripts found in the
"create_network_matrices_scripts". The "randmio_und_bipartite" is the randomization function that we 
call in the "randomization_scrip". This creates 100 randomized gene network matrices

- next we want to compute the non-randomness scores for each of the 100 randomizations matrices. This
script can be found in the "hypothesis testing" directory under the name "compute_GNM_non_randomness_scores". 
You can modify the script to just give you the non-randomness score for the original gene network matrix

- The R script under the name "plots.Rmd" creates the plot showing the null distribution which represents the non-randomness scores
of the 100 randomized matrices and the non-randomness score of the original gene network matrix. It also performes the statistical test
to see whether the non-randomness score of the original gene network matrix is significantly different from the null distribution

Make Gene t-SNE network:
- Within the "T-SNE_scripts" directory, the "T_SNE_Genes.py" creates the trait t-SNE network map
- To run DBSCAN on the t-SNE data, run the "DBSCAN.py" script
- To plot the t-SNE embeddings run the "plot.Rmd" script found in the "R_code" directory

Hypothesis testing for Gene t-SNE network:
- To test how the topology of our gene t-SNE network deviates from random, we first need to build the null distriubtion. This can be done by using
the 100 randomized Gene network matrices created previously and running t-SNE on them followed by DBSCAN to compute the number of clusters. This can be done
running the "t_sne_dbscan_randomized_GNMs.py".
- To build the alternate hypothesis distriubtion run the "t_sne_dbscan_GNM.py" script
- To perform the hypothesis testing, run the "t_sne_dbscan_t_test.py"

Gene t-SNE network cluster analysis:
- To count the number of times each trait is associated to each cluster, run the "analyse_clusters.py" script within the "cluster_analysis_scripts"
directory
- To perform the GO enrichment analysis, run the "Cluster_Analysis.Rmd" script found in "R_code" directory
- To adjust the P value of the GO terms, run the "multiple_testing_correcte_p.py" script


