import pandas as pd
import networkx as nx
from collections import Counter
from create_network_map_scripts.pyvis_for_networkx import draw_graph3
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt

"""
Here the trait network map is built using the connected traits dictionary
where traits (keys) contain the traits they are connected (values) to by 
shared associated gene. The color of the nodes represent the 
categories that the trait belongs to

"""

trait_number_of_associated_genes = pickle.load(
    open('data_files/trait_gene_dictionaries/trait_number_of_associated_genes.pickle', 'rb'))
connected_traits = pickle.load(open('data_files/trait_gene_dictionaries/connected_traits.pickle', 'rb'))

# I make a categorical trait dictionary so that the node colors are labeled accroding to the trait
# category they fall under
categorical_data = pd.read_csv(r"Data_Input_Files/gwasATLAS_v20191115.txt", sep='\t')
categorical_data = categorical_data[['Domain', 'uniqTrait']]
categorical_data = categorical_data.set_index("Domain")
categorical_data = categorical_data.drop_duplicates(keep='first')
categorical_dict = dict(zip(categorical_data.uniqTrait, categorical_data.index))

# In this section I create the trait network map

G = nx.Graph()

for trait, trait_list in connected_traits.items():
    if categorical_dict[trait] == 'Cognitive':
        G.add_node(trait, color='red', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Cardiovascular':
        G.add_node(trait, color='blue', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Infection':
        G.add_node(trait, color='purple', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Metabolic':
        G.add_node(trait, color='yellow', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Ear, Nose, Throat':
        G.add_node(trait, color='brown', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Dermatological':
        G.add_node(trait, color='magenta', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Respiratory':
        G.add_node(trait, color='white', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Reproduction':
        G.add_node(trait, color='black', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Neoplasms':
        G.add_node(trait, color='green', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Mortality':
        G.add_node(trait, color='cyan', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Activities':
        G.add_node(trait, color='grey', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Skeletal':
        G.add_node(trait, color='teal', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Social Interactions':
        G.add_node(trait, color='greenyellow', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Immunological':
        G.add_node(trait, color='darkblue', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Nutritional':
        G.add_node(trait, color='lime', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Neurological':
        G.add_node(trait, color='pink', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Muscular':
        G.add_node(trait, color='darkred', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Psychiatric':
        G.add_node(trait, color='gold', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Endocrine':
        G.add_node(trait, color='turquoise', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Ophthalmological':
        G.add_node(trait, color='chocolate', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Environment':
        G.add_node(trait, color='sandybrown', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Body Structures':
        G.add_node(trait, color='orange', size=trait_number_of_associated_genes[trait])
    elif categorical_dict[trait] == 'Gastrointestinal':
        G.add_node(trait, color='darkgreen', size=trait_number_of_associated_genes[trait])

    for connected_trait in trait_list:
        C = Counter(trait_list)
        G.add_edge(trait, connected_trait, width=C[connected_trait])

draw_graph3(G, output_filename='data_files/network_data/network_map/trait_network_map.html', notebook=False)

print("# of edges: {}".format(G.number_of_edges()))
print("# of nodes: {}".format(G.number_of_nodes()))

print(nx.info(G))

# MAKE FIGURE LEGEND - this will go along with my trait network map figure
# create color palette
trait_colour_dict = {'Cognitive': 'red', 'Cardiovascular': 'blue', 'Infection': 'purple', 'Metabolic': 'yellow',
                     'Ear, Nose, Throat': 'brown', 'Dermatological': 'magenta', 'Respiratory': 'white',
                     'Reproduction': 'black', 'Neoplasms': 'green',
                     'Mortality': 'cyan', 'Activities': 'grey', 'Skeletal': 'teal',
                     'Social Interactions': 'greenyellow',
                     'Immunological': 'darkblue', 'Nutritional': 'lime', 'Neurological': 'pink', 'Muscular': 'darkred',
                     'Psychiatric': 'gold', 'Endocrine': 'turquoise', 'Ophthalmological': 'chocolate',
                     'Environment': 'sandybrown',
                     'Body Structures': 'orange', 'Gastrointestinal': 'darkgreen'}

# Create legend handles manually
handles = [mpl.patches.Patch(color=trait_colour_dict[x], label=x) for x in trait_colour_dict.keys()]
# Create legend
plt.legend(handles=handles)
# Get current axes object and turn off axis
plt.gca().set_axis_off()
plt.show()

for key in trait_colour_dict:
    curLabel = key
    c = next(color)
    for item in dictMy[key]:
        x = item[0]
        y = item[1]
        plt.scatter(x, y, c=c)
    plt.legend(str(curLabel))
plt.show()

""" In this section I perform the network analysis"""

# NETWORK ANALYSIS

# I want to examine the relationship between:
# - Node Size and Degree
# - Node Size and Degree Centrality

# Lets begin with node size and degree
node_size = nx.get_node_attributes(G, "size")


# Because the degree list I get from network x comes in the following
# tuple format [(gene, number of edges), (gene_1, number of edges_1)...]
# I used the bellow function to unpack this and stored it in a dictionary format
def Convert(tup, di):
    for a, b in tup:
        di.setdefault(a, []).append(b)
    return di


# create empty dictionary
degrees_dict = {}

# store node degree in dictionary
Convert(G.degree(), degrees_dict)

# In the next few lines - I check that node size and degree size are of the same length
node_size_keys = []
for key, val in node_size.items():
    node_size_keys.append(key)

degrees_dict_keys = []
for key, val in degrees_dict.items():
    degrees_dict_keys.append(key)

# apparently one trait was missing from node size
# I believe removing one sample will not affect my results
# 1 trait is missing from the node-size dict
list(set(degrees_dict_keys) - set(node_size_keys))
del degrees_dict["Rheumatoid Arthritis"]

node_size_keys = []
for key, val in node_size.items():
    node_size_keys.append(key)

degrees_dict_keys = []
for key, val in degrees_dict.items():
    degrees_dict_keys.append(key)

node_size_list = [val for val in node_size.values()]

degrees_list = []
for node, degree_list in degrees_dict.items():
    for degree in degree_list:
        degrees_list.append(degree)

# Here I make the dataframe that contains node size and degree to be plotted in R
degree_node_size_pd = pd.DataFrame({'Node Size': node_size_list, 'Degrees': degrees_list})
degree_node_size_pd.to_csv("data_files/network_data/network_analysis/degree_node_size.csv", index=False)

# Now I want to look at the relationship between node size and degree centrality
# obtain degree_centrality scores and node size of each trait

degree_centrality = nx.degree_centrality(G)

del degree_centrality["Rheumatoid Arthritis"]

degree_centrality_list = []
for key, val in degree_centrality.items():
    degree_centrality_list.append(val)

degree_centrality_node_size = pd.DataFrame.from_dict(degree_centrality, orient='index')
degree_centrality_node_size["nose_size"] = node_size_list
degree_centrality_node_size.to_csv("data_files/network_data/network_analysis/degree_centrality_node_size.csv")
