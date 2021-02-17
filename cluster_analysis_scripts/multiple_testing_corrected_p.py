import pandas as pd
import glob
from statsmodels.sandbox.stats.multicomp import multipletests

""" Because of the GO enrichment analysis does not correct for multiple testing, this is 
done here. P values are adjusted usng bonferonni correction followed by multiple clustering
correction."""

path = "C://Users//jaime//PycharmProjects//Defining_Novel_Gene_Sets//cluster_analysis//Cluster_Go_analysis//"
all_files = glob.glob(path + "*.csv")

GO_over_analysis = []
GO_under_analysis = []

for file in all_files:
    if "over" in file:
        print(file)
        df = pd.read_csv(file)
        GO_over_analysis.append(df)
    elif "under" in file:
        print(file)
        df = pd.read_csv(file)
        GO_under_analysis.append(df)

for i, data in enumerate(GO_over_analysis):
    print(i)
    p_adjusted = multipletests(data["Pvalue"], method='bonferroni')
    overall_p_adjusted = p_adjusted[1] * 23
    data["adjusted_Pvalue"] = overall_p_adjusted
    data = data[data['adjusted_Pvalue'] < 0.05]
    data.to_csv("cluster_analysis/Cluster_Go_analysis/Go_over_cluster_" + str(i) + ".csv")

for i, data in enumerate(GO_under_analysis):
    print(i)
    p_adjusted = multipletests(data["Pvalue"], method='bonferroni')
    overall_p_adjusted = p_adjusted[1] * 23
    data["adjusted_Pvalue"] = overall_p_adjusted
    data = data[data['adjusted_Pvalue'] < 0.05]
    data.to_csv("cluster_analysis/Cluster_Go_analysis/Go_under_cluster_" + str(i) + ".csv")



