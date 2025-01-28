#!/usr/bin/env python3

#Load packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gprofiler import GProfiler

#Read data sets 
A_C = pd.read_table("A_vs_C.deseq2.results.tsv")
A_E = pd.read_table("A_vs_E.deseq2.results.tsv")

print("File read") # Confirm data has been read in

#pre processing data
A_C = A_C.dropna() #Drop rows containing NA values for both datasets
A_E = A_E.dropna()

A_C = A_C.drop_duplicates() #remove genes that may have been recorded multiple times
A_E = A_E.drop_duplicates()


print("Data Cleaned") #confirm cleaning step has been completed

#print both datasets containing only the adjusted p - value and log change
A_C[["padj","log2FoldChange"]]

A_E[["padj","log2FoldChange"]]

#keep rows that have padj < 0.05 and LogFoldChange > 1 since that indicates up regulation of gene
#The double square brackets at end allow for only columns of interest to be shown
Up_A_C = A_C[(A_C["padj"] < 0.05) & (A_C["log2FoldChange"] > 1)][["gene_id", "log2FoldChange", "padj"]]

#Since number of rows = number of upregulated genes in the Up_A_C df len can be used to count number of rows
print("Number of significantly upregulated genes in A vs C" , len(Up_A_C))

#Repeat the above steps however log2FoldChange less than 1 indicates downregulated genes
Dn_A_C = A_C[(A_C["padj"] < 0.05) & (A_C["log2FoldChange"] < 1)][["gene_id", "log2FoldChange", "padj"]]

print("Number of significantly downregulated genes in A vs C", len(Dn_A_C))

#Above steps once again repeated but with the A vs E dataset
Up_A_E = A_E[(A_E["padj"] < 0.05) & (A_E["log2FoldChange"] > 1)][["gene_id", "log2FoldChange", "padj"]]

print("Number of significantly upregulated genes in A vs E", len(Up_A_E))

Dn_A_E = A_E[(A_E["padj"] < 0.05) & (A_E["log2FoldChange"] < 1)][["gene_id", "log2FoldChange", "padj"]]

print("Number of sigificantly downregulated genes in A vs E", len(Dn_A_E))

#Describe function gives a varity of summary stastistics 
#This was repeated for Up_A_C Dn_A_C Up_A_E and Dn_A_E
print("Summary stastistics of Adjusted P_value for A vs C", A_C[["padj"]].describe())
#describe() however does not provide the median so this was done seperatly 
print("Median Adjusted P-value for A vs C = ", A_C["padj"].median())
print("Summary stastistics of Log fold change for A vs C", A_C[["log2FoldChange"]].describe())
print("Medium Log fold change for A vs C = ", A_C[["log2FoldChange"]].median())

print("Summary stastistics of Adjusted P_value for A vs E", A_E[["padj"]].describe())
print("Median Adjusted P-value for A vs E = ", A_E["padj"].median())
print("Summary stastistics of Log fold change for A vs E", A_E[["log2FoldChange"]].describe())
print("Medium Log fold change for A vs E = ", A_E[["log2FoldChange"]].median())

#Scatter plot of Log2FoldChange vs  adjusted p value
#-np.log(x) typical of volcano plot 
#.apply lambda allows for each of the individuals points to be coloured differently depending on corresponding p-value 
plt.scatter(x = A_C["log2FoldChange"], y = (A_C["padj"]).apply(lambda x:-np.log10(x)),s=1,label= "No change")

plt.scatter(x = Up_A_C["log2FoldChange"], y = (Up_A_C["padj"]).apply(lambda x: -np.log10(x)), s = 1,label = "Up regulated", color = "red")

# s = size of points

plt.scatter(x = Dn_A_C["log2FoldChange"], y = (Dn_A_C["padj"]).apply(lambda x: -np.log10(x)), s = 1 ,label = "Down regulated", color = "blue")

plt.xlabel("Log2 Fold Change")
plt.ylabel("- log10 Adjusted P-value")
plt.title("A vs C differentially expressed genes")
plt.legend()
#Save png
plt.savefig("Scatter_AvsC.png")
#confirms graph was made
plt.close()
print("Scatter plot of A vs C is made")

#Same as the last graph but with A vs E data
plt.scatter(x = A_E["log2FoldChange"], y = (A_E["padj"]).apply(lambda x:-np.log10(x)),s=1,label= "No change")

plt.scatter(x = Up_A_E["log2FoldChange"], y = (Up_A_E["padj"]).apply(lambda x: -np.log10(x)), s = 1,label = "Up regulated", color = "red")
# s = size of points

plt.scatter(x = Dn_A_E["log2FoldChange"], y = (Dn_A_E["padj"]).apply(lambda x: -np.log10(x)), s = 1 ,label = "Down regulated", color = "blue")

plt.xlabel("Log2 Fold Change")
plt.ylabel("- log10 Adjusted P-value")
plt.title("A vs E differentially expressed genes")
plt.legend()
#Save png of A vs E graph
plt.savefig("Scatter_AvsE.png")
#confirms A vs E graph was made
plt.close()
print("Scatter plot of A vs E is made")

#log2 of base mean is being used as a proxy for mean expression
plt.scatter(x= np.log2(A_C["baseMean"]), y = A_C["log2FoldChange"] , s = 1)
plt.xlabel("Mean log fold change")
plt.ylabel("Log 2 fold change")
plt.title("Relationship between mean expression of log 2 fold change in A vs C")
#Save A vs C MA plot
plt.savefig('MA_AvsC.png')
plt.close()
#confirms graph was made
print("MA plot of A vs C is made")

#same as last MA plot however with A_E data
plt.scatter(x= np.log2(A_E["baseMean"]), y = A_E["log2FoldChange"] , s = 1)

plt.xlabel("Mean log fold change")

plt.ylabel("Log 2 fold change")

plt.title("Relationship between mean expression of log 2 fold change in A vs E")
#Saves A vs E MA plot
plt.savefig('MA_AvsE.png')
plt.close()
#Confirms graph was made
print("MA plot of A vs E is made")

#distribution of p-values for the A vs C dataset
#30 bins chosen since it captured the distribution of p-values the best without the data being too clumped or too spread out
plt.hist(A_C["pvalue"], bins = 30)
plt.xlabel("P values")
plt.ylabel("Frequency")
plt.title("Distribution of P-values in A vs C")
#Saves A vs C histogram distribution of p values
plt.savefig('Histogram_pvalues_AvsC.png')
#confirms graph was made
plt.close()
print("Histogram of A vs C is made")

#Same as last histogram but with A vs E dataset
plt.hist(A_E["pvalue"], bins = 30)
plt.xlabel("P values")
plt.ylabel("Frequency")
plt.title("Distribution of P-values in A vs E")
#Saves histogram of A vs E distribution of p-values
plt.savefig('Histogram_pvalues_AvsE.png')
plt.close()
print("Histogram of A vs E made")


#Creation of heatmap comparing the top differentially expressed genes for each dataset

#Keeps data for A vs E and A vs C which have corresponding P-values of < 0.001
#0.05 is the standard value however it is important only strongest p-values are kept
diff_A_C = A_C[(A_C["padj"] < 0.001)][[ "log2FoldChange","padj","gene_id"]]
diff_A_E = A_E[(A_E["padj"] < 0.001)][["log2FoldChange", "padj", "gene_id"]]


#.head(50) takes the top 10 adjusted p-values
#The default for sort_value is in accending order so lowesest p-values at the top of the dataframe
#because of this .head(10) can be used to take the top 50 lowest p-values
#repeated for both datasets
Top_A_C = diff_A_C.sort_values("padj").head(10)
Top_A_E = diff_A_E.sort_values("padj").head(10)
#A_E log2FoldChange renamed to help seperate the two datasets when the heatmap is made
Top_A_E=Top_A_E.rename(columns = {"log2FoldChange" : "log2FoldChange A vs E"})
#creates list of lists from the dataframes
diff = [Top_A_C, Top_A_E]
#Combines the two lists into single dataframe with the concat function
#keys from concat allow for deliniation of the two datasets
Fdiff = pd.concat(diff, keys = ["A vs C", "A vs E"])

#creates dataframe that only keeps the columns of interest since padj values have now been sorted
Finaldiff = Fdiff[["gene_id", "log2FoldChange","log2FoldChange A vs E"]]

#Setting gene_id to index prevents error of being unable to convert the gene_ids from str to float
Finaldiff.set_index('gene_id', inplace  = True)
#produces heat map from the combind dataframe comparing the top differentially expressed genes from each of the datasets 
sns.heatmap(Finaldiff)
plt.xlabel("Fold change of A vs C and A vs E")
plt.ylabel("Gene id")
plt.title("top differentially Expressed genes for each dataset")
#Saves heatmap of AvC and AvsE combined heat map
plt.savefig('heatmap.png')
plt.close()
print("heatmap created")

#creates dataframes of significantly differentially expressed genes with corresponding important columns
#intermediate step of creating the list of these values
inter_A_C = A_C[(A_C["padj"] < 0.05)][["gene_id","log2FoldChange","pvalue","padj"]]
inter_A_E = A_E[(A_E["padj"] < 0.05)][["gene_id","log2FoldChange","pvalue","padj"]]

#takes the created intermediate dataframe containing the important information and converts the df into lists
list_A_C = inter_A_C.values.tolist()
list_A_E = inter_A_C.values.tolist()

#Functional analysis of the top differentially expressed gene for each dataset

#similar to above steps data filtered so only significantly differentially expressed genes are present
#sort values organises data into the top differentially expressed genes with head keeping only the top 50 genes
inter2A_C = A_C[(A_C["padj"] < 0.05)][["gene_id", "padj"]]

inter2A_C = inter2A_C.sort_values("padj").head(50)

inter2A_E = A_E[(A_E["padj"] < 0.05)][["gene_id", "padj"]]

inter2A_E= inter2A_E.sort_values("padj").head(50)

#"adjusted p-values removed so only relevant information is present
#converted to list so Gprofiler can understand the input
ProfilerA_C = inter2A_C[["gene_id"]]

ProfilerA_C = ProfilerA_C.values.tolist()

ProfilerA_E = inter2A_E[["gene_id"]]

ProfilerA_E = ProfilerA_E.values.tolist()

#Gprofiler used to peform funcitonal analysis 
gp = GProfiler(return_dataframe=True)
print(gp.convert(organism='scerevisiae',            #sacchromyces used since that is the species data came from
            query= ProfilerA_C[0],            #takes the tops most differentially expressed gene from the dataset
            target_namespace='ENSG'))         # converts input from ENTREZGENE_ACC numeric input to ENSG

#above steps repeated for A vs E data set
gp = GProfiler(return_dataframe=True)
print(gp.convert(organism='scerevisiae',
            query= ProfilerA_E[0],
            target_namespace='ENSG'))
