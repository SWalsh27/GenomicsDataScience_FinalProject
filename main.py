import numpy as np
from scipy.spatial.distance import jaccard
from sklearn.metrics import jaccard_score
import warnings
from rdkit import Chem, DataStructs
import networkx as nx
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MultiLabelBinarizer
import pandas as pd
from networkx.algorithms.bipartite.matrix import biadjacency_matrix
import matplotlib as mt
warnings.filterwarnings('ignore')

intersect_data = pd.read_csv("Datasets/intersection_Fingerprint_mat.tsv", sep='\t')
side_effect_data = pd.read_csv("Datasets/side-effect-and-drug_name_upper.tsv",sep='\t')

df = pd.read_csv("Datasets/side-effect-and-drug_name_upper.tsv",sep='\t')
drug_id = df["drugbank_id"]
drug_name = df["drugbank_name"]
side_effect = df["side_effect_name"]

edgelist1 = zip(side_effect, drug_name)

B = nx.DiGraph()
B.add_nodes_from(side_effect,bipartite=0)
B.add_nodes_from(drug_name,bipartite=1)
B.add_edges_from(edgelist1)

drug_nodes = {n for n, d in B.nodes(data=True) if d['bipartite'] == 1}
side_effect_nodes = {n for n, d in B.nodes(data=True) if d['bipartite'] == 0}
drug_nodes = list(drug_nodes)
drug_nodes.sort()
side_effect_nodes = list(side_effect_nodes)
side_effect_nodes.sort()

matrix_all = biadjacency_matrix(B, row_order= side_effect_nodes,column_order = drug_nodes)
matrix_A = matrix_all.A
m_all,n_all = matrix_all.shape

print(matrix_A)
print(matrix_all.shape)
def Logistic_Regression_Model(dataset):
    data = dataset.drop(dataset.columns[[0,1,2]],axis=1)
    MultiLabelBinarizer().fit(data)
    side_effect = MultiLabelBinarizer().fit_transform(data.transpose().values)

    print(side_effect)
    print(side_effect.shape)
    #data.to_csv('Datasets/new_data.csv',index=False)

#Logistic_Regression_Model(side_effect_data)

matrix_a = {1,0,0,1,1,1}
matrix_b = {0,0,1,1,1,0}
matrix_binaryA = np.array([1,0,0,1,1,1])
matrix_binaryB = np.array([0,0,1,1,1,0])

# def jaccard_similarity(A,B):
#
#     # Find intersection of two sets
#     nominator = A.intersection(B)
#
#     # Find union of two sets
#     denominator = A.union(B)
#
#     # Take the ratio of sizes
#     similarity = len(nominator) / len(denominator)
#
#     return similarity

# Jaccard similarity function
def jaccard_similarity(A,B):
    similarity_coefficient = jaccard_score(matrix_A, matrix_A, average='micro')
    return similarity_coefficient



def jaccard_distance(A, B):
    # Find symmetric difference of two sets
    nominator = A.symmetric_difference(B)

    # Find union of two sets
    denominator = A.union(B)

    # Take the ratio of sizes
    distance = len(nominator) / len(denominator)

    return distance

Jaccard_similarity = jaccard_similarity(matrix_A,matrix_A)
distance_metric = jaccard_distance(matrix_a,matrix_b)
binary_distance = jaccard(matrix_binaryA,matrix_binaryB)
similarity = jaccard_similarity(matrix_a,matrix_b)

print("real number distance ", distance_metric)
print("Jaccard similarity real number R", similarity)
print("Jaccard similarity is", Jaccard_similarity)
print("Binary Distance for jaccard_distance is", binary_distance)





#Tanimoto coefficient



#Pearson's Correlation











