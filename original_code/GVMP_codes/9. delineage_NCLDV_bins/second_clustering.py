#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy import genfromtxt, matrix
from typing import List
from sklearn import preprocessing
from sklearn.metrics.pairwise import euclidean_distances as d
from sklearn import manifold
from scipy.cluster.hierarchy import linkage, to_tree, ClusterNode
import sys
import os
import glob

def dendrogram2newick(
    node: ClusterNode, parent_dist: float, leaf_names: List[str], newick: str = ""
) -> str:
    """
    Convert scipy dendrogram tree to newick format tree

    Args:
        node (ClusterNode): Tree node
        parent_dist (float): Parent distance
        leaf_names (List[str]): Leaf names
        newick (str): newick format string (Used in recursion)

    Returns:
        str: Newick format tree
    """
    if node.is_leaf():
        return f"{leaf_names[node.id]}:{(parent_dist - node.dist):.2f}{newick}"
    else:
        if len(newick) > 0:
            newick = f"):{(parent_dist - node.dist):.2f}{newick}"
        else:
            newick = ");"
        newick = dendrogram2newick(node.left, node.dist, leaf_names, newick)
        newick = dendrogram2newick(node.right, node.dist, leaf_names, f",{newick}")
        newick = f"({newick}"
        return newick

binid = str(sys.argv[1])

# Input coverage matrix files
labels = genfromtxt(f'/path/to/coverage_matrix/{binid}.depth_modified.txt', usecols=0, dtype=str)
vectors = genfromtxt(f'/path/to/coverage_matrix/{binid}.depth_modified.txt', dtype=np.float64)[:,1:]
#make contig id be index for each row
data = {label: row for label, row in zip(labels, vectors)}

#modify zero values
zero_filter = vectors == 0
vectors[zero_filter] = 0.0000000000000001

normalizer = preprocessing.Normalizer(norm='l1')
normalized_vectors = normalizer.fit_transform(vectors)

similarities = d(normalized_vectors)

mds = manifold.MDS(max_iter=300, eps=1e-10, dissimilarity="precomputed", n_jobs=1)
coverage_scaled_vectors = mds.fit(similarities).embedding_
# Input tetra matrix files
labels = genfromtxt(f'/path/to/tetra_matrix/{binid}_tetra_count_normalized.csv', delimiter=',', usecols=0, dtype=str)
vectors = genfromtxt(f'/path/to/tetra_matrix/{binid}_tetra_count_normalized.csv', delimiter=',', dtype=np.float64)[:,1:]
#make contig id be index for each row
data = {label: row for label, row in zip(labels, vectors)}

#modify zero values
zero_filter = vectors == 0
vectors[zero_filter] = 0.0000000000000001

normalizer = preprocessing.Normalizer(norm='l1')
normalized_vectors = normalizer.fit_transform(vectors)

similarities = d(normalized_vectors)

mds = manifold.MDS(max_iter=300, eps=1e-10, dissimilarity="precomputed", n_jobs=1)
tetra_scaled_vectors = mds.fit(similarities).embedding_
# Combine coverage and tetranucleotide matrix for clustering
combined_matrix = np.hstack([coverage_scaled_vectors, tetra_scaled_vectors])

dict = {label: row for label, row in zip(labels, combined_matrix)}

name_list = [k for k in dict.keys()]
linkage = linkage(combined_matrix, metric='euclidean', method='ward')
tree = to_tree(linkage, rd=False)
# dendrogram transform to newick
newick = dendrogram2newick(tree, tree.dist, name_list)

with open(f'/path/to/newick/{binid}.newick', 'w') as F:
    F.write(newick)