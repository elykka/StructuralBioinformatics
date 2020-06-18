from sklearn.cluster import AffinityPropagation
import pandas as pd
import os
import numpy as np
import tmscoring
import re
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as shc


def make_dendrogram(array, labels, output_path, method, optimal_ordering, orientation, distance_sort,
                    show_leaf_counts, truncation, truncate_mode, color_threshold, protein_name):
    plt.figure(figsize=(20, 10))
    shc.dendrogram(shc.linkage(array, method=method, optimal_ordering=optimal_ordering),
                   orientation=orientation,
                   labels=labels,
                   distance_sort=distance_sort,
                   show_leaf_counts=show_leaf_counts,
                   p=truncation,
                   truncate_mode=truncate_mode,
                   color_threshold=color_threshold,
                   leaf_font_size=8
                   )
    plt.savefig(os.path.join(output_path, protein_name + '_RMSD_dendrogram.png'))

def getElements(labels, size, indexes):
    elements = {}
    for i in range(size):
        elements[i] = []
    for index in indexes:
        elements[labels[index]].append(index)
    return elements

output_path = '/home/elykka/Documents/CMC/output'
input_path = '/home/elykka/Documents/CMC/data/antibody/pdb'

list_names = [os.path.join(input_path, f) for f in listdir(input_path) if isfile(join(input_path, f))]
list_names.sort(key=lambda f: int(re.sub('\D', '', f)))
print(list_names)
 # the total number of edges file
number_files = len(list_names)
protein_name = input_path.split('/')[-2]

# NXN matrix with all zeros
distance_matrix = np.zeros((number_files, number_files), dtype=float)
#
# insert Jaccard distance into each element of the matrix dist
for i in range(1, number_files):
    for j in range(0, i):
        alignment = tmscoring.TMscoring(list_names[i], list_names[j])
        alignment.optimise()
        distance_matrix[i, j] = alignment.rmsd(**alignment.get_current_values())
        distance_matrix[j, i] = distance_matrix[i, j]
# create DataFrame
distance_df = pd.DataFrame(distance_matrix, columns=list_names, index=list_names)
# export DataFrame into file csv
distance_df.to_csv(os.path.join(output_path, protein_name + '_RMSD_distance_matrix.csv'))

distance_matrix = pd.read_csv(os.path.join(output_path, protein_name + '_RMSD_distance_matrix.csv'), index_col=0)
distance_matrix = np.array(distance_matrix)

affinity_matrix = 1 - distance_matrix
indexes = [i for i in range(number_files)]
labels = [i+1 for i in range(number_files)]

make_dendrogram(distance_matrix, labels, output_path, 'single', True, 'top',
                'descending', True, 50, 'none', 3.0, protein_name)


ap = AffinityPropagation(damping=0.5, affinity='precomputed', random_state=1)
labels = ap.fit_predict(affinity_matrix)
centers = ap.cluster_centers_indices_
size = len(centers)
clusters = getElements(labels, size, indexes)

# output file with clustering results
file = open(os.path.join(output_path, protein_name + '_RMSD_cluster_results.txt'), 'w+')
file.write('CLUSTERING RESULTS \n')
for index in clusters:
    file.write('CLUSTER {}:\n'.format(index))
    for cluster in clusters[index]:
        file.write('Element: ' + str(list_names[cluster] + '\n'))
file.close()