############################# IMPORTS ####################################################################
import re
import numpy as np
import pandas as pd
import argparse
import logging
import sys
import os
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as shc
import shutil
import time
from scipy.spatial.distance import squareform
import concurrent.futures
from sklearn.cluster import AffinityPropagation
from buildcontactmaps import list_contacts_per_snapshot, list_residues, list_matrix
from scipy.spatial import distance


############################# K-MEANS IMPLEMENTATION #######################################################

# class used for clustering
class affinity_propagation:
    start = time.time()

    def __init__(self, indexes, damping, max_iter, toll, random_state, number_residues, matrices, similarity_matrix):
        self.indexes = indexes
        self.damping = damping
        self.max_iter = max_iter
        self.toll = toll
        self.random_state = random_state
        self.number_residues = number_residues
        self.matrices = matrices
        self.similarity_matrix = similarity_matrix

    # list of labels indicating which cluster the element in that position belongs to
    def getElements(self, labels, size):
        # empty dictionary
        clusters = {}
        # creating size lists to represent the clusters
        for i in range(size):
            clusters[i] = []
        # for every element add the element to the cluster it belongs to
        for index in self.indexes:
            clusters[labels[index]].append(index)
        # return cluster list
        return clusters

    def clusterize(self):
        #TODO: remove if useless
        #pref = np.min(self.similarity_matrix)

        ap = AffinityPropagation(damping=self.damping, affinity='precomputed', max_iter=self.max_iter,
                                 random_state=self.random_state, convergence_iter=self.toll)
        # running clustering algorithm and getting a list indicating which cluster the element
        # in that position belongs to
        labels = ap.fit_predict(self.similarity_matrix)
        # getting the indexes indicating the elements that are the cluster centers
        centers = ap.cluster_centers_indices_
        # number of clusters
        size = len(centers)
        # repartition elements in the cluster
        clusters = self.getElements(labels, size)

        return clusters, centers

############################# MAIN #########################################################################


def main():
    ############################# HELPER FUNCTIONS  ###########################################################

    # function for parsing the values in the config file, path is the path to the config file
    def parseConfigFile(path):
        List = []
        # open file
        with open(path) as f:
            # for every line
            for line in f:
                # ignore comments or empty lines
                if not line.startswith(('#', '\n')):
                    # eliminate spaces, new lines and split values from names
                    x, y = line.replace(' ', '').replace('\n', '').split('=')
                    List.append((x, y))
        # return all values
        return np.array(List)[:, 1]

    # Given an array of points or a distance matrix and the labels associated to the elements,
    # this function creates a dendrogram.
    # We can also specify:
    # - the type of link
    # - the orientation for the dendrogram
    # - the truncation and truncation mode
    # - color threshold
    # - labels font size
    # - optimal sorting for leaves
    # - which child to plot first in the dendrogram
    # - if the dendrogram is truncated we can see how many elements are not shown
    def make_dendrogram(array, labels, output_path, method, optimal_ordering, orientation, distance_sort,
                        show_leaf_counts, truncation, truncate_mode, color_threshold, protein_name, font_size):
        plt.figure(figsize=(20, 10))
        shc.dendrogram(shc.linkage(array, method=method, optimal_ordering=optimal_ordering),
                       orientation=orientation,
                       labels=labels,
                       distance_sort=distance_sort,
                       show_leaf_counts=show_leaf_counts,
                       p=truncation,
                       truncate_mode=truncate_mode,
                       color_threshold=color_threshold,
                       leaf_font_size=font_size
                       )
        plt.savefig(os.path.join(output_path, protein_name + '_dendrogram.png'))

    # given two indexes we get the corresponding csv file, we compare them and save in a list where
    # they differ
    def get_differences(indx_1, indx_2, list_paths, number_residues):
        snap_1 = pd.read_csv(os.path.join(list_paths[indx_1]))
        snap_1 = np.array(snap_1)
        snap_2 = pd.read_csv(os.path.join(list_paths[indx_2]))
        snap_2 = np.array(snap_2)
        residues = snap_1[:, [0]]
        snap_1 = np.delete(snap_1, 0, 1)
        snap_2 = np.delete(snap_2, 0, 1)
        list = []
        for i in range(0, number_residues-1):
            for j in range(i, number_residues):
                if snap_1[i, j] != snap_2[i, j]:
                    list.append((residues[i][0], residues[j][0]))
        return list

    # given the indexes of the clusters centers, we confront them to see in which residues
    # they differ
    def confront_clusters(centers, list_paths, number_residues):
        elements = len(centers)
        differences = {}
        for el_1 in range(0, elements - 1):
            for el_2 in range(1, elements):
                if el_1 != el_2:
                    differences[(el_1, el_2)] = get_differences(el_1, el_2, list_paths, number_residues)
        return differences

    ############################# ARGUMENT PARSER ############################################################

    inizio = time.time()
    parser = argparse.ArgumentParser(description='Arguments for clustering contact maps.')
    parser.add_argument('input_path', metavar='contact_maps_file', type=str,
                        help='the directory which contains the paths to RING contact map files (edge files), one '
                             'file per line')
    parser.add_argument('output_path', metavar='output_path', type=str,
                        help='the path of the output directory in which to save the results')
    parser.add_argument('config_file', metavar='config_file', help='configuration file for clustering algorithm',
                        type=str)
    parser.add_argument('temporary_path', metavar='temporary_path', help='the path of the temporary directory in '
                                                                         'which to save the contact matrixes',
                        type=str)
    args = parser.parse_args()

    # show message to user
    print('parsing program variables...')

    # if output directory doesn't already exist, we create it
    if not os.path.exists(args.output_path):
        os.mkdir(args.output_path)
    if not os.path.isdir(args.output_path):
        logging.error('output_path not a directory path, please check path', 'output_path', args.output_path)
        sys.exit(1)

    # if temporary directory doesn't already exist, we create it
    if not os.path.exists(args.temporary_path):
        os.mkdir(args.temporary_path)
    if not os.path.isdir(args.temporary_path):
        logging.error('temporary_path not a directory path, please check path', 'temporary_path', args.temporary_path)
        sys.exit(1)

    # TODO: substitute path with file containing residues paths
    # if temporary directory doesn't already exist we show an error
    if not os.path.exists(args.input_path):
        logging.error('input_path directory does not exists, please check path', 'input_path', args.input_path)
        sys.exit(1)

    # if config file isn't found or does not exist we show an error
    if not os.path.isfile(args.config_file):
        logging.error('configuration file not found, please check path', 'config_file', args.config_file)
        sys.exit(1)

    ############################# CONFIG FILE PARSING #########################################################

    # saving protein name (name on folder before edges folder)
    protein_name = args.input_path.split('/')[-2]

    # reading all the variables
    try:
        contact_treshold, energy_treshold, method, optimal_ordering, truncation, truncate_mode, \
        color_threshold, orientation, distance_sort, show_leaf_counts, \
        font_size, damping, random_state, toll, max_iter = parseConfigFile(args.config_file)
    except ValueError:
        logging.error("Couldn't parse configuration file, please check",
                      ValueError.args, 'config_file', args.config_file)
        sys.exit(1)

    # checking if values are correct, otherwise show user an error
    try:
        contact_treshold = float(contact_treshold)
    except ValueError:
        logging.error("In config file contact_treshold value is not convertible to float, please check",
                      ValueError.args, 'contact_treshold', contact_threshold)
        sys.exit(1)

    try:
        energy_treshold = float(energy_treshold)
    except ValueError:
        logging.error("In config file energy_treshold value is not convertible to float, please check",
                      ValueError.args, 'energy_treshold', energy_treshold)
        sys.exit(1)

    if method not in ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']:
        logging.error("In config file method value is not correct, please check", 'method', method)
        sys.exit(1)

    try:
        optimal_ordering = bool(optimal_ordering)
    except ValueError:
        logging.error("In config file optimal_ordering value is not convertible to bool, please check",
                      ValueError.args, 'optimal_ordering', optimal_ordering)
        sys.exit(1)

    try:
        truncation = int(truncation)
    except ValueError:
        logging.error("In config file truncation value is not convertible to int, please check",
                      ValueError.args, 'truncation', truncation)
        sys.exit(1)

    if truncate_mode not in ['none', 'lastp', 'level']:
        logging.error("In config file truncate_mode value is not correct, please check", 'truncate_mode', truncate_mode)
        sys.exit(1)

    try:
        color_threshold = float(color_threshold)
    except ValueError:
        logging.error("In config file color_threshold value is not convertible to float, please check",
                      ValueError.args, 'color_threshold', color_threshold)
        sys.exit(1)

    if orientation not in ['top', 'bottom', 'left', 'right']:
        logging.error("In config file orientation value is not correct, please check", 'orientation', orientation)
        sys.exit(1)

    if orientation not in ['ascending', 'descending']:
        try:
            distance_sort = bool(distance_sort)
        except ValueError:
            logging.error(
                "In config file distance_sort value is not convertible to bool or is not correct, please check",
                ValueError.args, 'distance_sort', distance_sort)
            sys.exit(1)

    try:
        show_leaf_counts = bool(show_leaf_counts)
    except ValueError:
        logging.error("In config file show_leaf_counts value is not convertible to bool, please check",
                      ValueError.args, 'show-leaf-counts', show_leaf_counts)
        sys.exit(1)

    try:
        font_size = int(font_size)
    except ValueError:
        logging.error("In config file font_size value is not convertible to int, please check",
                      ValueError.args, 'font_size', font_size)
        sys.exit(1)

    try:
        damping = float(damping)
    except ValueError:
        logging.error("In config file damping value is not convertible to float, please check",
                      ValueError.args, 'damping', damping)
        sys.exit(1)

    try:
        random_state = int(random_state)
    except ValueError:
        logging.error("In config file random_state value is not convertible to int, please check",
                      ValueError.args, 'random_state', random_state)
        sys.exit(1)

    try:
        toll = int(toll)
    except ValueError:
        logging.error("In config file toll value is not convertible to int, please check",
                      ValueError.args, 'toll', toll)
        sys.exit(1)

    try:
        max_iter = int(max_iter)
    except ValueError:
        logging.error("In config file max_iter value is not convertible to int, please check",
                      ValueError.args, 'max_iter', max_iter)
        sys.exit(1)

    ############################# READ INPUT ###############################################################

    print('reading contact maps...')
    # all files names contained in input_path are added in a list
    # not ordered
    list_names = [f for f in listdir(args.input_path) if isfile(join(args.input_path, f))]

    # the total number of edges file
    number_files = len(list_names)

    # get a list of dictionaries that contain the contats for all files
    # each cell contains a dictionary that contains all residues contacts in a file
    list_contacts = list_contacts_per_snapshot(contact_treshold, energy_treshold, list_names, args.input_path)

    # get a list that contains all residues contained in a file
    # the list is ordered
    listRes = list_residues(number_files, list_contacts)
    number_residues = len(listRes)

    # matrices rappresent a contact map
    listMatrix = list_matrix(number_files, list_contacts, listRes)

    # construisco i dataframe e li salvo in una cartella
    for i in range(0, number_files):
        ((pd.DataFrame(listMatrix[i], columns=listRes, index=listRes)).fillna(0)) \
            .to_csv((os.path.join(args.temporary_path, list_names[i] + '.csv')), index=listRes, header=True)

    ############################# DISTANCE MATRIX ################################################

    print('making distance matrix...')
    # ordered crescently by snapshot number
    list_names.sort(key=lambda f: int(re.sub('\D', '', f)))
    listPaths = [os.path.join(args.temporary_path, list_names[i] + '.csv') for i in range(number_files)]

    # matrices are also ordered by snapshot
    matrices = []
    for path in listPaths:
        m = pd.read_csv(path, index_col=0)
        matrices.append(m.values)

    # NXN matrix with all zeros
    distance_matrix = np.zeros((number_files, number_files), dtype=float)

    # # insert Jaccard distance into each element of the matrix dist, we need a symmetric matrix
    # for i in range(1, len(listPaths)):
    #     f_i = matrices[i].flatten()
    #     for j in range(0, i):
    #         f_j = matrices[j].flatten()
    #         distance_matrix[i, j] = distance.jaccard(f_i, f_j)
    #         distance_matrix[j, i] = distance_matrix[i, j]
    #
    #
    # # create DataFrame for distance_matrix
    # distance_df = pd.DataFrame(distance_matrix, columns=list_names, index=list_names)
    # # export DataFrame into file csv
    # distance_df.to_csv(os.path.join(args.output_path, protein_name + '_distance_matrix.csv'))

    # # TODO:remove
    # READ MATRIX FROM FILE
    distance_matrix = pd.read_csv(os.path.join(args.output_path, protein_name + '_distance_matrix.csv'), index_col=0)
    distance_matrix = np.array(distance_matrix)

    # making similarity matrix from distance matrix. The similarity matrix is also called affinity matrix
    affinity_matrix = 1 - distance_matrix

    # indexes of snapshots (0 -> snapshot 1, 1-> snapshot 2, ...)
    indexes = [i for i in range(number_files)]
    # labels of snapshots (1 -> snapshot 1, 2-> snapshot 2, ...) for dendrogram
    labels = [i + 1 for i in range(number_files)]

    ############################# DENDROGRAM ##################################################################

    # show message to user
    print('making dendrogram...')

    # call function to create dendrogram
    make_dendrogram(squareform(distance_matrix), labels, args.output_path, method, optimal_ordering, orientation,
                    distance_sort,
                    show_leaf_counts, truncation, truncate_mode, color_threshold, protein_name, font_size)

    ############################# CLUSTERING ############################################################

    # message to user
    print('clustering...')

    # creating instance of clustering class
    clustering = affinity_propagation(indexes, damping, max_iter, toll, random_state, number_residues,
                                      matrices, affinity_matrix)

    # running clustering algorithm and getting results
    final_clusters, final_centers = clustering.clusterize()

    ############################# PROCESSING RESULTS ############################################################

    print('working on outputs...')

    # # extract samples from the clusters
    # for index in final_clusters:
    #     sample = final_centers[index]
    #     shutil.copy(listPaths[sample], os.path.join(args.output_path, protein_name + '_Cluster_{}_'
    #                                                 .format(index) + 'extract.csv'))

    # output file with clustering results
    file = open(os.path.join(args.output_path, protein_name + '_cluster_results.txt'), 'w+')
    file.write('CLUSTERING RESULTS \n')
    # dictionary indices
    for index in final_clusters:
        file.write('CLUSTER {}:\n'.format(index))
        # elements in dictionary
        for cluster in final_clusters[index]:
            mat = matrices[cluster].flatten()
            cent = matrices[final_centers[index]].flatten()
            file.write('Element: ' + str(list_names[cluster]) + ', distance from center: ' +
                       str(distance.jaccard(mat, cent)) + '\n')
    file.close()

    # output file with clusters differences
    file = open(os.path.join(args.output_path, protein_name + '_cluster_differences.txt'), 'w+')
    file.write('CLUSTERS DIFFERENCES \n')
    differences = confront_clusters(final_centers, listPaths, number_residues)
    # dictionary indices
    for index in differences:
        file.write('CLUSTER {} - CLUSTER {}:\n'.format(index[0], index[1]))
        # elements in dictionary
        for elements in differences[index]:
            file.write(str(elements) + '\n')
    file.close()



    # TODO: remove tempo
    print('tempo totale = ' + str(time.time() - inizio))
    sys.exit(0)


if __name__ == "__main__":
    main()
