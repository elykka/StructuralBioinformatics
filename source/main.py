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
import math
import random
import time
from buildcontactmaps import list_contacts_per_snapshot, list_residues, list_matrix
from builddistancematrix import distanceJaccard


############################# K-MEANS IMPLEMENTATION #######################################################


class k_means:
    def __init__(self, metric, average_func, toll, max_iter):
        self.toll = toll
        self.max_iter = max_iter
        self.metric = metric
        self.average_func = average_func

    def clusterize(self, data, k, distance_matrix, array_elements):

        self.k = k

        self.centers = {}
        for i in range(self.k):
            self.centers[i] = array_elements[data[i]]

        for i in range(self.max_iter):
            self.clusters = {}

            for i in range(self.k):
                self.clusters[i] = []

            # Find the distance of each point in the data set with the k centers
            for element in data:
                distances = []
                for center in self.centers:
                    # Compute the distances from the element to the center
                    distances.append(self.metric(element, self.centers.get(center), distance_matrix, array_elements))
                    # Assign the element to the closest center, index returns the center
                cluster_center = distances.index(min(distances))
                self.clusters[cluster_center].append(element)

            all_centers = dict(self.centers)

            # Find the new centers by taking the average of the points in each cluster group.
            for cluster in self.clusters:
                self.centers[cluster] = self.average_func(self.clusters[cluster], array_elements)

            stop = True

            for center in self.centers:
                original_center = all_centers[center]
                current_center = self.centers[center]
                # if the centers difference is higher than the tollerance treshold we continue the execution
                if np.sum(current_center - original_center) > self.toll:
                    stop = False

            if stop:
                break

    def optimal(self, k_range, user_function, average, indexes, distance_matrix, data, toll, max_iter):
        elbow_values = []
        optimal_clusters = []
        optimal_centers = []
        KMeans = k_means(user_function, average, toll, max_iter)
        for k in k_range:
            KMeans.clusterize(indexes, k, distance_matrix, data)
            clust = KMeans.clusters
            sum = 0.0
            for index in clust:
                value = 0.0
                for i in range(len(clust[index]) - 1):
                    for j in range(i + 1, len(clust[index])):
                        value += user_function(clust[index][i], clust[index][j], distance_matrix, data)
                sum += value / len(clust[index])
            elbow_values.append(sum)
            optimal_clusters.append(KMeans.clusters)
            optimal_centers.append(KMeans.centers)
        return elbow_values, optimal_clusters, optimal_centers


############################# MAIN #########################################################################


def main():


############################# HELPER FUNCTIONS  ###########################################################

    def parseConfigFile(path):
        List = []
        with open(path) as f:
            for line in f:
                if not line.startswith(('#', '\n')):
                    x, y = line.replace(' ', '').replace('\n', '').split('=')
                    List.append((x, y))
        return np.array(List)[:, 1]


    # given an array of points or a distance matrix and the labels associated to the elements,
    # this funtion creates a dendrogram.
    # we can also specify:
    # the type of link
    def make_dendrogram(array, labels, output_path, method, optimal_ordering, orientation, distance_sort,
                        show_leaf_counts, truncation, truncate_mode, color_threshold):
        plt.figure(figsize=(50, 30))
        shc.dendrogram(shc.linkage(array, method=method, optimal_ordering=optimal_ordering),
                       orientation=orientation,
                       labels=labels,
                       distance_sort=distance_sort,
                       show_leaf_counts=show_leaf_counts,
                       p=truncation,
                       truncate_mode=truncate_mode,
                       color_threshold=color_threshold
                       )
        plt.savefig(os.path.join(output_path, 'dendrogram.png'))

    def distancePointLine(x1, y1, a, b, c):
        return abs((a * x1 + b * y1 + c)) / math.sqrt(a * a + b * b)

    def findOptimalValue(elbow_values, k_values):
        lst = len(k_values) - 1
        a = elbow_values[0] - elbow_values[lst]
        b = k_values[lst] - k_values[0]
        c1 = k_values[0] * elbow_values[lst]
        c2 = k_values[lst] * elbow_values[0]
        c = c1 - c2
        distances = []
        for k in range(len(k_values)):
            distances.append(distancePointLine(k_values[k], elbow_values[k], a, b, c))
        return distances.index(max(distances))


############################# ARGUMENT PARSER ############################################################


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

    # if temporary directory doesn't already exist we show an error
    if not os.path.exists(args.input_path):
        logging.error('input_path directory does not exists, please check path', 'input_path', args.input_path)
        sys.exit(1)

    # if config file isn't found or does not exist we show an error
    if not os.path.isfile(args.config_file):
        logging.error('configuration file not found, please check path', 'config_file', args.config_file)
        sys.exit(1)


############################# CONFIG FILE PARSING #########################################################


    # reading all the variables
    contact_treshold, energy_treshold, method, optimal_ordering, truncation, truncate_mode, \
    color_threshold, orientation, distance_sort, show_leaf_counts, \
    rangeK, toll, max_iter = parseConfigFile(args.config_file)

    # checking if they are correct


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

    if 'to' in rangeK:
        x, y = rangeK.split('to')
        try:
            x = int(x)
            y = int(y)
            if x != 0 and x < y:
                rangeK = range(x, y)
            else:
                logging.error(
                    "In config file rangeK value is not correct, x cannot be 0 and x must be inferior to y value, please check",
                    '(x,y)', (x, y))
                sys.exit(1)

        except ValueError:
            logging.error("In config file rangeK value is not convertible to a range please check",
                          ValueError.args, '(x,y)', (x, y))
            sys.exit(1)
    else:
        try:
            rangeK = int(rangeK)
            if rangeK == 0:
                logging.error("In config file rangeK value is not correct, x cannot be 0, please check",
                              'rangeK', rangeK)
                sys.exit(1)
        except ValueError:
            logging.error("In config file rangeK value is not convertible to int please check",
                          ValueError.args, 'rangeK', rangeK)
            sys.exit(1)

    try:
        toll = float(toll)
    except ValueError:
        logging.error("In config file toll value is not convertible to float, please check",
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
        ((pd.DataFrame(listMatrix[i], columns=listRes, index=listRes)).fillna(0))\
            .to_csv((os.path.join(args.temporary_path, list_names[i]+'.csv')), index=listRes, header=True)




############################# DISTANCE MATRIX ################################################

    print('making distance matrix...')
    # ordered crescently by snapshot
    list_names.sort(key=lambda f: int(re.sub('\D', '', f)))
    listPaths = [os.path.join(args.temporary_path, list_names[i]+'.csv') for i in range(number_files)]

    # matrices are also ordered by snapshot
    matrices = []
    for path in listPaths:
        m = pd.read_csv(path, index_col=0)
        matrices.append(m.values)

    # # NXN matrix with all zeros
    # distance_matrix = np.zeros((number_files, number_files), dtype=float)
    #
    # # insert Jaccard distance into each element of the matrix dist
    # for i in range(1, number_files):
    #     for j in range(0, i):
    #         distance_matrix[i, j] = distanceJaccard(matrices[i], matrices[j])
    #
    # # create DataFrame
    # distance_df = pd.DataFrame(distance_matrix, columns=list_names, index=list_names)
    # # export DataFrame into file csv
    # distance_df.to_csv(os.path.join(args.output_path, 'distance_matrix.csv'))

    distance_matrix = pd.read_csv(os.path.join(args.output_path, 'distance_matrix.csv'), index_col=0)
    distance_matrix = np.array(distance_matrix)
    distance_matrix = np.tril(distance_matrix)



############################# DENDROGRAM ##################################################################

    # Nota: quando viene passata la matrice e non un vettore, la metrica viene ignorata, non serve farne una custom
    # show message to user
    print('making dendrogram...')
    make_dendrogram(distance_matrix, list_names, args.output_path, method, optimal_ordering, orientation,
                    distance_sort,
                    show_leaf_counts, truncation, truncate_mode, color_threshold)


############################# CLUSTERING ############################################################

    start = time.time()
    # message to user
    print('clustering...')

    # custom metric for calculations between indexes and centers
    # we chose to use indexes and not the matrices directly to sped up the calculations during the
    # elbow method by accessing the distance matrix
    def user_function(point1, point2, distance_matrix, array_elements):

        if type(point1) is int and type(point2) is int:
            return distance_matrix[point1][point2]
        if type(point1) is int:
            return distanceJaccard(array_elements[point1], point2)
        if type(point2) is int:
            return distanceJaccard(point1, array_elements[point2])
        return distanceJaccard(point1, point2)

    # list of indexes, have to return an contact map
    def average(cluster, array_elements):
        map = np.zeros((number_residues, number_residues))
        for el in cluster:
            map = np.add(map, array_elements[el])
        return np.rint(map / len(cluster))

    clustering = k_means(user_function, average, toll, max_iter)
    indexes = [i for i in range(number_files)]

    opt_k = None
    if type(rangeK) is range:
        elbow_values, opt_clusters, opt_centers = clustering.optimal(rangeK, user_function,
                                                                     average, indexes, distance_matrix,
                                                                     matrices, toll, max_iter)
        opt_k = findOptimalValue(elbow_values, rangeK)
        final_clusters = opt_clusters[opt_k]
        final_centers = opt_centers[opt_k]
    else:
        clustering.clusterize(indexes, rangeK, distance_matrix, matrices)
        final_clusters = clustering.clusters
        final_centers = clustering.centers
    print(time.time() - start)

############################# PROCESSING RESULTS ############################################################

    print('working on outputs...')

    # extract samples from the clusters
    for index in final_clusters:
        sample = random.choice(final_clusters[index])
        shutil.copy(listPaths[sample], os.path.join(args.output_path, 'Cluster_{}_'.format(index) +
                                                    str(list_names[sample]) + '.csv'))
        # file = open(os.path.join(args.output_path, 'Cluster_{}_'.format(index) + str(list_names[sample])), 'w+')
        # file.write(str(matrices[sample]))
        # file.close()

    # output file with clustering results
    file = open(os.path.join(args.output_path, 'cluster_results.txt'), 'w+')
    file.write('CLUSTERING RESULTS ')
    file.write('for K {}\n'.format(rangeK))
    if type(rangeK) is range:
        file.write('The optimal number of clusters found was {}:\n'.format(opt_k + rangeK[0]))
    for index in final_clusters:
        file.write('CLUSTER {}:\n'.format(index))
        for cluster in final_clusters[index]:
            file.write('Element: ' + str(list_names[cluster]) + ', distance from center: ' +
                       str(user_function(cluster, final_centers[index], distance_matrix, matrices)) + '\n')

    file.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
