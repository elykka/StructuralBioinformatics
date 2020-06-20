############################# IMPORTS ####################################################################
import re
import numpy as np
import pandas as pd
import argparse
import logging
import sys
import os
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as shc
import shutil
import configparser
from distutils.util import strtobool
from scipy.spatial.distance import squareform
from sklearn.cluster import AffinityPropagation
from matrices_snapshots import list_contacts_per_snapshot, list_residues, list_matrix
from scipy.spatial.distance import cdist


############################# K-MEANS IMPLEMENTATION #######################################################

# Class used for clustering
class affinity_propagation:

    def __init__(self, indexes, damping, max_iter, toll, random_state, number_residues, matrices, similarity_matrix):
        self.indexes = indexes
        self.damping = damping
        self.max_iter = max_iter
        self.toll = toll
        self.random_state = random_state
        self.number_residues = number_residues
        self.matrices = matrices
        self.similarity_matrix = similarity_matrix

    # Given a list of labels indicating which cluster the element in that position belongs to and the number of
    # clusters, this function returns a dictionary that contains lists, where each list contains all the indexes
    # of elements that belong to a certain cluster. The labelling of clusters starts from 0.
    def getElements(self, labels, n_clusters):
        # empty dictionary
        clusters = {}
        # creating size lists to represent the clusters
        for i in range(n_clusters):
            clusters[i] = []
        # for every element add the element to the cluster it belongs to
        for index in self.indexes:
            clusters[labels[index]].append(index)
        # return cluster list
        return clusters

    # This function performs the clustering and returns a dictionary containing the clusters as lists and
    # and a list containing the clusters centers.
    # Cluster and centers do not contain the elements, but an index associated to them.
    def clusterize(self):
        ap = AffinityPropagation(damping=self.damping, affinity='precomputed', max_iter=self.max_iter,
                                 random_state=self.random_state, convergence_iter=self.toll)
        # running clustering algorithm and getting a list indicating which cluster the element
        # in that position belongs to
        labels = ap.fit_predict(self.similarity_matrix)
        # getting the indexes indicating the elements that are the cluster centers
        centers = ap.cluster_centers_indices_
        # number of clusters
        size = len(centers)
        # repartition elements in the clusters
        clusters = self.getElements(labels, size)

        return clusters, centers

############################# MAIN #########################################################################


def main():
    ############################# HELPER FUNCTIONS  ###########################################################

    # Given the energy threshold, the contact threshold and the list of all file edges, this function
    # returns the list of contacts that are to be considered for all snapshots.
    # The list returned is actually a list of dictionaries, where the key represents the residue, node_1, we found
    # involved in a contact and the value represents all the residues node_2 involved in a contact with node_1.
    def list_contacts_per_snapshot(contact_threshold, energy_threshold, listPaths):

        listCon = []
        # total number of file edges
        numFiles = len(listPaths)
        # for every edges file
        for i in range(0, numFiles):
            # dictionary containing all contacts that we will consider
            contacts = {}
            # open files for parsing
            with open(listPaths[i]) as f:
                next(f)
                # parse every line and get the residues distance and energy
                # every line is a contact
                for line in f:
                    node_1, edge, node_2, distance, _, energy, atom_1, atom_2 = line.split()[0:8]
                    distance = float(distance)
                    energy = float(energy)
                    # check if contact energy and distance values are inside the thresholds
                    if (distance <= contact_threshold) & (energy > energy_threshold):
                        edge_type, edge_loc = edge.split(":")
                        contacts.setdefault((node_1, edge_type), [])
                        contacts.setdefault((node_2, edge_type), [])

                        # check if the contact has already been added to the dictionary
                        if (node_2, edge_type) not in contacts[(node_1, edge_type)]:
                            contacts[(node_1, edge_type)].append((node_2, edge_type))

                        if (node_1, edge_type) not in contacts[(node_2, edge_type)]:
                            contacts[(node_2, edge_type)].append((node_1, edge_type))
                # add dictionary to the list of contacts
                listCon.append(contacts)
        return listCon

    # given the number of files and the list of contacts, this function returns a list with all the different
    # residues contained in the files.
    def list_residues(numFiles, listDizCont):
        res = []
        # for all files
        for i in range(0, numFiles):
            # for every residue in the first position of the contact
            for node in listDizCont[i].keys():
                numNode2 = len(listDizCont[i][node])
                # for every residue in the second position of the contact that we found
                for j in range(0, numNode2):
                    # if the contact is not already present in the list, add it
                    if listDizCont[i][node][j] not in res:
                        res.append(listDizCont[i][node][j])
        # removing duplicates
        res = list(dict.fromkeys(res))
        # numerical order
        res.sort(key=lambda x: int(x[0].split(':')[1]))
        # alphabetic order
        res.sort(key=lambda x: x[0].split(':')[0])
        return res

    # given a list of contacts we count how many of each type are present and we assign some penalties, called costs,
    # to make sure that the most common contacts don't overshadow the rest.
    def costants(list_diz_contacts):
        hid, vdw, pip, pic, ioc, ss, iac = 0, 0, 0, 0, 0, 0, 0
        # counting for each type of contact how many of are present
        for i in list_diz_contacts[7]:
            if i[1] == 'HBOND':
                hid = hid + 1
            if i[1] == 'VDW':
                vdw = vdw + 1
            if i[1] == 'PICATION':
                pic = pic + 1
            if i[1] == 'PIPISTACK':
                pip = pip + 1
            if i[1] == 'IONIC':
                ioc = ioc + 1
            if i[1] == 'SSBOND':
                ss = ss + 1
            if i[1] == 'IAC':
                iac = iac + 1
        # total number of contacts
        tot = hid + ss + ioc + pip + pic
        # calculating penalties
        hid = 1 - (hid / tot)
        ss = 1 - (ss / tot)
        ioc = 1 - (ioc / tot)
        pip = 1 - (pip / tot)
        pic = 1 - (pic / tot)

        return hid, ss, ioc, pip, pic

    # Given the index indicating a contact map, the list of contacts and the list of residues, this function
    # builds a matrix with a value in the cell[i][j] if there was a connection between residues i and j, otherwise
    # the cell will contain 0.
    def build_matrix(index, list_contacts, list_res):
        matrix = {}
        # for every contact
        for node1 in list_contacts[index].keys():
            # we calculate the weights based on the number and type of contacts
            costH, costS, costI, costPIP, costPIC = costants(list_contacts)
            matrix.setdefault(node1, [])
            # for every residue, we look if a contact present in the file
            for residue in list_res:
                find = False
                num_contacts = len(list_contacts[index][node1])
                # for every other contact, we search if there is a connection
                for node2 in range(0, num_contacts):
                    # if we find a match, the residue has a contact
                    if residue == list_contacts[index][node1][node2]:
                        find = True
                        # we assign a value based the type of bond using weights and a penalty
                        # and we put the value in the corresponding matrix cell
                        if list_contacts[index][node1][node2][1] == "HBOND":
                            val = 17.0000 * costH
                            matrix[node1].append(val)

                        if list_contacts[index][node1][node2][1] == "VDW":
                            val = 6.000
                            matrix[node1].append(val)

                        if list_contacts[index][node1][node2][1] == "SSBOND":
                            val = 167.000 * costS
                            matrix[node1].append(val)

                        if list_contacts[index][node1][node2][1] == "IONIC":
                            val = 20.000 * costI
                            matrix[node1].append(val)

                        if list_contacts[index][node1][node2][1] == "PIPISTACK":
                            val = 9.400 * costPIP
                            matrix[node1].append(val)

                        if list_contacts[index][node1][node2][1] == "PICATION":
                            val = 9.600 * costPIC
                            matrix[node1].append(val)

                        if list_contacts[index][node1][node2][1] == "IAC":
                            val = 0
                            matrix[node1].append(val)
                # if the contact combination was not present, we put 0  in the corresponding matrix cell
                # to signify that the contact was not present
                if not find:
                    matrix[node1].append(0)
        return matrix

    # Given the total number of files, the contacts list and the list of residues considered, this
    # function returns the list containing the matrices representing the contact maps.
    def list_matrix(num_files, list_contacts, list_res):
        mat_list = []
        # for every file
        for i in range(0, num_files):
            # make matrix
            matrix = build_matrix(i, list_contacts, list_res)
            # add matrix to list
            mat_list.append(matrix)
        return mat_list

    # Function for parsing the file edges paths written in the input file, path is the path to the input file.
    def parseInputFile(path):
        list = []
        names = []
        # open file
        with open(path) as f:
            # for every line, meaning edges file path
            for line in f:
                # ignore comments or empty lines
                if not line.startswith(('#', '\n')):
                    line = line.replace('\n', '')
                    list.append(line)
                    # eliminate spaces, new lines and split values from names
                    names.append(line.split('/')[-1])
        # return all values
        protein_name = list[0].split('/')[-3]
        return protein_name, names, list

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

    # Given two indexes we get the corresponding csv files, compare them and save in a list all the contacts for which
    # they differ from one another.
    # The csv files are temporary files and represent the contact maps as matrices.
    def get_differences(indx_1, indx_2, list_paths, number_residues):
        snap_1 = pd.read_csv(os.path.join(list_paths[indx_1]))
        snap_1 = np.array(snap_1)
        snap_2 = pd.read_csv(os.path.join(list_paths[indx_2]))
        snap_2 = np.array(snap_2)
        residues = snap_1[:, [0]]
        snap_1 = np.delete(snap_1, 0, 1)
        snap_2 = np.delete(snap_2, 0, 1)
        list = []
        for i in range(0, number_residues - 1):
            for j in range(i, number_residues):
                if snap_1[i, j] != snap_2[i, j]:
                    list.append((residues[i][0], residues[j][0]))
        return list

    # Given the indexes of all the clusters centers, we confront them to see in which contacts
    # they differ.
    # This function returns a dictionary of lists that has as key a tuple containing the indexes of two of the centers
    # confronted and as a value a list containing all the residues for which they differ.
    def confront_clusters(centers, list_paths, number_residues):
        elements = len(centers)
        differences = {}

        # comparing all elements avoiding useless and repeated confrontations
        for el_1 in range(0, elements - 1):
            for el_2 in range(el_1 + 1, elements):
                # adding the differences between centers el_1 and el_2 to the dictionary
                differences[(el_1, el_2)] = get_differences(el_1, el_2, list_paths, number_residues)
        return differences

    ############################# ARGUMENT PARSER ############################################################

    # setting up parser
    parser = argparse.ArgumentParser(description='Arguments for clustering contact maps.')
    parser.add_argument('input_file', metavar='contact_maps_file', type=str,
                        help='the file which contains the paths to RING contact map files (edge files), one '
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

    # if input file doesn't exist we show an error
    if not os.path.isfile(args.input_file):
        logging.error('input file does not found, please check path', 'input_file', args.input_file)
        sys.exit(1)

    # if config file isn't found or does not exist we show an error
    if not os.path.isfile(args.config_file):
        logging.error('configuration file not found, please check path', 'config_file', args.config_file)
        sys.exit(1)

    ############################# CONFIG FILE PARSING #########################################################

    # reading all the variables
    try:
        config = configparser.ConfigParser()
        config.read(args.config_file)
        contact_threshold = config['CONTACT MAPS FILTERING CONFIGURATIONS']['contact_threshold']
        energy_threshold = config['CONTACT MAPS FILTERING CONFIGURATIONS']['energy_threshold']
        method = config['DENDROGRAM CONFIGURATIONS']['method']
        optimal_ordering = config['DENDROGRAM CONFIGURATIONS']['optimal_ordering']
        truncation = config['DENDROGRAM CONFIGURATIONS']['truncation']
        truncate_mode = config['DENDROGRAM CONFIGURATIONS']['truncate_mode']
        color_threshold = config['DENDROGRAM CONFIGURATIONS']['color_threshold']
        orientation = config['DENDROGRAM CONFIGURATIONS']['orientation']
        distance_sort = config['DENDROGRAM CONFIGURATIONS']['distance_sort']
        show_leaf_counts = config['DENDROGRAM CONFIGURATIONS']['show_leaf_counts']
        font_size = config['DENDROGRAM CONFIGURATIONS']['font_size']
        damping = config['OPTIMAL CLUSTERING CONFIGURATIONS']['damping']
        random_state = config['OPTIMAL CLUSTERING CONFIGURATIONS']['random_state']
        toll = config['OPTIMAL CLUSTERING CONFIGURATIONS']['toll']
        max_iter = config['OPTIMAL CLUSTERING CONFIGURATIONS']['max_iter']
        delete_temporary = config['OUTPUT CONFIGURATIONS']['delete_temporary']

    except ValueError:
        logging.error("Couldn't parse configuration file, please check", ValueError.args)
        sys.exit(1)

    # checking if values are correct, otherwise show user an error
    try:
        contact_threshold = float(contact_threshold)
    except ValueError:
        logging.error("In config file contact_threshold value is not convertible to float, please check",
                      ValueError.args, 'contact_threshold', contact_threshold)
        sys.exit(1)

    try:
        energy_threshold = float(energy_threshold)
    except ValueError:
        logging.error("In config file energy_threshold value is not convertible to float, please check",
                      ValueError.args, 'energy_threshold', energy_threshold)
        sys.exit(1)

    if method not in ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']:
        logging.error("In config file method value is not correct, please check", 'method', method)
        sys.exit(1)

    try:
        optimal_ordering = bool(strtobool(optimal_ordering))
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

    if distance_sort not in ['ascending', 'descending']:
        try:
            distance_sort = bool(strtobool(distance_sort))
        except ValueError:
            logging.error(
                "In config file distance_sort value is not convertible to bool or is not correct, please check",
                ValueError.args, 'distance_sort', distance_sort)
            sys.exit(1)

    try:
        show_leaf_counts = bool(strtobool(show_leaf_counts))
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
        if damping < 0.5 or damping > 1:
            logging.error("In config file damping value is not correct, it must be between 0.5 and 1, please check",
                          ValueError.args, 'damping', damping)
            sys.exit(1)
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

    try:
        delete_temporary = bool(strtobool(delete_temporary))
    except ValueError:
        logging.error("In config file delete_temporary value is not convertible to bool, please check",
                      ValueError.args, 'delete_temporary', delete_temporary)
        sys.exit(1)

    ############################# READ INPUT ###############################################################

    print('reading contact maps...')

    # the protein name, all files names and all files paths contained in input_file are returned
    # not ordered
    try:
        protein_name, list_names, list_paths = parseInputFile(args.input_file)
    except ValueError:
        logging.error("Couldn't parse input_file, please check", ValueError.args)
        sys.exit(1)

    # the total number of edges files
    number_files = len(list_names)

    # get a list of dictionaries that contain the contacts for all files
    # each cell contains a dictionary that contains all residues contacts in a file
    try:
        list_contacts = list_contacts_per_snapshot(contact_threshold, energy_threshold, list_paths)
    except ValueError:
        logging.error("An error occurred while parsing contacts in edges files, please check", ValueError.args)
        sys.exit(1)

    # get a list that contains all residues contained in a file
    # the list is ordered
    try:
        listRes = list_residues(number_files, list_contacts)
    except ValueError:
        logging.error("An error occurred while parsing residues in edges files, please check", ValueError.args)
        sys.exit(1)

    # total number of residues
    number_residues = len(listRes)

    # matrices rappresent a contact map
    try:
        listMatrix = list_matrix(number_files, list_contacts, listRes)
    except ValueError:
        logging.error("An error occurred while transforming edges files in matrices, please check", ValueError.args)
        sys.exit(1)

    # creating dataframes for each of the matrices representing the snapshots and saving them in the temporary folder
    try:
        for i in range(0, number_files):
            ((pd.DataFrame(listMatrix[i], columns=listRes, index=listRes)).fillna(0)) \
                .to_csv((os.path.join(args.temporary_path, list_names[i] + '.csv')), index=listRes, header=True)
    except ValueError:
        logging.error("An error occurred while saving contact maps as csv files, please check", ValueError.args)
        sys.exit(1)

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

    # insert distances into each element of the distance_matrix, we need a symmetric matrix
    for i in range(0, len(listPaths)):
        f_i = matrices[i].reshape(1, -1)
        for j in range(0, len(listPaths)):
            f_j = matrices[j].reshape(1, -1)
            val = cdist(f_i, f_j, metric='cityblock')
            if (val != 0):
                distance_matrix[i, j] = np.sqrt(val)
            else:
                distance_matrix[i, j] = 0

    # create DataFrame for distance_matrix
    distance_df = pd.DataFrame(distance_matrix, columns=list_names, index=list_names)
    # export DataFrame into file csv
    distance_df.to_csv(os.path.join(args.output_path, protein_name + '_distance_matrix.csv'))

    # making a similarity matrix from distance matrix. The similarity matrix is also called affinity matrix
    # this is a technique used for transforming euclidean distances in a similarity matrix.
    affinity_matrix = distance_matrix * -1

    # indexes of snapshots (0 -> first snapshot, 1-> second snapshot, ...)
    indexes = [i for i in range(number_files)]

    ############################# DENDROGRAM ##################################################################

    # show message to user
    print('making dendrogram...')

    # call function to create dendrogram
    try:
        make_dendrogram(squareform(distance_matrix), indexes, args.output_path, method, optimal_ordering, orientation,
                        distance_sort,
                        show_leaf_counts, truncation, truncate_mode, color_threshold, protein_name, font_size)
    except ValueError:
        logging.error("An error occurred while making the dendrogram, please check", ValueError.args)
        sys.exit(1)

    ############################# CLUSTERING ############################################################

    # message to user
    print('clustering...')

    try:
        # creating instance of clustering class
        clustering = affinity_propagation(indexes, damping, max_iter, toll, random_state, number_residues,
                                      matrices, affinity_matrix)

        # running clustering algorithm and getting results
        final_clusters, final_centers = clustering.clusterize()
    except ValueError:
        logging.error("An error occurred while clustering with Affinity Propagation, please check", ValueError.args)
        sys.exit(1)

    ############################# PROCESSING RESULTS ############################################################

    print('working on outputs...')

    # show an example for each cluster. In this case, given the nature of Affinity Propagation, the examples will be
    # the clusters centers.
    for index in final_clusters:
        sample = final_centers[index]
        shutil.copy(list_paths[sample], os.path.join(args.output_path, protein_name + '_Cluster_{}_'
                                                    .format(index) + 'extract_'+list_names[sample]))

    # create output file with the clustering results
    file = open(os.path.join(args.output_path, protein_name + '_cluster_results.txt'), 'w+')
    file.write('CLUSTERING RESULTS \n')
    # dictionary indices representing clusters
    for index in final_clusters:
        file.write('CLUSTER {}:\n'.format(index))
        # elements in dictionary
        for cluster in final_clusters[index]:
            # write element and distance from center
            file.write('Element: ' + str(cluster) + ', distance from center: ' +
                       str(distance_matrix[cluster][final_centers[index]]) + '\n')
    file.close()

    # output file with clusters differences
    file = open(os.path.join(args.output_path, protein_name + '_cluster_differences.txt'), 'w+')
    file.write('CLUSTERS DIFFERENCES \n')
    differences = confront_clusters(final_centers, listPaths, number_residues)
    # dictionary indices representing the centers confronted
    for index in differences:
        file.write('CLUSTER {} - CLUSTER {}:\n'.format(index[0], index[1]))
        # for all contacts in dictionary
        for contacts in differences[index]:
            file.write(str(contacts) + '\n')
    file.close()

    # output file with legend for relation between indexes and snapshots
    file = open(os.path.join(args.output_path, protein_name + '_legend.txt'), 'w+')
    file.write('LEGEND: \n')
    for index in indexes:
        file.write('Snapshot {} -> index {}\n'.format(list_names[index], index))
    file.close()

    # check if the temporary folder has to be deleted
    if delete_temporary:
        # delete temporary folder and all of it's content
        shutil.rmtree(args.temporary_path)

    # terminate program
    sys.exit(0)


if __name__ == "__main__":
    main()
