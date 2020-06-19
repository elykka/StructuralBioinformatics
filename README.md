# USER MANUAL

**SCOPE**

This program wants to test a different approach to try and provide a clusterization of the structures obtained in Molecular Dynamics trajectories, this is done using the structure's contact map, one for each of them.
The reason behind the use of contact maps is because they provide more informations about the type of contacts that each residue in the structure has.
This program uses the contact maps returned by RING, a Residue Interaction Network, that saves them in files called edges.
The lists will also be referred as snapshots in this document.

**Please note that this program is just an experimental test, assigned during the course of Structural Bioinformatics offered by the University of Padova.
The scope of this project was to choose, motivate, test and evaluate a different approach for the clustering of structures and is not to be viewed as a finished and valid product.**
 
**USE**

To launch the program, the following instruction must be used in the terminal:

python3 	nome_file.py 	input_file_path	 output_directory_path 	config_file_path 	temporary_directory_path

Please note that the directories paths must be in this format *'/home/path/output_folder'*, devoid of a final backslash, and that the files paths must contain the full file name, extension included if present.

At the end of the program the temporary directory and all it's files will be deleted.

**If the output or the temporary directory don't already exist, they will be created.**

**OUTPUTS**

The program will produce some files regarding a general clusterization, such as a dendrogram and a distance matrix, but also some others regarding what could be an optimal clusterization for all provided snapshots.
The metric used in the program was custom made by the team, whereas the clustering algorithm used to try and find an optimal clusterization was the Affinity Propagation algorithm.
**More informations about the metric used and how Affinity propagation works can be found in the report**

After the program has ended, in the output directory the user can find:
- a file containing the dendrogram, the snapshots are represented as indexes, see legend for the association.
- a file containing the distance matrix between all snapshots, calculated by using the custom defined metric.
- a file containing the results of the optimal clusterization, indicating the clusters, which elements belong to which clusters and their distance from the cluster center. The elements are represented as indexes, see legend to refer back to the snapshot associated.
- a file containing the differences from cluster to cluster, in terms of residues and contacts, all the contacts that were considered in the program are shown. The differences were found by confronting the clusters' centers since Affinity Propagation tries to identify the center of the cluster as the element that is more suited to represent all other elements in the cluster.
- an example from each cluster in the optimal clusterization, the example chosen were the centers of the clusterization.
- a file containing a legend to show the relation between indexes and snapshots.

**CONFIGURATION FILE**
The program uses a configuration file to give the user the possibility of changing the algorithm parameters.
The parameters contained are:

*CONTACT MAPS FILTERING CONFIGURATIONS*

- contact_treshold: residues that have a bond with a distance above the threshold are not considered
- energy_treshold: residues with an energy bond below the threshold are not considered
