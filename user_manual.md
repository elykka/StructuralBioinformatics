# USER MANUAL

### PURPOSE

snapshotClusterer wants to test a different approach to try and provide a clusterization of the structures obtained in Molecular Dynamics trajectories, this is done using the structure's contact map, one for each of them.

The reason behind the use of contact maps is because they provide more informations about the type of contacts that each residue in the structure has.

This program uses the contact maps returned by RING, a Residue Interaction Network, that saves them in files called edges.
This maps will also be referred as snapshots in this document.

**Please note that snapshotClusterer is just an experimental test, assigned during the course of Structural Bioinformatics offered by the University of Padova.
The scope of this project was to choose, motivate, test and evaluate a different approach for the clustering of structures and is not to be viewed as a finished and valid product.**
 
### USE

To launch the program, the following instruction must be used in the terminal:

`python3  	snapshotClusterer.py  	input_file_path.txt	  output_directory_path  	config_file_path/config.ini 	 temporary_directory_path`

where:
* *input_file_path.txt* is the path to the file txt that contains all of the paths to the edges files of a protein that are to be clustered. An example of the file is shown below.
![Example of input file](/img/input_example.png)

* *output_directory_path* is the path to the folder where the results are to be saved.
* *config_file_path/config.ini* is the path to the file ini that contains all of the parameters that are to be used in the program and is provided by the team.
* *temporary_directory_path* is the path to the folder where the temporary files are to be saved.

Please note that the directories paths must be in this format *'/home/path/output_folder'*, devoid of a final backslash, and that the files paths must contain the full file name, extension included if present.

**If the output or the temporary directory don't already exist, they will be created, however the rest of the path leading up to the directory must exist.**

### OUTPUTS

The program will produce some files regarding a general clusterization, such as a dendrogram and a distance matrix, but also some others regarding what could be an optimal clusterization for all provided snapshots.
The metric used in the program was custom made by the team, whereas the clustering algorithm used to try and find an optimal clusterization was the Affinity Propagation algorithm.

**More informations about the metric used and how Affinity propagation works can be found in the report**

After the program has ended, in the output directory the user can find:
* *__protein-name__\_legend.txt*, a file containing a legend to show the relation between indexes and snapshots, since in the dendrogram and in some of the output files the snapshots were represented by indexes.
* *__protein-name__\_dendrogram.png*, a file containing the snapshots dendrogram calculated using the metric defined by the team.
  * NOTE: the elements are represented as indexes, see legend to refer back to the snapshot associated.
* *__protein-name__\_distance_matrix.csv*, a file containing the distance matrix between all snapshots, calculated by using the metric defined by the team.
* *__protein-name__\_cluster_results.txt*, a file containing the results of the optimal clusterization, indicating:
  * the clusters;
  * which elements belong to which cluster;
  * their distance from the cluster center;
  * NOTE: the elements are represented as indexes, see legend to refer back to the snapshot associated.
* *__protein-name__\_cluster_differences*, a file containing the differences from cluster to cluster, in terms of residues and contacts. All of the contacts that were considered during the clustering are shown. 
  * the differences were found by confronting the clusters centers since Affinity Propagation tries to identify the center of the cluster as the element that is more suited to represent all other elements in the cluster.
  * the files contains an indication of which clusters were confronted, "CLUSTER i - CLUSTER j", followed by the list of contacts for which the two clusters differ.
  * the file does not contain repetitions, if the user wants to find the differences between cluster i and j and i < j, the file will contain the comparison under "CLUSTER i - CLUSTER j". 
* a series of edges files used as an example, one for each cluster in the optimal clusterization, the examples chosen were the centers of the clusterization.
  * the file name are going to be *__protein-name__\_Cluster_**x**\_extract_**edges-file-name***, where **x** is the number of cluster and **edges-file-name** is the name of the edges file used as an example.

**protein-name is obtained from the paths of edges files, usually the path is .../protein-name/edges/edges-file-name should that folder name be different then the resulting output files will have that name**

### CONFIGURATION FILE
The program uses a configuration file to give the user the possibility of changing the algorithm parameters.
The parameters contained are the following:

**CONTACT MAPS FILTERING CONFIGURATIONS**

* contact_treshold: residues that have a bond with a distance above the threshold are not considered, the distance considered is the one reported in the edges files.
  * acceptable values are floats (with a .).
* energy_treshold: residues with an energy bond below the threshold are not considered.
  * acceptable values are floats (with a .), the energy bonds are all floats ranging from 5 to 167 and change depending on the type of bond:
    * Hidrogen bonds (HBOND) have energy 17.000;
    * Van der Waals interactions (VDW) have energy 6.000;
    * Disulfide bridges (SBOND) have energy 167.000;
    * Salt bridges (IONIC) have energy 20.000;
    * Pipistacking (PIPISTACK) have energy 9.400;
    * Pication (PICATION) have energy 9.600;
    * Inter atomic contacts (IAC, meaning that two residues are close but we were not able to identify the type of bond) have energy 0.

**DENDROGRAM CONFIGURATIONS**

* method: type of linkage.
  * acceptable values: single, complete, average, weighted, centroid, median, ward.
* optimal_ordering: if True, the linkage matrix will be reordered so that the distance between successive leaves is minimal, it might become computationally expensive.
  * acceptable values are: True, False.
* truncation: value for truncate mode.
  * acceptable values: int.
* truncate_mode: condenses the dendrogram if it's difficult to read.
  * acceptable values: 
    * 'none', meaning no truncation;
    * 'lastp', meaning that the last 'truncation' non-singleton clusters formed in the linkage are going to be the only non-leaf nodes in the linkage;
    * 'level', meaning that no more than 'truncation' levels of the dendrogram are displayed. A “level” includes all nodes with 'truncation' merges from the last merge.
* color_threshold: colors all the descendent links below a cluster node the same color if is the first node below the threshold.   
  * acceptable values: float (with .).
* orientation: the direction in which to plot the dendrogram.
  * acceptable values: top, bottom, left, right.
* distance_sort: for each node n, the order (visually, from left-to-right) n’s two descendent links are plotted is determined by this parameter.
  * acceptable values: False, ascending or True, descending.
* show_leaf_counts: when True, leaf nodes representing original observation are labeled with the number of observations they contain in parentheses.
  * acceptable values: True or False.
* font_size: change dendrogram's labels size.
  * acceptable values: int.

**OPTIMAL CLUSTERING CONFIGURATIONS**

* damping: is the extent to which the current value is maintained relative to incoming values (weighted 1 - damping). This in order to avoid numerical oscillations.
  * acceptable values: float (with .), between 0.5 and 1.
* random_state: generator to control the starting state of the algorithm, an int is used for reproducible results across function calls.
  * acceptable values: int.
* toll: number of iterations with no change in the number of estimated clusters that stops the convergence.
  * acceptable values: int.
* max_iter: maximum number of iterations.
  * acceptable values: int.
  
**OUTPUT CONFIGURATIONS**

* delete_temporary: if True, at the end of the program the temporary directory and all it's files will be deleted.
  * acceptable values: True or False.

**In the program we have used scipy dendrogram implementation and scikit-learn Affinity Propagation implementation, for doubts or curiosities on paramenters use, refer to the official documentation:**
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html
- https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AffinityPropagation.html
