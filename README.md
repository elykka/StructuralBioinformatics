# USER MANUAL

**SCOPE**
This program uses the list of contacts in a protein complex returned by RING, a Residue Interaction Network.
The program wants to provide a clusterization of the structures in Molecular Dynamics trajectories, using 
RING list of contacts, one per structure.
This lists will also be referred as snapshots in this document.
 
**USE**
To launch the program, the following instruction must be used in the terminal:

python3 nome_file.py  input_file_path  output_directory_path  config_file_path  temporary_directory_path

Please note that the directories paths must be in this format *'/home/path/output_folder'*, devoid of a final backslash, and that the files paths must contain the full file name, extension included.

If the output or the temporary directory doesn't already exist, it will be created.

**OUTPUTS**
The program will provide informations about a general clusterization, such as a dendrogram and a distance matrix, but also some regarding what he believes to be the optimal number of clusters for all the snapshots provided.

After the program has ended, in the output directory the user can find:
- a file containing the dendrogram
- a file containing the distance matrix
- a file containing the results of the optimal clusterization
- a file containing the differences from cluster to cluster, in terms of residues and contacts.
- an example from each cluster in the optimal clusterization
