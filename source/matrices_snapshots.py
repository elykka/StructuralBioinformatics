# Given the energy threshold, the contact threshold and the list of all file edges, this function
# returns the list of contacts that are to be considered for all snapshots.
# The list returned is actually a list of dictionaries, where the key represents the residue, node_1,that we found
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
# residues contained in the files
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

# given a list of contacts we count how many of each type are present and we assign some weights,
# to make sure that the most common contacts don't overshadow the rest
def weights(list_diz_contacts):
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
    # total number of weights
    tot = hid + ss + ioc + pip + pic
    # calculating penalties
    w_hid = 1 - (hid / tot)
    w_ss = 1 - (ss / tot)
    w_ioc = 1 - (ioc / tot)
    w_pip = 1 - (pip / tot)
    w_pic = 1 - (pic / tot)

    return w_hid, w_ss, w_ioc, w_pip, w_pic

# Given the index indicating a contact map, the list of contacts and the list of residues, this function
# builds a matrix with a value in the cell[i][j] if there was a connection between residues i and j, otherwise
# the cell will contain 0.
def build_matrix(index, list_contacts, list_res):
    matrix = {}
    # for every contact
    for node1 in list_contacts[index].keys():
        # we calculate the weights based on the number and type of contacts
        wH, wS, wI, wPIP, wPIC = weights(list_contacts)
        matrix.setdefault(node1, [])
        # for every residue, we look if it is present in the list of contacts of the node1
        for residue in list_res:
            find = False
            num_contacts = len(list_contacts[index][node1])
            # for every contact of node1, we search if the residue is equal to node2
            for node2 in range(0, num_contacts):
                # if we find a match, the residue is in the list of contacts of node1
                if residue == list_contacts[index][node1][node2]:
                    find = True
                    # we assign a value based on the energy bound and the weight
                    # and we put the value in the corresponding matrix cell
                    if list_contacts[index][node1][node2][1] == "HBOND":
                        val = 17.0000 * wH
                        matrix[node1].append(val)

                    if list_contacts[index][node1][node2][1] == "VDW":
                        val = 6.000
                        matrix[node1].append(val)

                    if list_contacts[index][node1][node2][1] == "SSBOND":
                        val = 167.000 * wS
                        matrix[node1].append(val)

                    if list_contacts[index][node1][node2][1] == "IONIC":
                        val = 20.000 * wI
                        matrix[node1].append(val)

                    if list_contacts[index][node1][node2][1] == "PIPISTACK":
                        val = 9.400 * wPIP
                        matrix[node1].append(val)

                    if list_contacts[index][node1][node2][1] == "PICATION":
                        val = 9.600 * wPIC
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
