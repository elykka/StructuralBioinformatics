def list_contacts_per_snapshot(contact_threshold,energy_treshold,listPaths):
        listCon = []  # ogni elemento contiene il dizionario dei contatti di un file
        numFiles = len(listPaths)
        for i in range(0, numFiles):
            contacts = {}
            with open(listPaths[i]) as f:
                next(f)
                for line in f:
                    node_1, edge, node_2, distance, _,energy, atom_1, atom_2 = line.split()[0:8]
                    distance = float(distance)
                    energy = float(energy)
                    if (distance <= contact_threshold) & (energy > energy_treshold):
                        edge_type, edge_loc = edge.split(":")  # edge_loc values are LIG, MC, SC
                        contacts.setdefault((node_1,edge_type), [])  # for both node_1 end node_2 I have a key in the contacts
                        contacts.setdefault((node_2,edge_type), [])
                        contacts[(node_1,edge_type)].append((node_2, edge_type))
                        contacts[(node_2,edge_type)].append((node_1, edge_type))
                listCon.append(contacts)
        return listCon

def list_residues(numFiles,listDizCont):
    res = []  # contiene tutti i residui di tutti i filenella lista
    for i in range(0, numFiles):
        for node in listDizCont[i].keys():
            numNode2= len(listDizCont[i][node])
            for j in range(0,numNode2) :
                if listDizCont[i][node][j] not in res:
                    res.append(listDizCont[i][node][j])
    res = list(dict.fromkeys(res))  # elimino i duplicati
    res.sort(key=lambda x: int(x[0].split(':')[1]))  # ordine numerico
    res.sort(key=lambda x: x[0].split(':')[0])  # ordine alfabetico
    return res

def costants(node1, listDizContacs):
    hid, vdw ,pip ,pic ,ioc ,ss ,iac = 0,0,0,0,0,0,0
    for i in listDizContacs[7]:
        if (i[1] == 'HBOND'):
            hid = hid + 1
        if (i[1] == 'VDW'):
            vdw = vdw + 1
        if (i[1] == 'PICATION'):
            pic = pic + 1
        if (i[1] == 'PIPISTACK'):
            pip = pip + 1
        if (i[1] == 'IONIC'):
            ioc = ioc + 1
        if (i[1] == 'SSBOND'):
            ss = ss + 1
        if (i[1] == 'IAC'):
            iac = iac + 1
    tot = hid + ss + ioc + pip +pic
    hid=1-(hid/tot)
    ss = 1-(ss/tot)
    ioc =1-(ioc/tot)
    pip =1-(pip/tot)
    pic =1-(pic/tot)

    return hid,ss,ioc,pip,pic


def build_matrix(i, listDizCont, listRes):#Per ogni file costruisco una matrice
    matrix= {}
    for node1 in listDizCont[i].keys():   # per tutti i residui di un file, metto 1 per i contatti e metto 0 per gli altr
        costH, costS,costI,costPIP,costPIC = costants(node1, listDizCont)
        matrix.setdefault(node1, [])
        for re in listRes:
            #guardo se è presente
            find = False
            numCon = len(listDizCont[i][node1])
            for node2 in range(0, numCon):
                if re == listDizCont[i][node1][node2]:
                    find = True
                    if listDizCont[i][node1][node2][1] == "HBOND":
                        val = 17.0000*costH
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "VDW":
                        val = 6.000
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "SSBOND":
                        #* 35.5
                        val = 167.000*costS
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "IONIC":
                        #* 12.69
                        val = 20.000*costI
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "PIPISTACK":
                        #* 88.75
                        val = 9.400*costPIP
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "PICATION":
                        #* 88.75
                        val = 9.600*costPIC
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "IAC":
                        val = 0
                        matrix[node1].append(val)
            if find == False:
                matrix[node1].append(0)
    return matrix


def list_matrix(numFiles, listDizCont, listRes):
    listMatrix= []
    for i in range(0, numFiles):  # per ogni file,
        matrix = build_matrix(i, listDizCont, listRes)
        listMatrix.append(matrix)
    return listMatrix
