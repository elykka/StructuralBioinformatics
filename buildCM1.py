import pandas as pd
import os,os.path


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
        residuiFile= []
        for node in listDizCont[i].keys():
            numNode2= len(listDizCont[i][node])
            for j in range(0,numNode2) :
                if listDizCont[i][node][j] not in res:
                    res.append(listDizCont[i][node][j])
    res = list(dict.fromkeys(res))  # elimino i duplicati
    res.sort(key=lambda x: int(x[0].split(':')[1]))  # ordine numerico
    res.sort(key=lambda x: x[0].split(':')[0])  # ordine alfabetico
    return res

def build_matrix(i,listDizCont,listRes):#Per ogni file costruisco una matrice
    matrix= {}
    for node1 in listDizCont[i].keys():   # per tutti i residui di un file, metto 1 per i contatti e metto 0 per gli altr
        matrix.setdefault(node1, [])
        for re in listRes:
            #guardo se è presente
            find = False
            numCon = len(listDizCont[i][node1])
            for node2 in range(0, numCon):
                if re == listDizCont[i][node1][node2]:
                    find = True
                    if (listDizCont[i][node1][node2][1] == "HBOND" or
                        listDizCont[i][node1][node2][1] == "SSBOND" or
                        listDizCont[i][node1][node2][1] == "IONIC" or
                        listDizCont[i][node1][node2][1] == "PIPISTACK"or
                        listDizCont[i][node1][node2][1] == "PICATION"or
                        listDizCont[i][node1][node2][1] == "IAC"):
                        matrix[node1].append(1)
            if find == False:
                matrix[node1].append(0)
    return matrix


def list_matrix(numFiles,listDizCont,listRes):
    listMatrix= []
    for i in range(0, numFiles):  # per ogni file,
        matrix = build_matrix(i,listDizCont,listRes)
        listMatrix.append(matrix)
    return listMatrix

###################################
#########   CODICE  ##############
###################################

dir='antibody/edges/'

#ottengo il numero di file nella cartella edges
numFiles=len([name for name in os.listdir(dir) if os.path.isfile(os.path.join(dir, name))])
print(numFiles)

#inserisco tutti i path in una lista
listPaths = []
for i in range(0,numFiles):
    listPaths.append("antibody/edges/6J6Y_1_ms_1K0.{}.pdb_ringEdges".format(i+1))

#ottengo una lista di dizionari di contatti per tutti i file
#ogni cella contiene un dizionario che descrive tutti i contatti dei residui in un file
contact_treshold =5
energy_treshold =7
listDizContacs=list_contacts_per_snapshot(contact_treshold,energy_treshold,listPaths)
print(listDizContacs[0])

#ottengo la lista che contieni tutti i residui presenti in tutti i file
#lista di residui ordinata
listRes=list_residues(numFiles,listDizContacs)
print(len(listRes))
#ogni elemento rappresenta la matrice dei contatti di un file
listMatrix=list_matrix(numFiles,listDizContacs,listRes)

#print(listMatrix[0]['A:152:_:GLN', 'HBOND'])
print(len(listMatrix))

#construisco la lista di dataframe
dfListMatrix = []
for i in range(0,numFiles):
    dfListMatrix.append((pd.DataFrame(listMatrix[i],columns=listRes,index=listRes)).fillna(0))

#salvo tutti i DF in una cartella
for i in range(0,numFiles):
    dfListMatrix[i].to_csv(("dataframe_a/6J6Y_1_ms_1K0.{}.csv".format(i+1)),index=listRes,header=True)



