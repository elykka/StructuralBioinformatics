import pandas as pd
import os,os.path
from os import listdir
import re


def list_contacts_per_snapshot(contact_threshold,energy_treshold,listPaths):
        listCon = []
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
                        edge_type, edge_loc = edge.split(":")
                        contacts.setdefault((node_1,edge_type), [])
                        contacts.setdefault((node_2,edge_type), [])

                        if (node_2,edge_type) not in contacts[(node_1,edge_type)]:
                            contacts[(node_1,edge_type)].append((node_2, edge_type))

                        if (node_1, edge_type) not in contacts[(node_2, edge_type)]:
                            contacts[(node_2,edge_type)].append((node_1, edge_type))

                listCon.append(contacts)
        return listCon

def list_residues(numFiles,listDizCont):
    res = []
    for i in range(0, numFiles):
        for node in listDizCont[i].keys():
            if node not in res:
                res.append(node)
    res = list(dict.fromkeys(res))  # delete duplicate
    res.sort(key=lambda x: int(x[0].split(':')[1]))  # ordine numerico
    res.sort(key=lambda x: x[0].split(':')[0])  # alphabetical order
    return res

def costants(node1):
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
    #print(hid, ss, ioc, pip, pic)

    hid=1-(hid/tot)
    ss = 1-(ss/tot)
    ioc =1-(ioc/tot)
    pip =1-(pip/tot)
    pic =1-(pic/tot)

    return hid,ss,ioc,pip,pic


def build_matrix(i,listDizCont,listRes):#Per ogni file costruisco una matrice
    matrix= {}
    for node1 in listDizCont[i].keys():#node1 is a residue of the snapshot
        costH, costS,costI,costPIP,costPIC = costants(node1)


        matrix.setdefault((node1), []) #per ogni residuo devo crearmi una lista che identifica i contatti
        for re in listRes:

            find = False
            numCon = len(listDizCont[i][node1])

            for node2 in range(0, numCon):#node2 is one of the contacts of node1
                if re == listDizCont[i][node1][node2]:
                    find = True
                    if listDizCont[i][node1][node2][1] == "HBOND":
                        val = 17.0000*costH
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "VDW":
                        val = 6.000
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "SSBOND":
                        val = 167.000*costS
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "IONIC":
                        val = 20.000*costI
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "PIPISTACK":
                        val = 9.400*costPIP
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "PICATION":
                        val = 9.600*costPIC
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "IAC":
                        val = 0
                        matrix[node1].append(val)
                    if listDizCont[i][node1][node2][1] == "IAC":
                        val = 0
                        matrix[node1].append(val)
            if find == False:
                matrix[node1].append(0)
    return matrix


def list_matrix(numFiles,listDizCont,listRes):
    listMatrix= []
    for i in range(0, numFiles):
        matrix = build_matrix(i,listDizCont,listRes)
        listMatrix.append(matrix)

    return listMatrix

###################################
#########   CODICE  ##############
###################################

dir='FileProtein/frataxin/edges/'

listfiles = listdir(dir)

listfiles.sort(key=lambda f: int(re.sub('\D', '', f)))

listPaths = [dir + listfiles[i] for i in range(len(listfiles))]
numFiles=len(listPaths)


contact_treshold =5
energy_treshold =7

#each dictionary contains the contacts of each snapshot
listDizContacs=list_contacts_per_snapshot(contact_treshold,energy_treshold,listPaths)


#list of all the residues of the protein
listRes=list_residues(numFiles,listDizContacs)


#each matrix rapresent the contact map of each snapshot
listMatrix=list_matrix(numFiles,listDizContacs,listRes)

dfListMatrix = []
for i in range(0,numFiles):
    dfListMatrix.append((pd.DataFrame(listMatrix[i],columns=listRes,index=listRes)).fillna(0))

#save
for i in range(0,numFiles):
    dfListMatrix[i].to_csv(("DataFrame/{}.csv".format(i)),index=listRes,header=True)

