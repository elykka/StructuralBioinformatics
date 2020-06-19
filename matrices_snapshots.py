import pandas as pd
import os,os.path
from os import listdir
import re


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
    hid=1-(hid/tot)
    ss = 1-(ss/tot)
    ioc =1-(ioc/tot)
    pip =1-(pip/tot)
    pic =1-(pic/tot)

    return hid,ss,ioc,pip,pic


def build_matrix(i,listDizCont,listRes):#Per ogni file costruisco una matrice
    matrix= {}
    for node1 in listDizCont[i].keys():   # per tutti i residui di un file, metto 1 per i contatti e metto 0 per gli altr
        costH, costS,costI,costPIP,costPIC = costants(node1)
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


def list_matrix(numFiles,listDizCont,listRes):
    listMatrix= []
    for i in range(0, numFiles):  # per ogni file,
        matrix = build_matrix(i,listDizCont,listRes)
        listMatrix.append(matrix)
    return listMatrix

###################################
#########   CODICE  ##############
###################################

dir='FileProtein/cdk6_p16ink4a/edges/'


listfiles = listdir(dir)

listfiles.sort(key=lambda f: int(re.sub('\D', '', f)))
print(listfiles)

print(len(listfiles))
listPaths = [dir + listfiles[i] for i in range(len(listfiles))]
#list= qua ho in una lista tutti i cammini
numFiles=len(listPaths)
print(numFiles)
#ottengo una lista di dizionari di contatti per tutti i file
#ogni cella contiene un dizionario che descrive tutti i contatti dei residui in un file
contact_treshold =5
energy_treshold =7

listDizContacs=list_contacts_per_snapshot(contact_treshold,energy_treshold,listPaths)

print(len(listDizContacs))

#ottengo la lista che contieni tutti i residui presenti in tutti i file
#lista di residui ordinata
listRes=list_residues(numFiles,listDizContacs)
print(len(listRes))
#ogni elemento rappresenta la matrice dei contatti di un file
listMatrix=list_matrix(numFiles,listDizContacs,listRes)


for i in range(0,numFiles):
    print(listMatrix[i])



#construisco la lista di dataframe

#dfListMatrix = []
#for i in range(0,numFiles):
#    dfListMatrix.append((pd.DataFrame(listMatrix[i],columns=listRes,index=listRes)).fillna(0))

#salvo tutti i DF in una cartella
#for i in range(0,numFiles):
#    dfListMatrix[i].to_csv(("DataFrame/{}.csv".format(i+1)),index=listRes,header=True)

