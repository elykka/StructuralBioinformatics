import pandas as pd
from os import listdir
import re
import numpy as np


def distanceJaccard(A,B):
    M11,M01_M10=0,0
    #index row: i
    for i in range(1,A.shape[0]):
        #index column: j
        for j in range(0,i):
            # if there is a match between the 2 contact maps, (1,1)
            if A[i,j]+B[i,j]==2:
                M11+=1
            #if there is a mismatch, (0,1) or (1,0)
            elif (A[i,j]==0 and B[i,j]==1) or (A[i,j]==1 and B[i,j]==0):
                M01_M10+=1
    # return (a+b)/(a+b+c)
    return (M01_M10)/(M01_M10+M11)



#listafile
dir='DataFrame/'

listfiles = listdir(dir)
print(listfiles)
listfiles.sort(key=lambda f: int(re.sub('\D', '', f)))
print(listfiles)
listPaths = [dir + listfiles[i] for i in range(len(listfiles))]
#list= qua ho in una lista tutti i cammini

#print(listPaths)
l=[]  #ho una lista che contiene tutte le la matrici
for path in listPaths:
    m=pd.read_csv(path,index_col=0)
    m= m.to_numpy() #TRASFORMO GIA' IN NP ARRAY
    l.append(m)


#NXN matrix with all zeros
dist=np.zeros((len(listPaths),len(listPaths)),dtype=float)

from scipy.spatial import distance
#insert Jaccard distance into each element of the matrix dist
for i in range(1,len(listPaths)):
    f_i = l[i].flatten()
    for j in range(0,i):
        f_j= l[j].flatten()
        dist[i,j]=distance.jaccard(f_i,f_j)

#print(dist)

#create DataFrame
distance_df=pd.DataFrame(dist,columns=listfiles,index=listfiles)
#export DataFrame into file csv
distance_df.to_csv('dist_mat_ant.csv')