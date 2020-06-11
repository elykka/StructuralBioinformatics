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

first_matrix= pd.read_csv("dataframe_a/6J6Y_1_ms_1K0.1.csv",index_col=0)
second_matrix= pd.read_csv("dataframe_a/6J6Y_1_ms_1K0.2.csv",index_col=0)

#Trasformo in array
first_matrix= pd.DataFrame.to_numpy(first_matrix)
second_matrix= pd.DataFrame.to_numpy(second_matrix)


dis= distanceJaccard(first_matrix,second_matrix)
print(dis)


#listafile
dir='dataframe_a/'
listfiles = listdir(dir)
listfiles.sort(key=lambda f: int(re.sub('\D', '', f)))
listPaths = [dir + listfiles[i] for i in range(len(listfiles))]

print(listPaths)
l=[]
for path in listPaths:
    m=pd.read_csv(path,index_col=0)
    l.append(m.values)

#NXN matrix with all zeros
dist=np.zeros((len(listPaths),len(listPaths)),dtype=float)

#insert Jaccard distance into each element of the matrix dist
for i in range(1,len(listPaths)):
    for j in range(0,i):
        dist[i,j]=distanceJaccard(l[i],l[j])

print(dist)

#create DataFrame
distance_df=pd.DataFrame(dist,columns=listfiles,index=listfiles)
#export DataFrame into file csv
distance_df.to_csv('dist_mat/dist_mat_ant.csv')