
import pandas as pd
from os import listdir
import re
import numpy as np
from scipy.spatial.distance import cdist


#listafile
dir='DataFrame/'

listfiles = listdir(dir)
print(listfiles)
listfiles.sort(key=lambda f: int(re.sub('\D', '', f)))
print(listfiles)
listPaths = [dir + listfiles[i] for i in range(len(listfiles))]


#print(listPaths)
l=[]
for path in listPaths:
    m=pd.read_csv(path,index_col=0)
    m= m.to_numpy()
    l.append(m)


#NXN matrix with all zeros
dist=np.zeros((len(listPaths),len(listPaths)),dtype=float)


for i in range(0,len(listPaths)):
    f_i = l[i].reshape(1,-1)
    for j in range(0,len(listPaths)):
        f_j = l[j].reshape(1, -1)
        val = cdist(f_i, f_j, metric='cityblock')
        if (val != 0):
            dist[i,j] = np.sqrt(val)
        else:
            dist[i, j] = 0



#create DataFrame
distance_df=pd.DataFrame(dist,columns=listfiles,index=listfiles)
#export DataFrame into file csv
distance_df.to_csv('dist_mat_MANH.csv')