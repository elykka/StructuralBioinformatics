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
