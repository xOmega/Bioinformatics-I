import numpy as np

def iteration(A,labels,iter):

    # Calculate divergergence of each OTU from all other OTU's
    # Sum of all the columns
    diverg = np.sum(A,axis=1) 
    n = A.shape[0] 

    if n == 2:
        d = A[1][0]
        nod = node_dictionary[labels[0]]
        nod[labels[1]] = d
        return None,labels
 
    i=0
    j=0
    min_distance = A[i][j]
    for r,row in enumerate(A):
        if r == 0:
        	continue
        for c,col in enumerate(row):
            if c >= r:  
            	continue

            Mrc = A[r][c] - ((diverg[c] + diverg[r])/(n-2))
            if Mrc < min_distance:
                i,j,min_distance = r,c,Mrc
                
  
    print 'Cluster together {', labels[i], labels[j], '} to', iter
   
   
    d_i = A[i][j]/2.0 + (diverg[i] - diverg[j])/(2*(n-2))
    d_j = A[i][j] - d_i
    
    node = { labels[i] : d_i,
             labels[j] : d_j } 
    node_dictionary[iter] = node
    
    t = []
    ij_d = A[i][j]
    for k in range(len(A[0])):
        if k == i or k == j:  continue
        d = (A[i][k] + A[j][k] - A[i][j])/2
        t.append(d)
	   
    
    dele = range(n)
    for k in [j,i]:
        dele.remove(k)
        A1 = A[dele,:]
        A2 = A1[:,dele]

    A = A2
 
    labels = [iter] + labels[:j] + labels[j+1:i] + labels[i+1:]

    new_col = np.array(t)
    new_col.shape = (n-2,1)
    A = np.hstack([new_col,A])
    new_row = np.array([0] + t)
    new_row.shape = (1,n-1)
    A = np.vstack([new_row,A])

    return A,labels
   


node_dictionary = {}

chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

filename = "dataset.txt"
FH = open(filename,'r')
data = FH.read().strip()
FH.close()
data = data.split('\n')

A = []
for i in data:
	A.append([float(n) for n in i.split()])
labels = list(chars[:len(data)]) # [A,B,C,D,E,F,G]

A = np.array(A)

iter=0

while A is not None:
	print '\n############################Iteration : ', iter, ' #########################################\n'
	A,labels = iteration(A,labels,iter)
	iter+=1


print "\nFinal tree \n"

for i in sorted(node_dictionary.keys()):
    print i, ':  '
    node = node_dictionary[i]
    print node


