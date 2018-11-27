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
    print min_distance
   
   
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
filename = "distances_jukes_cantour.txt"
FH = open(filename,'r')
data = FH.read().strip()
FH.close()
data = data.split('\n')

A = []
for i in data:
	A.append([float(n) for n in i.split()])

print A


'''
A = [
 [0,          0.07361198, 0.05302558, 0.13038054, 0.06331878, 0.03867748,  0.13038054, 0.0205864 ],
 [0.02183406, 0,          0.05583281, 0.12944479, 0.05988771, 0.0371179,  0.12944479, 0.02183406],
 [0.02432938, 0.0658141,  0,          0.02432938, 0.06363069, 0.0371179,  0.13131628, 0.02432938],
 [0.02276981, 0.0679975,  0.0555209,  0,          0.05895197, 0.04148472, 0.12476606, 0.02276981],
 [0.02121023, 0.07049283, 0.05489707, 0.13006862, 0,          0.0371179,  0.13006862, 0.02121023],
 [0.02089832, 0.06924517, 0.05396132, 0.13100437, 0.06113537, 0,         0.13100437, 0.02089832],
 [0.02089832, 0.08328135, 0.05240175, 0.13006862, 0.0595758,  0.0349345,  0,         0.02089832],
 [0.02027449, 0.07361198, 0.05271366, 0.13038054, 0.06331878, 0.03867748, 0.13038054, 0        ]]
'''

'''
# Kimura protein distance

A = [[0.,         0.07763265, 0.0550772,  0.14361677, 0.06626869, 0.03975659,  0.14361677, 0.0208878 ],
 [0.02217343, 0.,         0.05811257, 0.14248101, 0.06251925, 0.03811052,  0.14248101, 0.02217343],
 [0.02475157, 0.06900758, 0.,         0.02475157, 0.0666105,  0.03811052,  0.14475422, 0.02475157],
 [0.02313916, 0.07141248, 0.05777469, 0.,         0.06149998, 0.04272893,  0.13682754, 0.02313916],
 [0.02153033, 0.07417055, 0.05709939, 0.143238,   0.,         0.03811052,  0.143238,   0.02153033],
 [0.02120899, 0.07279023, 0.0560876,  0.14437488, 0.06388047, 0.,          0.14437488, 0.02120899],
 [0.02120899, 0.08846899, 0.05440438, 0.143238,   0.06217933, 0.03581225,  0.,         0.02120899],
 [0.02056675, 0.07763265, 0.05474071, 0.14361677, 0.06626869, 0.03975659,  0.14361677, 0.        ]]

'''

'''
#p-distance

A = [[0.,         0.08671241, 0.29101684, 0.17186525, 0.09357455, 0.09326263, 0.11166563, 0.00904554],
 [0.08671241, 0.,         0.29507174, 0.17747973, 0.1110418,  0.1110418, 0.10480349, 0.08296943],
 [0.29101684, 0.29507174, 0.,         0.26232065, 0.29351216, 0.29600749,0.30411728, 0.28789769],
 [0.17186525, 0.17747973, 0.26232065, 0.,         0.17654398, 0.17841547, 0.1933874,  0.16687461],
 [0.09357455, 0.1110418,  0.29351216, 0.17654398, 0.,         0.0246413, 0.12788522, 0.08889582],
 [0.09326263, 0.1110418,  0.29600749, 0.17841547, 0.0246413,  0., 0.13006862, 0.08951965],
 [0.11166563, 0.10480349, 0.30411728, 0.1933874,  0.12788522, 0.13006862, 0.,         0.10979414],
 [0.00904554, 0.08296943, 0.28789769, 0.16687461, 0.08889582, 0.08951965, 0.10979414, 0.        ]]

'''

A = np.array(A)

print A

iter=0
labels = ['Homo_sapiens','Dog','Zebrafish','Chicken','Rat','Mouse','Cattle','Chimpanzee']

while A is not None:
	print '\n############################Iteration : ', iter, ' #########################################\n'
	A,labels = iteration(A,labels,iter)
	iter+=1


print "\nTree obtained :\n"

for i in sorted(node_dictionary.keys()):
    print i, ':  ',
    node = node_dictionary[i]
    print node
