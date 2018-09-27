import random
import itertools
import numpy as np
#from dis import dis
n = 4
coords = list(itertools.product(range(1, 5000), range(1, 5000)))
L = random.sample(coords,n)  #生成随机坐标
L.insert(0,[0,0])
coor = {}   # 存坐标
E = [1,1.2,0.7,2]
B = [1.2*10**6,0.5*10**6,2.1*10**6,1.1*10**6]
#print(random.sample(random_n_list, n))
for i in range(5):
    coor[i] = L[i]
def dis(L1,L2):
    distance = np.sqrt((L1[0]-L2[0])**2+(L1[1]-L2[1])**2)
    return distance
#print(coor)
a = [[0]*(n+1)]*(n+1)
dist = np.array(a)
#print(a)
for i in range(5):
    for j in range(5):
        dist[i][j] = int(dis(coor[i],coor[j]))
print(dist)

