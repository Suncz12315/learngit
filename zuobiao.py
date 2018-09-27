import random
import itertools
import numpy as np
import sys


W = 20000 #带宽
H = 100  #无人机高度
v_max = 26 #无人机最大速度
err = 0.00001
beta = 10**8
alpha = 2
En = [1,1.2,0.7,1.4]#每个节点的最大能量
B = [1.2*10**6,0.5*10**6,2.1*10**6,1.1*10**6]#每个节点需要传输的信息量

node_num = 4 #待采集节点数量
coords = list(itertools.product(range(1, 5000), range(1, 5000)))
L = random.sample(coords,node_num)  #生成随机坐标
L.insert(0,[0.0,0.0])
coor = {}   # 存坐标
#print(random.sample(random_n_list, n))
for i in range(5):
    coor[i] = L[i]
def dis(L1,L2):
    distance = np.sqrt((L1[0]-L2[0])**2+(L1[1]-L2[1])**2)
    return distance
#print(coor)
a = [[0]*(node_num+1)]*(node_num+1)
dist = np.array(a,dtype = np.float)
t_min = np.array(a,dtype = np.float)
#print(a)
for i in range(node_num + 1):
    for j in range(node_num + 1):
        dist[i][j] = dis(coor[i],coor[j])
print(dist)

#计算满足条件的最小速度
def v_m(x,d,E):
    v = 2*(d-x)**3/(3*beta*E)
    return v

#计算在某一速度下的信息量
def fomu_B(x,v,d,E):
    gamma_0 = fomugamma_0(x,v,d,E)
    max_B = W/(2*v)*((d-x)*np.log2(beta*gamma_0/((d-x)**2+H*H)**(alpha/2)+alpha*(d-x)/np.log(2)-alpha*H/np.log(2)*np.arctan((d-x)/H)))
    return max_B

def fomugamma_0(x,v,d,E):
    gamma_0=v*E/(d-x)+(d-x)**2/(3*beta)+H**2/beta
    return gamma_0

#寻找最优的v
def find_v(v,v_max,B,x,d):
    v_l = v
    v_r = v_max
    while(v_l<v_r):
        mid_v = v_l+(v_r-v_l)/2
        #print(mid_v)
        B_op = fomu_B(x,mid_v,d,E)
        #print(B_op)
        if(abs(B_op-B)<err):
            return mid_v
        elif(B_op<B):
            v_r = mid_v
        else:
            v_l = mid_v


for i in range(node_num+1):
    for j in range(1,node_num+1):
      t_min_m = sys.maxsize
      if(i != j):
        for x in np.arange(0,dist[i][j],1):
            if(i != j):
                d = dist[i][j]
                #print(d)
                E = En[j-1]
                v_min = v_m(x,d,E)
            #print(dist[i][j])
                B_max = fomu_B(x,v_min,d,E)
                if(B_max>B[j-1]):
                    v_opt  = find_v(v_min,v_max,B[j-1],x,d)
                    #print(v_opt)
                    if(v_opt):
                        t_op = (d-x)/v_opt+x/v_max
                        #print(t_op)
                        t_min_m = min(t_min_m,t_op)
        t_min[i][j] = t_min_m
for i in range(1,node_num + 1):
    t_min[i][0] = dist[i][0]/v_max
print(t_min)

