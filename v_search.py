import numpy as np
import sys
W = 20000
beta = 10**8
H = 100
v_max = 26
err = 0.00001
d = 3000
E = 1
B = 2.5*10**6
alpha = 2
#v = 2*d**3/(3*beta*E)

def v_m(x):
    v = 2*(d-x)**3/(3*beta*E)
    return v
def fomu_B(x,v):
    max_B = W/(2*v)*((d-x)*np.log2(beta*fomugamma_0(x,v)/((d-x)**2+H*H)**(alpha/2)+alpha*(d-x)/np.log(2)-alpha*H/np.log(2)*np.arctan((d-x)/H)))
    return max_B
#v = d**3/(12*beta*E)
#print(v)
def fomugamma_0(x,v):
    gamma_0=v*E/(d-x)+(d-x)**2/(3*beta)+H**2/beta
    return gamma_0
                     
def find_v(v,v_max,B):
    v_l = v
    v_r = v_max
    while(v_l<v_r):
        mid_v = v_l+(v_r-v_l)/2
        B_op = fomu_B(x,mid_v)
        if(abs(B_op-B)<err):
            return mid_v
        elif(B_op<B):
            v_r = mid_v
        else:
            v_l = mid_v
t_min = sys.maxsize
dict1 = {}                  
for x in range(0,d,10):
    v_min = v_m(x)
    #print(v_min)
    B_max = fomu_B(x,v_min)
    #print(B_max)
    if(B_max>B):
        v_opt = find_v(v_min,v_max,B)
        #print(v_opt)
        if(v_opt):
            t_op = (d-x)/v_opt+x/v_max
            #print(t_op)
            dict1[t_op] = x
            t_min = min(t_min,t_op)

print(t_min)
print(dict1[t_min])
def fomu_B_H(t,E):
    B_H = 1/2*t*W*np.log2(1+beta*E/(t*H*H))
    return B_H
def find_t(t,B):
    t_l = 0
    t_r = t
    while(t_l<t_r):
        mid_t = t_l+(t_r-t_l)/2
        B_H = fomu_B_H(mid_t,E)
        if(abs(B_H-B)<err):
            return mid_t
        elif(B_H<B):
            t_l = mid_t
        else:
            t_r = mid_t
B_H = fomu_B_H(t_min,E)
print(B_H)
if(B_H>B):
    t_hover = find_t(t_min,B)+d/v_max
    print(t_hover)
else:
    print("not optimal")

     


