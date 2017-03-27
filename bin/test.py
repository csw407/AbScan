'''
Created on Dec 7, 2016

@author: s3cha
'''

class Borg:
    __shared_state = {}
    x = 0
    def __init__(self):
        self.__dict__ = self.__shared_state
    
        
    def setx(self,value):
        self.x = value
        pass
    
    def getx(self):
        return self.x
    

def compare(item1,item2):
    if item1[0] < item2[0]:
        return -1
    elif item1[0] > item2[0]:
        return 1
    else:
        return 0
    
# x = [(i,-i) for i in range(10,1,-1)]

# print sorted(x,cmp=compare)

import math
def compare_distance(item1,item2):
    '''compare the items for the distance list. first ordered by larged dif, second ordered by larged read depth 
    Note distance = [node index1, node index2, distance between them, max read depth]'''
    if math.fabs(item1[3]) > math.fabs(item2[3]):
        return -1
    elif math.fabs(item1[3]) < math.fabs(item2[3]):
        return 1
    elif item1[2] > item2[2]:
        return -1
    elif item1[2] < item2[2]:
        return 1
    else:
        return 0
    

# x = {1:[(0,1,2,3),(0,2,-2,4),(1,2,3,1),(0,3,1,1)]}
# for key in x:
#     print sorted(x[key],cmp=compare_distance)
    
from matplotlib import pyplot as plt
import math


f1 = open('/home/s3cha/data/Miin/dbsize/1DLC130826QE_UTI_#8_rep2_Homosapiens_unmatched_sorted_15.tsv','r')
f2 = open('/home/s3cha/data/Miin/dbsize/1DLC130826QE_UTI_#8_rep2_Homosapiens_unmatched_sorted_15_aaplus.tsv','r')

dist1 = [[],[]]
dist2 = [[],[]]

for line in f1:
    if line.startswith('#'):
        continue
    line = line.strip().split('\t')
    if line[10].startswith('XXX'):
        dist1[1].append(-math.log(float(line[13])))
    else:
        dist1[0].append(-math.log(float(line[13])))

for line in f2:
    if line.startswith('#'):
        continue
    line = line.strip().split('\t')
    if line[10].startswith('XXX'):
        dist2[1].append(-math.log(float(line[13])))
    else:
        dist2[0].append(-math.log(float(line[13])))

# plt.hist(dist1[0],alpha=0.3,normed = True,bins = 100,label='1target')
plt.hist(dist1[1],alpha=0.3,normed = True,bins = 50,label='1decoy')
# plt.hist(dist2[0],alpha=0.3,normed = True,bins = 100,label='2target')
plt.hist(dist2[1],alpha=0.3,normed = True,bins = 50,label='2decoy')
plt.legend()
plt.show()


    
