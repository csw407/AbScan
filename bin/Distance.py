'''
Created on Dec 6, 2016

@author: s3cha
'''
import sys

class Distance(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    @staticmethod
    def HammingDistance(s1,s2):
        count = 0
        if len(s1) != len(s2):
    #         print s1,s2
    #         raise ValueError("Undefined for sequences of unequal length")
            return sys.maxint
        for i in range(len(s1)):
            if s1[i] != s2[i]:
                count += 1
        return count
    
    @staticmethod
    def HammingDistance_dif(s1,s2):
        long = s1 if len(s1) > len(s2) else s2
        short = s2 if len(s1) > len(s2) else s1
        distance = 0
        for i in range(len(short)):
            if short[i] != long[i]:
                distance += 1
        return distance
        
    
    @staticmethod
    def hamming_distance_different_length(s1, s2):
        count = 0
        if len(s1) != len(s2):
            if len(s1) > len(s2):
                dif = len(s1)-len(s2)
                for i in range(len(s2)):
                    if s1[dif+i] != s2[i]:
                        count += 1
            else:
                dif = len(s2)-len(s1)
                for i in range(len(s1)):
                    if s1[i] != s2[i+dif]:
                        count += 1
        else:
            for i in range(len(s1)):
                if s1[i] != s2[i]:
                    count += 1
        return count
    
    @staticmethod
    def LevenshteinDistance(s,t):
        if (s == t):
            return 0
        if len(s) == 0:
            return len(t)
        if len(t) == 0:
            return len(s)
        v0 = [0]*(len(t)+1)
        v1 = [0]*(len(t)+1)
        for i in range(len(v0)):
            v0[i] = i
        for i in range(len(s)):
            v1[0] = i+1
            for j in range(len(t)):
                if s[i] == t[j]:
                    cost = 0
                else:
                    cost = 1
                v1[j+1] = min(v1[j]+1,v0[j+1]+1,v0[j]+cost)
            for j in range(len(v0)):
                v0[j] = v1[j]
        return v1[len(t)]

if __name__ == '__main__':
    x = Distance()
    x.HammingDistance('abc','abe')
    print Distance.HammingDistance('abc','abe')