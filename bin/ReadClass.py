'''
Created on Nov 15, 2016

@author: s3cha
'''
import sys


class ReadClass(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        
    def AddRead(self,read_seq,Graph):
        #abstract method. use different method if the read set is not same. ex) single end read, paired end read, etc.
        raise NotImplementedError("Subclass should implement this")
            

class SingleEnd(ReadClass):
    def __init__(self):
        #initialize the readcut parameter which indicate how much read part will the graph contain at the end of the reads. 
        #ex) ABCDE. l = 2,r = 2,readcut = 0.5 >>> keep until l = 2*0.5 = 1. (AB,CD),(BC,DE),(CD,E)
        self.readcut = 0.5
        pass
    
    #for each read, chop up the read in to lmer and rmer piece, and connect the edges, and combine the nodes based on the rmer sequence identity.
    def AddRead(self,read_seq,Graph):
        lmer = Graph.get_Llength()
        rmer = Graph.get_Rlength()
        read_length = len(read_seq)
        #error check in case of wrong parameter.
        if read_length < lmer+rmer:
            print "Err: read length is shorter than sum of l,r. read length: %d, l,r: %d,%d" %(read_length,lmer,rmer)
            sys.exit()
        rmer_seq = read_seq[0:0+rmer]
        lmer_seq = read_seq[rmer:rmer+lmer]
        prev_lnode = Graph.getRNode(rmer_seq).getLNode(lmer_seq)
        #make node if not sepcified. connect edges between lnodes.
        for index in range(1,read_length-(lmer+rmer)+int(lmer*self.readcut)+1):
            rmer_seq = read_seq[index:index+rmer]
            lmer_seq = read_seq[index+rmer:index+rmer+lmer]
            current_lnode = Graph.getRNode(rmer_seq).getLNode(lmer_seq)
            current_lnode.ConnectPrev(prev_lnode)
            prev_lnode.ConnectNext(current_lnode)
            prev_lnode = current_lnode
        
        pass
    