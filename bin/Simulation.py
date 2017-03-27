'''
Created on Dec 8, 2016

@author: s3cha
'''

import sys
sys.path.append(r'/home/s3cha/Dropbox/Workspace/IG/src/GraphConstruction/')
sys.path.append(r'/home/s3cha/Dropbox/Workspace/IG/pdbn')
sys.path.append(r'/home/s3cha/Dropbox/Workspace/IG/src/GraphConstruction/')
sys.path.append(r'/home/s3cha/Dropbox/Workspace/IG/src/Simulation/')
import os
import random
import cPickle as pickle
import copy
from matplotlib import pyplot as plt
from wgsim import WGSIM
from BuildMLgraph import *
from LevenshteinDistance import *
from Simulation1_wgsim_count_missing_true import *
from PDBn import PDBn
import ReadClass
import time
# from BuildMLgraph import AddRead
# from BuildMLgraph import AddRead3
# from BuildMLgraph import CleanMLGraph
# from BuildMLgraph import hamming_distance_different_length
# from BuildMLgraph import CheckGraphPath

class Simulation(object):
    '''
    classdocs
    '''
    wgsim = WGSIM()


    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def get_Reads(self,num_reads):
#         self.wgsim.
        return
    
    def sim(self):
        sequencing_error = 0.01
        coverage = .3
        read_length = 100
        reference = 'CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGAGGTACAACTGGAACGACGCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG'
    #     reference = 'CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGG' 
        lmer = 20
        rmer = 10
        
        num_reads = (len(reference)-read_length)*coverage
        print len(reference)
        print num_reads
        
        if 0:
            wgsim = WGSIM()
            wgsim.setErrorRate(sequencing_error)
            wgsim.setReadNumber(num_reads)
            wgsim.addRef(reference)
            READS = wgsim.MakeRead()
            pickle.dump(READS,open('test.p','wb'))
        else:
            READS = pickle.load(open('test.p','rb'))
            
        pdb = PDBn()
        pdb.__initialize__(lmer, rmer)
        read_class = ReadClass.SingleEnd()
        
    #     for read in READS:
    #         print read
            
        for read in READS:
            pdb.AddRead(read_class, read)
            
        pdb.NodeClustering()
        pdb.naiveCutThreshold(1)
    #     pdb.ShowGraph()
        print pdb.error_rate(reference, 10)
        print pdb.divergence()
        pdb.write()      
        pass
    

if __name__ == '__main__':
    plt.figure(1)
    simulation = Simulation()
#     simulation.sim()
    lmer = 40
    rmer = 20
    pdb = PDBn(rmer,lmer)
#     pdb.__initialize__(lmer,rmer)
    read_class = ReadClass.SingleEnd()
    
    f = open('/media/s3cha/MyBook/stefano/7_SAM/final_repertoire.clusters.fa','r')
    s = '/home/s3cha/1/testing/sample_splicegraph.txt'
    count = 0
    btime = time.time()
    for line in f:
        if line.startswith('>'):
            continue
        else:
            pdb.AddRead(read_class,line.strip())
        count += 1
        if count > 200:
            break
    plt.subplot(211)
    pdb.showNodeDepth()
    print 'Graph build time: ',time.time()-btime
    pdb.setLnodeSeq()################################
    btime = time.time()
    pdb.NodeClustering()#################################
    print 'Node clusting time:',time.time()-btime
    btime = time.time()
    plt.subplot(212)
    pdb.showNodeDepth()
    pdb.simplify()##############################
    pdb.tipClip()
    pdb.naiveCutThreshold(1)#########################
    pdb.simplify()##############################3
    print 'Node cut by threshold time:',time.time()-btime
    print pdb.divergence()
    pdb.writeSpliceGraph(s)
    print 222,pdb.test_sinkseq()
#     plt.show()
        
    
    
    


    
    
    
    
    