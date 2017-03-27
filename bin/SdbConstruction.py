'''
Created on Mar 23, 2017

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
import ReadTrackAndFilter


class SdbConstruction(object):
    def __init__(self,rmer,lmer):
        self.sdbn = PDBn(rmer,lmer)
        self.read_class = ReadClass.SingleEnd()
        self.bam_tracker = ReadTrackAndFilter.ReadTrackAndFilter()
        self.output_filename = None
        self.btime = self.setTime()
        pass
    
    def setBamfile(self,bam_filename):
        self.bam_tracker.set_Bamfilename(bam_filename)
        pass
    
    def setOutputfile(self,output_filename):
        self.output_filename = output_filename
        pass
    
    def graphConstruction_from_text(self,text_filename):
        self.setTime()
        f = open(text_filename)
        for line in f:
            if line.startswith('#') or line.startswith('>'):
                continue
            self.sdbn.AddRead(self.read_class,line.strip())
        self.graphSimplification()
        pass
    
    def graphConstruction_from_bam(self,bam_filename,kmer = 19):
        if not self.output_filename:
            raise 'Output filename is not specified.'
        self.setBamfile(bam_filename)
        self.bam_tracker.set_Kmer(kmer) #kmer = parameter to compare the reads with reference
        mapped_read_count, unmapped_read_count = self.bam_tracker.SdbConstruction(self.sdbn,self.read_class)
        self.graphSimplification()
        pass
    
    def setTime(self):
        self.btime = time.time()
        pass
        
    def getTime(self):
        begin = self.btime
        self.setTime()
        return time.time()-begin
    
    def graphSimplification(self):
        print 'Graph build time: %f (sec)'%(self.getTime())
        self.sdbn.setLnodeSeq()
        self.sdbn.NodeClustering()
        print 'Node clustering time: %f (sec)'%(self.getTime())
        self.sdbn.simplify()
        self.sdbn.tipClip()
        self.sdbn.naiveCutThreshold(1)#########################
        self.sdbn.simplify()
        print 'Graph simplification time: %f (sec)'%(self.getTime())
        self.sdbn.writeSpliceGraph(self.output_filename)
        print 'End.'
        
    def SaveBamToText(self,bam_filename,text_filename):
        self.setBamfile(bam_filename)
        self.bam_tracker.SaveReadToText(text_filename)
        sys.exit()
        pass
    
    def setMappedReadTrack(self,check):
        self.bam_tracker.set_mapreadtrack(check)
        pass
    
    def quick_test(self,x):
        input_bam = '/media/s3cha/MyBook/VU/bam/UNCID_1580578.35e535d8-7265-4271-a503-5e1728cd5209.sorted_genome_alignments.bam'
        output_graph_filename = '/home/s3cha/1/testing/test_bam/sample_splicegraph_bam.txt'
        output_text = '/home/s3cha/1/testing/test_bam/filtered_read_from_bam_woMapped.txt'
        rmer = 10
        lmer = 20
        
        
#         x = SdbConstruction(rmer,lmer)
        x.setOutputfile(output_graph_filename)
    #     x.setMappedReadTrack(False)
        x.graphConstruction_from_text(output_text)
    #     x.SaveBamToText(input_bam,output_text)
    #     x.graphConstruction_from_bam(input_bam)
    
    
if __name__ == '__main__': 
    input_bam = '/media/s3cha/MyBook/VU/bam/UNCID_1580578.35e535d8-7265-4271-a503-5e1728cd5209.sorted_genome_alignments.bam'
    output_graph_filename = '/home/s3cha/1/testing/test_bam/sample_splicegraph_bam.txt'
    output_text = '/home/s3cha/1/testing/test_bam/filtered_read_from_bam_woMapped.txt'
    rmer = 10
    lmer = 20
    
    
    x = SdbConstruction(rmer,lmer)
    x.setOutputfile(output_graph_filename)
#     x.setMappedReadTrack(False)
    x.graphConstruction_from_text(output_text)
#     x.SaveBamToText(input_bam,output_text)
#     x.graphConstruction_from_bam(input_bam)
    