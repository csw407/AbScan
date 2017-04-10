'''
Created on Mar 23, 2017

@author: s3cha
'''
import sys
import os
system_folder = os.path.split(sys.argv[0])[0]
# sys.path.append(r'/home/s3cha/Dropbox/Workspace/IG/src/GraphConstruction/')
sys.path.append(r'%s/pdbn'%(system_folder))
# sys.path.append(r'/home/s3cha/Dropbox/Workspace/IG/src/GraphConstruction/')
# sys.path.append(r'/home/s3cha/Dropbox/Workspace/IG/src/Simulation/')

import random
import cPickle as pickle
import copy
from matplotlib import pyplot as plt
# from wgsim import WGSIM
# from BuildMLgraph import *
# from LevenshteinDistance import *
# from Simulation1_wgsim_count_missing_true import *
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
        self.threshold = 1
        self.naivecut = True
        pass
    
    def setParameter(self,rmer,lmer):
        self.sdbn.__initialize__(rmer,lmer)
    
    def setBamfile(self,bam_filename):
        self.bam_tracker.set_Bamfilename(bam_filename)
        pass
    
    def setOutputfile(self,output_filename):
        self.output_filename = output_filename
        pass
    
    def setThreshold(self,_threshold,_boolean):
        self.threshold = _threshold
        self.naivecut = _boolean
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
        self.setTime()
        if not self.output_filename:
            raise 'Output filename is not specified.'
        self.setBamfile(bam_filename)
        self.bam_tracker.set_Kmer(kmer) #kmer = parameter to compare the reads with reference
        self.bam_tracker.SdbConstruction(self.sdbn,self.read_class)
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
        if self.naivecut:
            self.sdbn.naiveCutThreshold(self.threshold)#########################
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
    
    def quick_test(self):
        '''
        write read to text: textout true
        read from text: textout false, read_from_text true
        read from bam : textout false, read_from_text false
        '''
        textout = False
        read_from_text = True
        input_bam = '/media/s3cha/MyBook/VU/bam/UNCID_1804970.2ff064bb-9e08-497e-a83a-7360cf80f62d.sorted_genome_alignments.bam'#UNCID_1580578.35e535d8-7265-4271-a503-5e1728cd5209.sorted_genome_alignments.bam'
        output_graph_filename = '/home/s3cha/1/testing/test_bam/sample_splicegraph_bam.txt'
        output_text = '/home/s3cha/1/testing/test_bam/filtered_read_from_bam_woMapped.txt'
        rmer = 10
        lmer = 20
        
        self.setParameter(rmer,lmer)
        self.setOutputfile(output_graph_filename)
        if textout:
            self.setMappedReadTrack(True)
            self.SaveBamToText(input_bam,output_text)
        elif read_from_text:
            self.graphConstruction_from_text(output_text)
        else:
            self.graphConstruction_from_bam(input_bam)            
#         self = SdbConstruction(rmer,lmer)
#         self.setMappedReadTrack(False)
#         self.SaveBamToText(input_bam,output_text)
#         self.graphConstruction_from_bam(input_bam)
    def ms2db_to_fasta(self,ms2db_folder,fasta_folder,param):
        file_list = os.listdir(ms2db_folder)
        default_command = system_folder.rstrip('/')+'/ACGT03102013_split_IgGraphConvert.py '
        for file in file_list:
            current = os.path.splitext(file)
            input_filename = ms2db_folder.rstrip('/')+'/'+file
            output_filename = fasta_folder.rstrip('/')+'/'+current[0]+'.fa'
            command_line = '%s %s %s %s'%(default_command, input_filename, output_filename, str(param))
            os.system(command_line) 
        

    def tmp_graphConstruction_from_igseq(self,text_filename):
        text_filename = '/media/s3cha/MyBook/stefano/7_SAM/final_repertoire.clusters.fa'
        output_graph_filename = '/home/s3cha/data/AbScan/Database/IG_seq/ms2db/IG_seq_polyclonalAB.ms2db'
        ms2db_folder = '/home/s3cha/data/AbScan/Database/IG_seq/ms2db'
        fasta_folder = '/home/s3cha/data/AbScan/Database/IG_seq/fasta'
        param = '2'
        rmer = 20
        lmer = 40
        self.setOutputfile(output_graph_filename)
        self.setParameter(rmer,lmer)
#         self.setThreshold(0,False)
        f = open(text_filename)
        div = 100000
        num = 1
        count = 0
        output_filename = os.path.splitext(self.output_filename)
        self.setTime()
        for line in f:
            if count > div*num:
                self.graphSimplification()
                num += 1
                self.setOutputfile(output_filename[0]+'_'+str(num)+output_filename[1]) 
                self.setTime()
                self.sdbn = PDBn(rmer,lmer)
            if line.startswith('#') or line.startswith('>'):
                continue
            count += 1
            self.sdbn.AddRead(self.read_class,line.strip())
        self.graphSimplification()
        self.ms2db_to_fasta(ms2db_folder,fasta_folder,param)
        pass
    
if __name__ == '__main__': 
    if len(sys.argv) < 2:
        x = SdbConstruction(10,20)
#         x.quick_test()
        x.tmp_graphConstruction_from_igseq(None)
    else:
        input_bam = sys.argv[1]
        output_graph = sys.argv[2]
        rmer = int(sys.argv[3])
        lmer = int(sys.argv[4])
        x = SdbConstruction(rmer,lmer)
        x.setOutputfile(output_graph)
        if os.path.splitext(input_bam)[1].lower() == '.bam':
            x.graphConstruction_from_bam(input_bam)
        else:
            x.graphConstruction_from_text(input_bam)
            
            
        

    
    

    