'''
Created on Mar 11, 2016

@author: s3cha
'''
import os
import sys
import cPickle as pickle
import pysam

# input_bam_filename = "/media/s3cha/MyBook/VU/bam/UNCID_1580578.35e535d8-7265-4271-a503-5e1728cd5209.sorted_genome_alignments.bam"
# s = open('/home/s3cha/data/SpliceDB/NEW_IG_Graph/KmerTest/temp_unmapped_read_list2.txt','w')

class ReadTrackAndFilter(object):
    
    def __init__(self):
        self.bam_filename = ''
        self.read_tracked_filename = ''
        self.p_filename = ''
        self.read_length = 76
        self.average_quality_threshold = 25
        self.quality_threshold = 10
        self.is_trackingmappedread = True
        self.vdj_ref_folder = ''
        self.sam_file = None
        self.Kmer = 19
        
    def set_Kmer(self,kmer):
        self.Kmer = kmer
    def set_Bamfilename(self,filename):
        self.bam_filename = filename
    
    def set_readtrackfilename(self,read_tracked_filename):
        self.read_tracked_filename = read_tracked_filename
    
    def set_pfilename(self,p_filename):
        self.p_filename = p_filename
        
    def set_mapreadtrack(self,check):
        self.is_trackingmappedread = check
    
    def set_qualitythreshold(self,threshold):
        self.quality_threshold = threshold
    
    def set_averagequalitythreshold(self,threshold):
        self.average_quality_threshold = threshold
        
    def get_ReverseComplement(self,sequence):
        map = {'G':'C','C':'G','T':'A','A':'T','N':'N'}
        complement_seq = []
        for index in range(len(sequence))[::-1]:
            complement_seq.append(map.get(sequence[index]))
        return ''.join(complement_seq)
        
    def get_RefKmer(self,ref_filename,Kmer):
        ref_seq = []
        kmer_seq = {}
        kmer_rc_seq = {}
        for file in ref_filename:
            f = open(file,'r')
            for line in f:
                if line.startswith('>'):
                    continue
                ref_seq.append(line.strip().upper())
                for index in range(len(line.strip())-Kmer):
                    cur_kmer = line[index:index+Kmer].upper()
                    kmer_seq[cur_kmer] = 0
                rev_line = self.get_ReverseComplement(line.strip().upper())
                for index in range(len(rev_line)-Kmer):
                    cur_kmer = rev_line[index:index+Kmer].upper()
                    kmer_rc_seq[cur_kmer] = 0
        return [kmer_seq,kmer_rc_seq,ref_seq]
    
    def LoadReference(self,Kmer):
        if self.vdj_ref_folder == '':
            reference_vfilename = '/home/s3cha/data/SpliceDB/create_ig_db/data/imgt_data/ig/func/human/human_IGHV.fa'
            reference_dfilename = '/home/s3cha/data/SpliceDB/create_ig_db/data/imgt_data/ig/func/human/human_IGHD.fa'
            reference_jfilename = '/home/s3cha/data/SpliceDB/create_ig_db/data/imgt_data/ig/func/human/human_IGHJ.fa'
            reference_cfilename = '/home/s3cha/data/SpliceDB/create_ig_db/data/imgt_data/ig/func/human/human_IGHC_oneline.fa'
        ref_filename = [reference_vfilename,reference_dfilename,reference_jfilename,reference_cfilename]
        kmer_seq,kmer_rc_seq,ref_seq = self.get_RefKmer(ref_filename, Kmer)
        return [kmer_seq,kmer_rc_seq,ref_seq]
    
    def LoadReference_DNA(self,Kmer):
        if self.read_tracked_filename == '':
            raise ValueError('Unmapped read file name is not specified')
        dna_chr14_filename = '/home/s3cha/s3cha/data/dna/human70/chr14_Homo_sapiens.GRCh37.70.dna.chromosome_formatted_chr14.trie'
        f = open(dna_chr14_filename,'r')
        dna_string = f.readline()
#         ref_seq = dna_string[106032614:107288051]#[105566277:106879844]#106032614,107288051
        ref_seq = dna_string[105566277:106879844]#106032614,107288051
        ref_rc_seq = self.get_ReverseComplement(ref_seq[::-1])
        kmer_seq = {}
        kmer_rc_seq = {}
        for index in range(len(ref_seq)-Kmer):
            kmer_seq[ref_seq[index:index+Kmer]] = 0
            kmer_rc_seq[ref_rc_seq[index:index+Kmer]] = 0
        
        return [kmer_seq,kmer_rc_seq,ref_seq]
        
    
    def TrackFilteredfile(self,Kmer,picklefilename): # find the read from the read_tracked_filename which has at least Kmer mapping to the reference
        
        kmer_seq,kmer_rc_seq,ref_seq = self.LoadReference(Kmer)
        
        print '# of reference sequence: ',len(ref_seq),', # of Kmer in the reference: ',len(kmer_seq)
        ACGT = ['A','C','G','T']
        
        kmer = Kmer
        f = open(self.read_tracked_filename,'r')
        filtered_read = []
        for line in f:
            line = line.strip()
            for index in range(len(line)-kmer):
                cur_kmer = line[index:index+kmer]
                if kmer_seq.has_key(cur_kmer):
                    filtered_read.append(line)
                    break
                elif kmer_rc_seq.has_key(cur_kmer):
                    filtered_read.append(self.get_ReverseComplement(line))
                    break     
        if picklefilename != '':
            pickle.dump(filtered_read, open(picklefilename,'wb'))
        return filtered_read
    
    def CheckKmerMap(self,kmer_seq,kmer_rc_seq,read,Kmer):
        check = False
        read_seq = read.query_sequence
        for i in range(len(read_seq)-Kmer):
            cur_kmer = read_seq[i:i+Kmer]
            if kmer_seq.has_key(cur_kmer):
                return True
            elif kmer_rc_seq.has_key(cur_kmer):
                return True
        return False
    
    '''def TrackBamfile(self): #tracking the bam file (mapped read to reference region, and unmapped read) and save it to read_tracked_filename
        if self.bam_filename == '':
            raise ValueError("BAM file name is not specified")
        if self.read_tracked_filename == '':
            raise ValueError("Read output file name is not specified")
        samfile = pysam.AlignmentFile(self.bam_filename,"rb")
        ###
        write_samfile = pysam.AlignmentFile("/home/s3cha/data/SpliceDB/NEW_IG_Graph/KmerTest/StefanoCompR/filtered_read_kmermap_imgtref.bam","wb",template=samfile)
        ###
        mapped_read_count = 0
        unmapped_read_count = 0
        ### trimming the 3' end of the read based on the quality value. input: quality scores as a list, output: index of position need to be trimmed
        def Trim_read_index(quality,threshold):
            index = self.read_length
            if sum(quality[:5]) < 5*threshold: #filter the 3' if the first 5 entries has lower quality values.
                return 0
            for i in quality[::-1]: #trim the 5' end if the quality value is lower than threshold.
                if i < threshold: 
                    index -= 1
                else:
                    break
            return index    
        s = open(self.read_tracked_filename,'w')
        read_length = self.read_length
        threshold = self.average_quality_threshold * read_length
        if self.is_trackingmappedread:
            count = 0
#             for read in samfile.fetch('chr14',106032614,107288051,until_eof=True):
            for read in samfile.fetch('chr14',105566277,106879844,until_eof=True):
                count += 1
                quality = read.query_qualities
                if sum(quality) < threshold:
                    continue
                
                ####
                write_samfile.write(read)
                
                trim_index = Trim_read_index(quality,self.quality_threshold)
                if trim_index < read_length * 2 / 3:
                    continue
                else:
                    read_seq = read.query_sequence[:trim_index]
                s.write(read_seq+'\n')
                mapped_read_count += 1
            print count
        count = 0
        for read in samfile.fetch(until_eof=True):
            count += 1
            if read.is_unmapped:
                quality = read.query_qualities
                
                ####
                write_samfile.write(read)
                
                
                if sum(quality) < threshold:
                    continue
                trim_index = Trim_read_index(quality,self.quality_threshold)
                if trim_index < read_length * 2 / 3:
                    continue
                else:
                    read_seq = read.query_sequence[:trim_index]
                s.write(read_seq+'\n')
                unmapped_read_count += 1
        print count
        write_samfile.close()
        s.close()
        samfile.close()
        
        return [mapped_read_count,unmapped_read_count]'''
    
    '''def TrackBamfile_withFilter(self,Kmer,pfilename): #tracking the bam file (mapped read to reference region, and unmapped read) and save it to read_tracked_filename
        if self.bam_filename == '':
            raise ValueError("BAM file name is not specified")
        if self.read_tracked_filename == '':
            raise ValueError("Read output file name is not specified")
        samfile = pysam.AlignmentFile(self.bam_filename,"rb")
        ###
#         kmer_seq, kmer_rc_seq, ref_seq = self.LoadReference_DNA(Kmer)
        kmer_seq, kmer_rc_seq, ref_seq = self.LoadReference(Kmer)
        write_samfile = pysam.AlignmentFile("/home/s3cha/data/SpliceDB/NEW_IG_Graph/KmerTest/StefanoCompR/filtered_read_kmermap_imgtrefer.bam","wb",template=samfile)
        ###
        mapped_read_count = 0
        unmapped_read_count = 0
        ### trimming the 3' end of the read based on the quality value. input: quality scores as a list, output: index of position need to be trimmed
        def Trim_read_index(quality,threshold):
            index = self.read_length
            if sum(quality[:5]) < 5*threshold: #filter the 3' if the first 5 entries has lower quality values.
                return 0
            for i in quality[::-1]: #trim the 5' end if the quality value is lower than threshold.
                if i < threshold: 
                    index -= 1
                else:
                    break
            return index    
        s = open(self.read_tracked_filename,'w')
        read_length = self.read_length
        threshold = self.average_quality_threshold * read_length
        filtered_read = []
        if self.is_trackingmappedread:
            count = 0
#             for read in samfile.fetch('chr14',106032614,107288051,until_eof=True):
            for read in samfile.fetch('chr14',105566277,106879844,until_eof=True):
                quality = read.query_qualities
                if sum(quality) < threshold:
                    continue
                
                ####
                if self.CheckKmerMap(kmer_seq, kmer_rc_seq, read, Kmer):
                    count += 1
                    write_samfile.write(read)
                else:
                    continue
                
                trim_index = Trim_read_index(quality,self.quality_threshold)
                if trim_index < read_length * 2 / 3:
                    continue
                else:
                    read_seq = read.query_sequence[:trim_index]
                s.write(read_seq+'\n')
                filtered_read.append(read_seq)
                mapped_read_count += 1
            print count
        count = 0
        for read in samfile.fetch(until_eof=True):
            
            if read.is_unmapped:
                quality = read.query_qualities
                
                ####
                if self.CheckKmerMap(kmer_seq, kmer_rc_seq, read, Kmer):
                    count += 1
                    write_samfile.write(read)
                else:
                    continue
                
                
                if sum(quality) < threshold:
                    continue
                trim_index = Trim_read_index(quality,self.quality_threshold)
                if trim_index < read_length * 2 / 3:
                    continue
                else:
                    read_seq = read.query_sequence[:trim_index]
                s.write(read_seq+'\n')
                filtered_read.append(self.get_ReverseComplement(read_seq))
                unmapped_read_count += 1
        print count
        write_samfile.close()
        s.close()
        samfile.close()
        if pfilename != '':
            print len(filtered_read)
            pickle.dump(filtered_read,open(pfilename,'wb'))
        return [mapped_read_count,unmapped_read_count]'''
                
    
    '''def Test_input(self,Kmer,pfilename): #tracking the bam file (mapped read to reference region, and unmapped read) and save it to read_tracked_filename
        if self.bam_filename == '':
            raise ValueError("BAM file name is not specified")
        if self.read_tracked_filename == '':
            raise ValueError("Read output file name is not specified")
        samfile = pysam.AlignmentFile(self.bam_filename,"rb")
        ###
#         kmer_seq, kmer_rc_seq, ref_seq = self.LoadReference_DNA(Kmer)
        kmer_seq, kmer_rc_seq, ref_seq = self.LoadReference(Kmer)
        write_samfile = pysam.AlignmentFile("/home/s3cha/data/SpliceDB/NEW_IG_Graph/KmerTest/StefanoCompR/filtered_read_kmermap_test.bam","wb",template=samfile)
        ###
        mapped_read_count = 0
        unmapped_read_count = 0
        ### trimming the 3' end of the read based on the quality value. input: quality scores as a list, output: index of position need to be trimmed
        def Trim_read_index(quality,threshold):
            index = self.read_length
#             if sum(quality[:5]) < 5*threshold: #filter the 3' if the first 5 entries has lower quality values.
#                 return 0
            for i in quality[::-1]: #trim the 5' end if the quality value is lower than threshold.
                if i < threshold: 
                    index -= 1
                else:
                    break
            return index    
        s = open(self.read_tracked_filename,'w')
        read_length = self.read_length
        threshold = self.average_quality_threshold
        filtered_read = []
        if self.is_trackingmappedread:
            count = 0
#             for read in samfile.fetch('chr14',106032614,107288051,until_eof=True):
            for read in samfile.fetch('chr14',105566277,106879844,until_eof=True):
                if self.CheckKmerMap(kmer_seq, kmer_rc_seq, read, Kmer):
                    count += 1
#                     write_samfile.write(read)
                else:
                    continue
                quality = read.query_qualities
                
                if float(sum(quality))/read_length < threshold:
                    continue
                write_samfile.write(read)
                ####
                
                
                trim_index = Trim_read_index(quality,self.quality_threshold)
                if trim_index < read_length * 2 / 3:
                    continue
                else:
#                     if trim_index < read_length:
#                         print read.query_sequence
#                         sys.exit()
                    read.query_sequence = read.query_sequence[:trim_index]+'A'*(read_length-trim_index)
                    read_seq = read.query_sequence
#                     if trim_index < read_length:
#                         
#                         print read_seq
#                         print len(read_seq)
#                         sys.exit()
#                 write_samfile.write(read)
                s.write(read_seq+'\n')
                filtered_read.append(read_seq)
                mapped_read_count += 1
            print count
        count = 0
        for read in samfile.fetch(until_eof=True):
            
            if read.is_unmapped:
                quality = read.query_qualities
                
                ####
                if self.CheckKmerMap(kmer_seq, kmer_rc_seq, read, Kmer):
                    count += 1
#                     write_samfile.write(read)
                else:
                    continue
                
                write_samfile.write(read)
                if float(sum(quality))/read_length < threshold:
                    continue
                
                trim_index = Trim_read_index(quality,self.quality_threshold)
                if trim_index < read_length * 2 / 3:
                    continue
                else:
                    read.query_sequence = read.query_sequence[:trim_index] + 'A'*(read_length-trim_index)
                    read_seq = read.query_sequence
#                 write_samfile.write(read)
                s.write(read_seq+'\n')
                filtered_read.append(self.get_ReverseComplement(read_seq))
                unmapped_read_count += 1
        print count
        write_samfile.close()
        s.close()
        samfile.close()
        if pfilename != '':
            print len(filtered_read)
            pickle.dump(filtered_read,open(pfilename,'wb'))
        return [mapped_read_count,unmapped_read_count]
    '''
    
    
    def Trim_read_index(self,quality,threshold):
        index = self.read_length
        for i in quality[::-1]: #trim the 5' end if the quality value is lower than threshold.
            if i < threshold: 
                index -= 1
            else:
                break
        return index
    
    class ReadProcessor(object):
        def __init__(self):
            pass
        def processRead(self,read):
            raise NotImplementedError("Subclass should implement this")
            pass
        
    class SaveToText(ReadProcessor):
        def __init__(self,text_filename):
            self.s = open(text_filename,'w')
            pass
        
        def processRead(self,read):
            self.s.write(read+'\n')
            pass
    
    class SdbnAddRead(ReadProcessor):
        def __init__(self,sdbn,readclass):
            self.sdbn = sdbn
            self.readclass = readclass
            pass
        
        def processRead(self,read):
            self.sdbn.AddRead(self.readclass,read)
            pass
        
    def SdbConstruction(self,sdbn,readclass):
        processor = self.SdbnAddRead(sdbn,readclass)
        self.BamfileProcessor(processor)
        pass
    
    def SaveReadToText(self,text_filename):
        processor = self.SaveToText(text_filename)
        self.BamfileProcessor(processor)
        pass
    
    def BamfileProcessor(self,read_processor):
        if self.bam_filename == '':
            raise ValueError("BAM file name is not specified")
        samfile = pysam.AlignmentFile(self.bam_filename,"rb")
        kmer_seq, kmer_rc_seq, ref_seq = self.LoadReference(self.Kmer)
        mapped_read_count = 0
        unmapped_read_count = 0
        read_length = self.read_length
        threshold = self.average_quality_threshold
        
        if self.is_trackingmappedread:
            count = 0
            for read in samfile.fetch('chr14',105566277,106879844,until_eof=True):
                if self.CheckKmerMap(kmer_seq, kmer_rc_seq, read, self.Kmer):
                    count += 1
                else:
                    count += 1
                    continue
                quality = read.query_qualities
                if float(sum(quality))/read_length < threshold:
                    continue
                trim_index = self.Trim_read_index(quality,self.quality_threshold)
                if trim_index < read_length * 2 / 3:
                    continue
                else:
                    read.query_sequence = read.query_sequence[:trim_index]
                    read_seq = read.query_sequence
                    
#                 filtered_read.append(read_seq)
#                 sdbn.AddRead(readclass,read_seq)
                read_processor.processRead(read_seq)
                mapped_read_count += 1
            print 'Number of mapped reads: %d, Percent of reads filtered: %f '%(mapped_read_count,float(mapped_read_count)/count*100)
        count = 0
        for read in samfile.fetch(until_eof=True):
            if read.is_unmapped:
                quality = read.query_qualities
                if self.CheckKmerMap(kmer_seq, kmer_rc_seq, read, self.Kmer):
                    count += 1
                else:
                    count += 1
                    continue
                if float(sum(quality))/read_length < threshold:
                    continue
                trim_index = self.Trim_read_index(quality,self.quality_threshold)
                if trim_index < read_length * 2 / 3:
                    continue
                else:
                    read.query_sequence = read.query_sequence[:trim_index]
                    read_seq = read.query_sequence
#                 filtered_read.append(self.get_ReverseComplement(read_seq))
#                 sdbn.AddRead(readclass,self.get_ReverseComplement(read_seq))
                read_processor.processRead(self.get_ReverseComplement(read_seq))
                unmapped_read_count += 1
        print 'Number of unmapped reads: %d, Percent of reads filtered: %f '%(unmapped_read_count,float(unmapped_read_count)/count*100)
        samfile.close()
        return [mapped_read_count,unmapped_read_count]
        
#105566277 and POS < 106879844  
#106032614,107288051 
# input_bam = '/home/s3cha/data/SpliceDB/NEW_IG_Graph/KmerTest/StefanoCompR/filtered_read_kmermap_test.bam'
# x.TrackBamfile()
# x.TrackBamfile_withFilter(19,'/home/s3cha/data/SpliceDB/NEW_IG_Graph/KmerTest/uncid/filtered_reads.p')
# x.TrackFilteredfile(19, '/home/s3cha/data/SpliceDB/NEW_IG_Graph/KmerTest/uncid/filtered_reads.p')

if __name__ == '__main__':
    print 123
    input_bam = '/media/s3cha/MyBook/VU/bam/UNCID_1580578.35e535d8-7265-4271-a503-5e1728cd5209.sorted_genome_alignments.bam'
    x = ReadTrackAndFilter()
    x.set_Bamfilename(input_bam)
    
    
    
'''
input_bam = '/media/s3cha/MyBook/VU/bam/UNCID_1580578.35e535d8-7265-4271-a503-5e1728cd5209.sorted_genome_alignments.bam'
x = ReadTrackAndFilter()
x.set_Bamfilename(input_bam)
x.set_readtrackfilename('/home/s3cha/data/SpliceDB/NEW_IG_Graph/KmerTest/temp_unmapped_read_list3.txt')
x.Test_input(19,'/home/s3cha/data/SpliceDB/NEW_IG_Graph/KmerTest/uncid/filtered_reads.p')
'''



