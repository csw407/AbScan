'''
Created on Nov 15, 2016

@author: s3cha
'''
import math
from string import digits
import sys
import time

import Distance
import RNode
import ReadClass
import StackNQueue
# import SdbConstruction


# import Simulation
# import time
# from matplotlib import pyplot as plt
# import cPickle as pickle
# from ReadClass import SingleEnd
class PDBn(object):
    '''
    paired deBruijn graph class.
    Uses RNode and LNode class for a node information.
    Uses ReadClass to adjust different type of input reads format and structure. (For now it only allow to use the single ended reads)
    '''
#     __shared_state = {}
    is_initialized = False
    length_lmer = None
    length_rmer = None
    Rmer_node_map = None
    error_tolerance = None
    
    def __init__(self,length_rmer,length_lmer):
#         self.__dict__ = self.__shared_state
        self.length_lmer = length_lmer
        self.length_rmer = length_rmer
        self.Rmer_node_map = {}
        self.error_tolerance = 1
        
        self.t1 = float(0)
        self.t2 = float(0)
    
    def __initialize__(self,length_rmer,length_lmer):
        #Rmer node map connects the Rmer sequence to the corresponding object class
        self.length_lmer = length_lmer
        self.length_rmer = length_rmer
        self.Rmer_node_map = {}
        self.error_tolerance = 1
        
    def check(self):
        print self.length_lmer, self.length_rmer
        pass
    
    def get_Llength(self):
        return self.length_lmer
    
    def get_Rlength(self):
        return self.length_rmer
    
    def get_error_tolerance(self):
        return self.error_tolerance
    
    #add reads using the proper way defined by the read class. (single ended read, paired ended reads, or etc)
    def AddRead(self,read_class,read_seq):
        read_class.AddRead(read_seq,self)
        pass
        
    #return corresponding rnode based on the rmer sequence.
    def getRNode(self,rmer_seq):
        rnode = self.Rmer_node_map.get(rmer_seq)
        if rnode == None:
            rnode = RNode.RNode(rmer_seq,self)
            self.Rmer_node_map[rmer_seq] = rnode
        return rnode
    
    #write the graph to test the current graph structure. (For decoding)
    def ShowGraph(self):
        for rnode in self.Rmer_node_map.values():
            print 'Rnode: ',rnode.Rmer_seq
            for lnode in rnode.LNodes.values():
                print 'Lnode(prev,next): ',lnode.Lmer_seq,[(node.Lmer_seq,lnode.prev_nodes.get(node)) for node in lnode.prev_nodes],[(node.Lmer_seq,lnode.next_nodes.get(node)) for node in lnode.next_nodes]
            print
    
    #find source nodes
    def FindSourceNode(self):
        source = []
        for node in self.Rmer_node_map.values():
            is_source = True
            for lnode in node.getLNodeList():
                if lnode.GetPrev():
                    is_source = False
                    break
            if is_source:
                source.append(node)
        return source
    
    #cluster the l-nodes if their difference is less then error-tolerance
    def NodeClustering(self):
        st = StackNQueue.Queue()
        self.resetVisited()
        [st.enqueue(node) for node in self.FindSourceNode()]
        for node in self.Rmer_node_map.values():
            if not node.get_visit():
                st.enqueue(node)
            while st.peek():
                node = st.dequeue()
                if node == None:
                    break
                if node.get_visit():
                    continue
                node.set_visit()
                node.MergeLNodes()
                next_list = []                    
                for lnode in node.getLNodeList():
                    next_list += lnode.GetNext().keys()
                next_list = set(next_list)
                for next_node in next_list:
                    temp_rnode = self.Rmer_node_map.get(next_node.GetRmer())
#                     print temp_rnode
                    if temp_rnode.get_visit():
                        continue
                    st.enqueue(temp_rnode)
        pass
    
    #set all l-nodes to include only the last sequence of the node. (r-node + l-node sequence if the node is source
    def setLnodeSeq(self):
        for rnode in self.Rmer_node_map.values():
            for lnode in rnode.getLNodeList():
                if lnode.GetPrev():
                    lnode.setSeq(lnode.GetRmer()[-1])
                else:
                    lnode.setSeq(lnode.GetRmer())
                if lnode.GetNext():
                    continue
                else:
                    lnode.addSeqTail(lnode.GetLmer().translate(None,digits).strip())
        
        ###debuging. For some reasons, iterating graph failed to traverse all nodes
#         source = self.findSourceLnode()
#         st = StackNQueue.Queue()
#         self.resetVisited()
#         [st.enqueue(node) for node in source]
#         while st.peek():
#             node = st.dequeue()
#             if node.get_visit():
#                 continue
#             node.set_visit()
#             if 'NA' in node.getSeq():
#                 print 'err',node,node.GetRmer(),node.GetLmer()
#                 print node.GetLmer(),node.getSeq()
#                 for node1 in self.Rmer_node_map.get(node.GetRmer()).getLNodeList():
#                     print node1.GetLmer(),node1.getSeq()
#                 print node.GetPrev(),node.GetNext()
#                 sys.exit()
#             for next in node.GetNext():
#                 st.enqueue(next)
        pass
    
    def test_sequence(self):
        source = self.findSourceLnode()
        st = StackNQueue.Queue()
        self.resetVisited()
        [st.enqueue(node) for node in source]
        while st.peek():
            node = st.dequeue()
            if node.get_visit():
                continue
            node.set_visit()
            if node.getSeq() == 'NA':
                print 'err',node,node.GetRmer(),node.GetLmer()
                print node.GetLmer(),node.getSeq()
                for node1 in self.Rmer_node_map.get(node.GetRmer()).getLNodeList():
                    print node1.GetLmer(),node1.getSeq()
                print node.GetPrev(),node.GetNext()
                return False
            for next in node.GetNext():
                st.enqueue(next)
        return True
    
    def findSourceLnode(self):
        source = []
        for rnode in self.Rmer_node_map.values():
            for lnode in rnode.getLNodeList():
                if not lnode.GetPrev():
                    source.append(lnode)
        return source
    
    def findSinkLnode(self):
        sink = []
        for rnode in self.Rmer_node_map.values():
            for lnode in rnode.getLNodeList():
                if not lnode.GetNext():
                    sink.append(lnode)
        return sink
    
    def test_sinkseq(self):
        sink = self.findSinkLnode()
        for node in sink:
            if len(node.getSeq()) < len(node.GetLmer()):
                return False
        return True
    
    def temp_test(self):
        #self.setLnodeSeq()
        source = self.findSourceLnode()
        st = StackNQueue.Queue()
        [st.enqueue(node) for node in source]
        def test_connection(Lnode):
            for prev in Lnode.GetPrev():
                if Lnode not in prev.GetNext():
                    return False
            for next in Lnode.GetNext():
                if Lnode not in next.GetPrev():
                    return False
            return True
        while st.peek():
            node = st.dequeue()
            if len(node.getSeq()) < self.length_rmer:
                print 'Test-Source node length is shorter than r:',node.getSeq(),[x.GerRmer() for x in node.GetPrev()],node.GetRmer(),node.GetLmer()
                return False
        for Rnode in self.Rmer_node_map.values():
            for Lnode in Rnode.getLNodeList():
                if not test_connection(Lnode):
                    print 'Test-Connection is not matched:',Rnode.GetRmer(),Lnode.GetLmer()
                    return False
        del st
        return True
    
    def test_connection(self):
        source = self.findSourceLnode()
        st = StackNQueue.Queue()
        [st.enqueue(node) for node in source]
        def test_connection(Lnode):
            for prev in Lnode.GetPrev():
                if Lnode not in prev.GetNext():
                    return False
            for next in Lnode.GetNext():
                if Lnode not in next.GetPrev():
                    return False
            return True
        visited_node = {}
        while st.peek():
            node = st.dequeue()
            if node in visited_node:
                continue
            else:
                visited_node[node] = True
            if not test_connection(node):
                print 'Test_connection is not matched:', node.GetRmer(),node.GetLmer()
                return False
            for next_node in node.GetNext():
                st.enqueue(next_node)
        del visited_node, st
        return True
    
    #simplify the graph (combine all nodes if they include only a single path)
    '''def simplify_debug(self):
        self.setLnodeSeq()
#         temp = self.Rmer_node_map['GGCTCCAGGGAAGGG'].LNodes['CCTGGAGTGGGTTGGCCGCATTAGAAACAA']
#         print temp.getSeq(),temp.GetRmer(),temp.GetLmer()
#         sys.exit()
        
        source = self.findSourceLnode()
        st = StackNQueue.Queue()
        [st.enqueue(node) for node in source]
        def is_one_to_one(current):
            if len(current.GetPrev()) <= 1 and len(current.GetNext()) <= 1:
                return True
            else:
                return False
        def find_tail(current):
            pndepth = current.GetPrev().values()+current.GetNext().values()
            max_depth = max(pndepth) if pndepth else 0
            if is_one_to_one(current) and len(current.GetNext()) == 1:
                tail = current.GetNext().keys()[0]
                current.addSeqTail(''.join(tail.getSeq()))
                if is_one_to_one(tail):
                    while True:
                        if len(tail.GetNext()) == 0:
                            break
                        elif is_one_to_one(tail.GetNext().keys()[0]):
                            max_depth = max(max_depth,tail.GetNext().values()[0])
                            temp_tail = tail
                            tail = tail.GetNext().keys()[0]
                            temp_tail.remove()
                            current.addSeqTail(''.join(tail.getSeq()))
                        else:
                            break
                else:
                    return None
            else:
                return None
            return [tail,max_depth]
        def mergePath(head,tail,max_depth):
            if len(head.GetPrev()) != 0:
                prev_head = head.GetPrev().keys()[0]
                head.GetPrev()[prev_head] = max_depth
                prev_head.GetNext()[head] = max_depth
            if len(tail.GetNext()) != 0:
                next_tail = tail.GetNext().keys()[0]
                next_tail.GetPrev()[head] = max_depth
#                 del next_tail.GetPrev()[tail]
#                 del head.GetNext()[tail]
                head.GetNext()[next_tail] = max_depth
            pass
        while st.peek():
            node = st.dequeue()
            if node.get_visit():
                continue
#             print node.GetRmer(),node.GetLmer(),node.getSeq()
            node.set_visit()
            tails = find_tail(node)
            if tails != None:
                mergePath(node,tails[0],tails[1])
                tails[0].remove()
                node = tails[0]
            for nextnode in node.GetNext().keys():
                st.enqueue(nextnode)   
        pass'''
    
    def simplify_node(self,node):
        def is_one_to_one(current):
            if len(current.GetPrev()) <= 1 and len(current.GetNext()) <= 1:
                return True
            else:
                return False
        
        def find_trail(path,seq,maxdepth):
            c_node = path[-1]
            if not is_one_to_one(c_node):
                path.pop()
                seq.pop()
                return (path,''.join(seq),maxdepth)
            elif len(c_node.GetNext()) == 1:
                path.append(c_node.GetNext().keys()[0])
                seq.append(c_node.GetNext().keys()[0].getSeq())
                maxdepth = max(maxdepth,c_node.GetNext().values()[0])
                return find_trail(path,seq,maxdepth)
            else:
                if len(c_node.GetNext()) != 0:
                    print 'Err PDBn.py simplify - find_trail method: tail case is not well observed'
#                 else:
#                     print ''.join(seq),path[-1].getSeq(),path[-1].GetNext(),path[-1].GetLmer(),path[-1].GetRmer()
                return (path,''.join(seq),maxdepth)
        
        def merge_path(path,seq,maxdepth):
            if len(path) <= 1:
                return
            head = path.pop(0)
            tail = path[-1]
            if len(head.GetNext()) != 1:
                print 'Err PDBn.py simplify - merge_path method: ',head.GetRmer(),head.GetLmer()
            head.setSeq(seq)
            for next in tail.GetNext():
                head.addNext(next,maxdepth)
                next.addPrev(head,maxdepth)
            for prev in head.GetPrev():
                prev.addNext(head,maxdepth)
                head.addPrev(prev,maxdepth)
            for node in path:
#                 print node,node.GetRmer(),node.GetLmer()
#                 if node in removed:
#                     print node
#                     print 'fjs'
#                     sys.exit()
#                 else:
#                     removed[node] = True
                self.Rmer_node_map.get(node.GetRmer()).remove(node)
            pass
        
        temp_depth_list = node.GetPrev().values()+node.GetNext().values()
        if len(temp_depth_list) > 2:
            max_depth = max(temp_depth_list)
        elif len(temp_depth_list) == 1:
            max_depth = temp_depth_list[0]
        else:
            max_depth = 0
        
        [path,seq,maxdepth] = find_trail([node],[node.getSeq()],max_depth)            
        merge_path(path,seq,maxdepth)
        return
    
    def simplify(self):

        removed = {}
        source = self.findSourceLnode()
        st = StackNQueue.Queue()
        self.resetVisited()
        [st.enqueue(node) for node in source]
        while st.peek():            
            node = st.dequeue()
            if node.get_visit():
                continue
            node.set_visit()
            self.simplify_node(node)
            
            for nextnode in node.GetNext().keys():
                st.enqueue(nextnode)
        
        pass
    
    
    def writeSpliceGraph(self,filename):
        number_of_node = 0
        for Rnode in self.Rmer_node_map.values():
            for Lnode in Rnode.getLNodeList():
                number_of_node += 1
        s = open(filename,'w')
        s.write('<Database CreatedBy="s3cha@PDBn.py">\n')
        s.write('\t<Gene Name="IG@0" ExonCount="%d" Chromosome="IG" ForwardFlag="1">\n'%(number_of_node))
        source = self.findSourceLnode()
        st = StackNQueue.Queue()
        cycle_test = StackNQueue.Queue()
        self.resetVisited()
        [st.enqueue(node) for node in source]
        index = 0
        last_coordinate = 0
        node_index_map = {}
        while st.peek():
            node = st.dequeue()
            if node in node_index_map:
                continue
            
            testing = False
            for prev in node.GetPrev().keys():
                if prev not in node_index_map:
                    cycle_test.enqueue(node)
                    testing = True
        
            if testing:
                if st.peek() == None and cycle_test.peek() != None:
                    cyc_node = cycle_test.dequeue()
                    for prev in cyc_node.GetPrev().keys()[::-1]:
                        if prev not in node_index_map:
                            print 'Disconnecting cyclic path: ',cyc_node.GetPrev(),cyc_node.GetPrev().get(prev)
                            del cyc_node.GetPrev()[prev]
                            del prev.GetNext()[cyc_node]
                    st.enqueue(cyc_node)
                continue
            
            node_index_map[node] = index
            index += 1
            for next in node.GetNext():
                st.enqueue(next)
            if st.peek() == None and cycle_test.peek() != None:
                cyc_node = cycle_test.dequeue()
                for prev in cyc_node.GetPrev().keys()[::-1]:
                    if prev not in node_index_map:
                        print 'Disconnecting cyclic path: ',cyc_node.GetPrev(),cyc_node.GetPrev().get(prev)
                        del cyc_node.GetPrev()[prev]
                        del prev.GetNext()[cyc_node]
                st.enqueue(cyc_node)

            
        [st.enqueue(node) for node in source]
        def writeNode(node,last_coordinate):
            sequence = node.getSeq()
            exon_string = ['\t\t<Exon Index="%d" Start="%d" End="%d" RefInfo="2">\n'%(node_index_map[node],last_coordinate+1,last_coordinate+1+len(sequence))]
            exon_string.append('\t\t\t<ExonSequence Length="%d">%s</ExonSequence>\n'%(len(sequence),sequence))
            last_coordinate += 1+len(sequence)
            for prev in node.GetPrev():
                exon_string.append('\t\t\t<LinkFrom Index="%d"/>\n'%(node_index_map[prev]))
            exon_string.append('\t\t<Exon>\n')
            s.write(''.join(exon_string))
            return last_coordinate
        
        reverse_index = {}
        for node in node_index_map:
            reverse_index[node_index_map[node]] = node
        
        for index in sorted(reverse_index.keys()):
            node = reverse_index.get(index)
            last_coordinate = writeNode(node,last_coordinate)
        
        
        
#         while st.peek():
#             node = st.dequeue()
#             if node.get_visit():
#                 continue
#             node.set_visit()
#             last_coordinate = writeNode(node,last_coordinate)
#             for next in node.GetNext():
#                 st.enqueue(next)
        s.write('\t</Gene>\n')
        s.write('</Database>')
        s.close()
        pass
        
    
    #reset visited information in all nodes
    def resetVisited(self):
        for node in self.Rmer_node_map.values():
            node.reset_visit()
            for lnode in node.getLNodeList():
                lnode.reset_visit()
        pass
    
    #find the threshold to cut the low depth nodes
    def naiveCut(self):
        number_node = 0
        depth_total = 0
        for node in self.Rmer_node_map.values():
            for lnode in node.getLNodeList():
                number_node += 1
                depth_total += max(sum(lnode.GetPrev().values()),sum(lnode.GetNext().values()))
        threshold = depth_total/number_node
        threshold = threshold * 2/3
        self.naiveCutThreshold(threshold)
        pass
    
    def naiveCutRatio(self,denum):
        number_node = 0
        depth_total = 0
        for node in self.Rmer_node_map.values():
            for lnode in node.getLNodeList():
                number_node += 1
                depth_total += max(sum(lnode.GetPrev().values()),sum(lnode.GetNext().values()))
        threshold = depth_total/number_node
        threshold = max(float(threshold)/denum,1)
        self.naiveCutThreshold(threshold)
        pass
    
    #manually insert the threshold and cut the low depth nodes
    def naiveCutThreshold(self,threshold):
#         print 'Removing the l nodes based on their threshold: ',threshold
        num_del = 0
        total = 0
        for node in self.Rmer_node_map.values():
            for lnode in node.getLNodeList():
                total += 1
                if max(sum(lnode.GetPrev().values()),sum(lnode.GetNext().values())) <= threshold:
                    prev_nodes = lnode.GetPrev()
                    next_nodes = lnode.GetNext()
                    node.remove(lnode)
                    for prev in prev_nodes:
                        if len(prev.GetNext()) == 0:
#                             print prev.getSeq(),prev.GetLmer(),
                            prev.addSeqTail(prev.GetLmer())
#                             print prev.getSeq()
                    for next in next_nodes:
                        if len(next.GetPrev()) == 0:
                            next.addSeqHead(next.GetRmer()[:-1])
                            
                    
                    num_del += 1
#         print 'Number of l node removed: %d over %d'%(num_del,total)
        pass
    
    
    '''clipping the tips in case the read depth of the tip is less than threshold or 1% or 5% of the main stream read depth.
    When clipping the tips, remove the tips only if there is main stream path, and add the read depth to the main stream path.'''
    def tipClip(self,threshold = 3):
        #self.simplify()
        source = self.findSourceLnode()
        while source:
            node = source.pop()
            return_value = self.getMainPath(node,False)
            if return_value == None:
                continue
            m_path,m_cov,c_path,c_cov = return_value 
            
            if node in m_path:
                continue
            elif m_cov == 0: #depth of main path
                print 'test whether this case is possible'
                continue
            elif float(c_cov)/m_cov >= 0.1 or c_cov > threshold:
                continue
            else:
                seq1 = ''.join([t_node.getSeq() for t_node in c_path[:-1]])
                seq2 = ''.join([t_node.getSeq() for t_node in m_path[::-1]])
                dist = Distance.Distance.LevenshteinDistance(seq1,seq2)-math.fabs(len(seq1)-len(seq2))
                if float(dist)/len(seq2) > 0.2:
                    continue
#             print 112
#             print seq1,c_cov
#             print seq2,m_cov
#             print dist,len(seq2),float(dist)/len(seq2)
                        
            #add the node depth to main path
#             print 123, c_cov,m_cov,c_path[-1].GetPrev(),[node.GetPrev().values() for node in [c_path[-1]]+m_path]
            
            for index,c_node in enumerate([c_path[-1]]+m_path[:-1]):
                n_node = m_path[index]
                c_node.GetPrev()[n_node] += c_cov
                n_node.GetNext()[c_node] += c_cov
            #remove the tip
            
            for c_node in c_path[:-1]:
                self.Rmer_node_map.get(c_node.GetRmer()).remove(c_node)
            
            
            
        
        sink = self.findSinkLnode()
        while sink:
            node = sink.pop()
            return_value = self.getMainPath(node,True)
            if return_value == None:
                continue
            m_path,m_cov,c_path,c_cov = return_value
            
            if node in m_path:
                continue
            elif m_cov == 0:
                print 'WARNNING: test whether this case is possible: not intended case'
                continue
            elif float(c_cov)/m_cov >= 0.1 or c_cov > threshold:
                continue
            else:
                seq1 = ''.join([t_node.getSeq() for t_node in c_path[:-1][::-1]])
                seq2 = ''.join([t_node.getSeq() for t_node in m_path])
                dist = Distance.Distance.LevenshteinDistance(seq1,seq2)-math.fabs(len(seq1)-len(seq2))
                if float(dist)/len(seq2) > 0.2:
                    continue
#             print 321, c_cov,m_cov,c_path[-1].GetNext(),[node.GetNext().values() for node in [c_path[-1]]+m_path]
#             seq1 = ''.join([t_node.getSeq() for t_node in c_path[:-1][::-1]])
#             seq2 = ''.join([t_node.getSeq() for t_node in m_path])
#             print 312
#             print dist,float(dist)/len(seq2)
#             print seq1,c_cov
#             print seq2,m_cov
            for index,c_node in enumerate([c_path[-1]]+m_path[:-1]):
                n_node = m_path[index]
                c_node.GetNext()[n_node] += c_cov
                n_node.GetPrev()[c_node] += c_cov
            for c_node in c_path[:-1]:
                self.Rmer_node_map.get(c_node.GetRmer()).remove(c_node)
                                      
#             if node in mainPath:
#                 continue
            ##if next split node has only one input stream for source, or output stream for sink, continue
            ##if the read depth difference of main path and current path is not larger than 0.01 and depth > threshold, then continue.
            
            ##Merge process (Remove the current tip, add the )
            
            
        return
    
    def bulgeRemove(self):
        return
    
    def getDivergentNode(self,node,is_downstream):
        cov = 0
        length = 0
        path = [node]
        if not is_downstream:
            while len(node.GetPrev()) <= 1 and len(node.GetNext()) == 1:
                if not node.GetPrev():
                    length += len(node.getSeq()) - self.length_rmer + 1
                else:
                    length += len(node.getSeq())
                cov = max(cov,node.GetNext().values()[0])
                node = node.GetNext().keys()[0]
                path.append(node)
        else:
            while len(node.GetPrev()) == 1 and len(node.GetNext()) <= 1:
                if not node.GetNext():
                    length += len(node.getSeq()) - len(node.GetRmer())
                else:
                    length += len(node.getSeq())
                cov = max(cov,node.GetPrev().values()[0])
                node = node.GetPrev().keys()[0]
                path.append(node)
        return [node,length,cov,path]
    
    
    def getDominent(self,node,length,is_downstream,path = [],max_depth = 0):
        max_cov = 0
        cur_length = len(node.getSeq())
        if not is_downstream:
            if len(node.GetPrev()) == 0:
                cur_length += 1-len(node.GetRmer())
            for prev_node in node.GetPrev():
                cur_cov = node.GetPrev().get(prev_node)
                if cur_cov > max_cov:
                    n_node = prev_node
                    max_cov = cur_cov
        else:
            if len(node.GetNext()) == 0:
                cur_length += -len(node.GetLmer())
            for next_node in node.GetNext():
                cur_cov = node.GetNext().get(next_node)
                if cur_cov > max_cov:
                    n_node = next_node
                    max_cov = cur_cov
        if max_cov == 0:
            return [path,max_depth]
        elif length < 0:
            return [path,max_depth]
        else:
            if max_depth < max_cov:
                max_depth = max_cov
            path.append(n_node)
            length -= cur_length
            return self.getDominent(n_node,length,is_downstream,path,max_depth)
    
    '''Find the other branches with the highest read depth rather than current node. is_downstream: indicate the path is upstream or downstream'''
    def getMainPath(self,node,is_downstream):
        if not is_downstream: #if error exists in the beginning of the reads
            div_node,length,depth,path = self.getDivergentNode(node,is_downstream)
#             print node.getSeq(),div_node.getSeq()
#             sys.exit()
            n_next_node = None
            if len(div_node.GetPrev()) <= 1:
                if len(div_node.GetNext()) <= 1:
                    return
                min_node = None
                min_cov = sys.maxint
                for n_next_node in div_node.GetNext():
                    cur_cov = div_node.GetNext().get(n_next_node)
                    if min_cov > cur_cov:
                        min_cov = cur_cov
                        min_node = n_next_node
                n_next_node,n_length,n_depth,n_path = self.getDivergentNode(min_node,is_downstream)
                length += n_length
                depth = max(depth,n_depth)
                if len(n_next_node.GetNext()) <= 1:
                    return
                path = path + n_path
                mainPath = self.getDominent(n_next_node,length,is_downstream,[],0)
            else:
                mainPath = self.getDominent(div_node,length,is_downstream,[],0)
        else:
            div_node,length,depth,path = self.getDivergentNode(node,is_downstream)
            n_prev_node = None
            if len(div_node.GetNext()) <= 1:
                if len(div_node.GetPrev()) <= 1:
                    return
                min_node = None
                min_cov = sys.maxint
                for n_prev_node in div_node.GetPrev():
                    cur_cov = div_node.GetPrev().get(n_prev_node)
                    if min_cov > cur_cov:
                        min_cov = cur_cov
                        min_node = n_prev_node
                n_prev_node,n_length,n_depth,n_path = self.getDivergentNode(min_node,is_downstream)
                length += n_length
                depth = max(depth,n_depth)
                if len(n_prev_node.GetPrev()) <= 1:
                    return
                path = path+n_path
                mainPath = self.getDominent(n_prev_node,length,is_downstream,[],0)
            else:
                mainPath = self.getDominent(div_node,length,is_downstream,[],0)
        
        return mainPath+[path,depth]
    
    #remove r-node
    def removeRnode(self,rnode):
        if not rnode.getLNodeList():
            try:
                del self.Rmer_node_map[rnode.GetRmer()]
            except KeyError:
                print 'KeyError: PDBn class, removeRnode method',rnode.GetRmer()
            pass
        else:
            print 'Err: Node not Empty. PDBn class, removeRnode method',rnode.GetRmer()
            
    #return the false negative rate of the graph (compare to the ideal reference should have)
    def error_rate(self,reference,boundary):
        printout = False
        reference = reference[boundary:len(reference)-boundary]
        total_true = len(reference)-self.length_lmer-self.length_rmer
        num_missing = 0
        possible_lmer = []
        for index in range(len(reference)-self.length_lmer-self.length_rmer-1):
            rmerseq = reference[index:index+self.length_rmer]
            lmerseq = reference[index+self.length_rmer:index+self.length_rmer+self.length_lmer]
            nextrmerseq = reference[index+1:index+1+self.length_rmer]
            
            if rmerseq not in self.Rmer_node_map:
                if printout:
                    print 'F: ',rmerseq
                num_missing += 1
                possible_lmer = []
            else:
                if possible_lmer == []:
                    possible_lmer = self.Rmer_node_map.get(rmerseq).getLNodeList()
                check = False
                next_possible_lmer = []
                for lmer in possible_lmer:
                    for next in lmer.GetNext():
                        if next.GetRmer() == nextrmerseq:
                            next_possible_lmer.append(next)
                            check = True
                if not check:
                    if printout:
                        print 'F: ',rmerseq
                    num_missing += 1
                    possible_lmer = []
                else:
                    if printout:
                        print 'T: ',rmerseq
                    possible_lmer = list(set(next_possible_lmer))
        
        return float(num_missing)/total_true
    
    def write(self):
        source = self.FindSourceNode()
        seq = ''
        result = []
        def dfs(node,seq,visited):
            visited[node] = 0
            seq += node.GetRmer()[-1]
            for next_node in node.GetNext():
                if next_node in visited:
                    result.append(seq)
                    continue
                dfs(next_node,seq,visited)
            if not node.GetNext():
                result.append(seq)
            return result
        for node in source:
            for subnode in node.getLNodeList():
                if not subnode.GetPrev():
                    result = dfs(subnode,'',{})
        print result
        self.resetVisited()
        pass
    
    #plot the histogram of the depth of all l-nodes
    def showNodeDepth(self):
        from matplotlib import pyplot as plt
        import numpy as np 
        depth = []
        number_node = 0
#         depth_total = 0
        for node in self.Rmer_node_map:
            for lnode in self.Rmer_node_map.get(node).getLNodeList():
                number_node += 1
                depth.append(max(sum(lnode.GetPrev().values()),sum(lnode.GetNext().values())))
        print 'Number of nodes: ',number_node
        print 'Max read depth: ',max(depth)
        plt.hist(depth,bins = 10**np.linspace(np.log10(1),np.log10(max(depth)),50))
        plt.gca().set_xscale("log")
        plt.gca().set_yscale("log")
        plt.title('Node depth histogram')
#         plt.show()
        pass
    
    #return the number of edges making divergence in the graph
    def divergence(self):
        total_edge = 0
        num_div_edge = 0
        for node in self.Rmer_node_map.values():
            for lnode in node.getLNodeList():
                ledge = len(lnode.GetPrev().keys())
                redge = len(lnode.GetNext().keys())
                total_edge += ledge + redge
                num_div_edge += ledge -1 if ledge >= 1 else 0
                num_div_edge += redge -1 if redge >= 1 else 0
        if total_edge == 0:
            return 0
        else:
            return float(num_div_edge)#/total_edge
        
    def divergenceRef(self,reference):
        total_edge = len(reference)-self.length_lmer-self.length_rmer
        num_div_edge = 0

        for index in range(len(reference)-self.length_lmer-self.length_rmer-1):
            rmerseq = reference[index:index+self.length_rmer]
            lmerseq = reference[index+self.length_rmer:index+self.length_rmer+self.length_lmer]
            nextrmerseq = reference[index+1:index+1+self.length_rmer]
            nextlmerseq = reference[index+1+self.length_rmer:index+1+self.length_rmer+self.length_lmer]
            
            if rmerseq not in self.Rmer_node_map:
                continue
            for lnode in self.Rmer_node_map.get(rmerseq).getLNodeList():
                if Distance.Distance.HammingDistance(lmerseq,lnode.GetLmer()) <=1:
                    num_div_edge += len(lnode.GetPrev())-1 if len(lnode.GetPrev()) != 0 else 0
                    num_div_edge += len(lnode.GetNext())-1 if len(lnode.GetNext()) != 0 else 0
        return float(num_div_edge)/total_edge
    


if __name__ == '__main__':
    rmer = 10
    lmer = 20
#     x = SdbConstruction(rmer,lmer)
#     x.quick_test(x)

if __name__ == '__main__0':
    lmer = 6
    rmer = 5
    
    pdb = PDBn(lmer,rmer)
#     pdb.__initialize__(lmer, rmer)
    distance = Distance.Distance()
    read_class = ReadClass.SingleEnd()
#     pdb.AddRead(read_class,'SFSAAFFAASSAADDAABAAAAAHFIEJFJWIHAFJWHFIWJFJAKSDNCMOQWIFJ')
#     pdb.AddRead(read_class,'BAAAAAHFIEJFJWIHRFJWHFRWJFJAKSDNCMOQWIFJ')
#     pdb.AddRead(read_class,'BAAAAAHFIEJFJWIHAFJWHFIWJFJAKSDNCMOQWI')
#     pdb.AddRead(read_class,'BAAAAAHFIEJFJMIMAFJWHFIWJFJAKSDNCMOQGKFJ')
#     pdb.AddRead(read_class,'DKKSSJJFFLLSSJAHABAAAAAHFIEKFJWIMAFJWHFIWJFJAKSDNCMOQGKFJ')
#     pdb.AddRead(read_class,'BAAAAAHFIEKFJWIMAFJWHFIWJFJAKSDNCMOQWIKJ')
    pdb.AddRead(read_class,'PPPQXQPPWWWQQWWAABBAARRGTTATTEEFEAA')
    pdb.AddRead(read_class,'PPPQQQPPWWWQQWWAABBAARRGTTATTXEFEAA')
    pdb.AddRead(read_class,'PPPQQQPPWWWQQWWAABBAXRRGTTATTEEFEAA')
    pdb.AddRead(read_class,'PPPQQQPPWWXQQWWAABBAARRGTTATTEEFEAA')
#     pdb.AddRead(read_class,'BBAARRGTQATTEEFEAAendoftheline')
#     pdb.AddRead(read_class,'ATTEEFEAAendoftheline')
    pdb.AddRead(read_class,'ATTEEFEAAendoftheline')
    pdb.AddRead(read_class,'OOIIUXUUYYUUJJKAABBAARRGGAATTEEFFAA')
    pdb.AddRead(read_class,'OOIIUUUUYYUUJXKAABBAARRGGAATTEEFFAA')
    pdb.AddRead(read_class,'OOIIUUUUYYUUJJKAABBAXRRGGAATTEEFFAA')
    pdb.AddRead(read_class,'OOIIUUUUYYUUJJKAABBAARRGGAATTEEXFAA')
#     pdb.AddRead(read_class,'DKKSSJJFFLLSSJAHABAAAAAHFIEKFJWIMAFJWHFIWJFJAKSDNCMOQGKFJ')
#     pdb.AddRead(read_class,'BAAAAAHFIEKFJWIMAFJWHFIWJFJAKSDNCMOQWIKJ')
    pdb.NodeClustering()
    pdb.naiveCutThreshold(1)
    pdb.ShowGraph()
    print pdb.error_rate('OOIIUUUUYYUUJJKAABBAARRGGAATTEEFFAA', 0)
    print pdb.divergence()
    print pdb.divergenceRef('OOIIUUUUYYUUJJKAABBAARRGGAATTEEFFAA')
