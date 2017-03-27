'''
Created on Nov 15, 2016

@author: s3cha
'''
import PDBn
from Distance import Distance
import LNode
import math
from string import digits

class RNode(object):
    '''
    classdocs
    '''
    def __init__(self, _rmerseq,pdbn):
        '''
        Constructor
        '''
        self.pdbn = pdbn
        self.LNodes = {}
        self.Rmer_seq = _rmerseq
        self.visited = False
        
    def GetRmer(self):
        return self.Rmer_seq
        
    #visited status setting for the visiting function
    def set_visit(self):
        self.visited = True
        pass
    
    def reset_visit(self):
        self.visited = False
        pass
    
    def get_visit(self):
        return self.visited
    
    #make a distance table for each of l-nodes in this current r-node.
    def distanceTable(self):
        lnodes = [item for item in self.LNodes.keys() if len(item) == self.pdbn.get_Llength()]
        Table = [[] for x in xrange(len(lnodes))]
        for index1,lnode1 in enumerate(lnodes):
            for index2,lnode2 in enumerate(lnodes):
                if index2 < index1:
                    Table[index1].append(Table[index2][index1])
                    continue
                elif index1 == index2:
                    Table[index1].append(0)
                else:
                    Table[index1].append(self.getLNode(lnode1).distance(self.getLNode(lnode2)))
        return [Table,[self.LNodes[node_instance] for node_instance in lnodes]]
    
    
    def CombineTwo(self,node1,node2):
        for other_node in node2.GetPrev():
            if other_node in node1.GetPrev():
                node1.GetPrev()[other_node] += node2.GetPrev().get(other_node)
                other_node.GetNext()[node1] = node1.GetPrev().get(other_node) 
            else:
                node1.GetPrev()[other_node] = node2.GetPrev().get(other_node)
                other_node.GetNext()[node1] = node1.GetPrev().get(other_node)
#             try:
#                 del other_node.GetNext()[node2]
#             except KeyError:
#                 print 'ERR RNode class, 1.CombineTwo function: CombineTwo L-node doesn\'t have correct key pair.',node1,node2
#                 pass
        for other_node in node2.GetNext():
            if other_node in node1.GetNext():
                node1.GetNext()[other_node] += node2.GetNext().get(other_node)
                other_node.GetPrev()[node1] = node1.GetNext()[other_node] 
            else:
                node1.GetNext()[other_node] = node2.GetNext().get(other_node)
                other_node.GetPrev()[node1] = node1.GetNext()[other_node]
#             try:
#                 del other_node.GetPrev()[node2]
#             except KeyError:
#                 print 'ERR RNode class, 2.CombineTwo function: CombineTwo L-node doesn\'t have correct key pair.',node1,node2
#                 pass
        try:
            self.remove(node2) #del self.LNodes[node2.GetLmer()]
            del node2
        except KeyError:
            print 'ERR RNode class, 3.CombineTwo function: L-node not exist.',node1,node2
            pass
        pass
    
    #Combine l-node in the list of listLnode param. sort them with the length of l-mer seq (first key), and read depth (second key)
    def CombineLNode(self,listLnode):
        listLnode = sorted(listLnode,reverse = True)
        for lnode in listLnode:
            if lnode == listLnode[0]:
                continue
            else:
                self.CombineTwo(listLnode[0],lnode)
        pass
    
    def compare_distance(self,item1,item2):
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
    
    '''Merge the l-nodes if their difference is less than the error tolerance using the following order.
    1. merge the l-nodes if their previous l-node is same
    2. build distance table of remaining l nodes where distance table = {key = distance: value = (l node1, l node2, weight difference, max weight of the node)}
    3. merge the l node in a order of (a) smallest key, (b) largest node weight, (c) largest weight difference 
    >> from the most largest weight node merge the smallest weight node if distance is less than threshold'''
    def SimplifyLnodes_method1(self,table,lnodes):
        node_depth = []
        for lnode in lnodes:
            node_depth.append(max(sum(lnode.GetPrev().values()),sum(lnode.GetNext().values())))
        distance = {}
        for index1 in range(len(table)):
            for index2 in range(index1+1,len(table)):
                dif = table[index1][index2]
                if dif <= self.pdbn.get_error_tolerance():
                    if dif not in distance:
                        distance[dif] = []
                    distance[dif].append((index1,index2,node_depth[index1]-node_depth[index2],max(node_depth[index1],node_depth[index2])))
            
        merged_index = []
        for key in sorted(distance.keys()):
            distance[key] = sorted(distance[key],cmp=self.compare_distance)
            for value in distance[key]:
                main_index = value[0] if value[2]>0 else value[1]
                sub_index  = value[1] if value[2]>0 else value[0]
                if set((main_index,sub_index)).intersection(merged_index):# in merged_index:
                    continue
                merged_index.append(sub_index)
                self.CombineTwo(lnodes[main_index],lnodes[sub_index])
#         if distance:
#             print 111,distance
        pass
    
    '''Merge the l-nodes if one of their l node length is shorter than regular length, and full node has no previous node and sequence is exactly matched.'''
    def SimplifyLnodes_boundarycase(self):
        full_node = []
        cut_node = []
        for lnode in self.getLNodeList():
            if len(lnode.GetLmer().translate(None,digits)) == self.pdbn.get_Llength():
                full_node.append(lnode)
            else:
                cut_node.append(lnode)
        full_node.sort()
        cut_node.sort()
        for main_node in full_node:
            if not cut_node:
                break
            if not main_node.GetPrev():
                for sub_node in cut_node[::-1]:
                    if Distance.HammingDistance_dif(main_node.GetLmer(), sub_node.GetLmer().translate(None,digits)) == 0:
                        self.CombineTwo(main_node,sub_node)
                        cut_node.remove(sub_node)
        pass
    
    '''merge similar L-nodes based on the similarity of their sequence (clustering based on the distance matrix)
    call distanceTable and CombineLNode functions'''
    def MergeLNodes(self):
        #merge the l nodes if their previous l node are identical
        test_prev = {}
        def setPrev(prev_lnode,lnode,test_prev):
            if prev_lnode not in test_prev:
                test_prev[prev_lnode] = [lnode]
            else:
                test_prev[prev_lnode].append(lnode)
            return test_prev
        merge_node_list = []
        for lnode in self.LNodes.values():
            if len(lnode.GetPrev())>1:
                if len(lnode.GetLmer()) < self.pdbn.get_Llength():
                    for index,prev_lnode in enumerate(lnode.GetPrev()):
                        new_lnode = LNode.LNode(lnode.GetLmer()+str(index),self.Rmer_seq)
                        new_lnode.SplitLnode(lnode,prev_lnode)
                        self.LNodes[lnode.GetLmer()+str(index)] = new_lnode
                        new_lnode.setSeq(lnode.getSeq())
                        test_prev = setPrev(prev_lnode,new_lnode,test_prev)
                    self.remove(lnode)
                else:    
                    merge_node_list.append(lnode.GetPrev().keys())
                    for prev_lnode in lnode.GetPrev():
                        test_prev = setPrev(prev_lnode,lnode,test_prev)
            else:
                for prev_lnode in lnode.GetPrev():
                    test_prev = setPrev(prev_lnode,lnode,test_prev)
        if merge_node_list:
            for index1 in range(len(merge_node_list))[::-1]:
                for index2 in range(index1):
                    if len(set(merge_node_list[index1]).intersection(merge_node_list[index2])) != 0:
                        
                        merge_node_list[index2] = merge_node_list[index2] + merge_node_list[index1]
                        merge_node_list.pop(index1)
                        break
        for merge in merge_node_list:
            merge = set(merge)
            first = None
            for other in merge:
                if first == None:
                    first = other
                    continue
                test_prev[first] += test_prev[other]
                del test_prev[other]
            test_prev[first] = list(set(test_prev[first]))
        
        for prev in test_prev:
            if len(test_prev.get(prev)) > 1:
                self.CombineLNode(test_prev.get(prev))
        #make distance matrix of all remaining l nodes
        [Table,lnodes] = self.distanceTable()
        self.SimplifyLnodes_method1(Table,lnodes)
        self.SimplifyLnodes_boundarycase()
        pass
    
    #return corresponding Lnode based on the lmer sequence. create one if not exist.
    def getLNode(self,_lmerseq):
        Lnode = self.LNodes.get(_lmerseq)
        if Lnode == None:
            Lnode = LNode.LNode(_lmerseq,self.Rmer_seq)
            self.LNodes[_lmerseq] = Lnode
        return Lnode
    
    def peekLNode(self,_lmerseq):
        node = self.LNodes.get(_lmerseq)
        if node == None:
            for lnode in self.LNodes.values():
                if lnode.GetLmer() == _lmerseq:
                    return lnode
        return node
#     def getLNodeS(self,_lmerseq):
#         if len(_lmerseq) < self.pdbn.get_Llength():
#             Lnode = LNode.LNode(_lmerseq,self.Rmer_seq)
#             self
        

    def getLNodeList(self):
        return self.LNodes.values()
    
    def remove(self,lnode):
        lnode.remove()
        try:
            for lmer_seq in self.LNodes:
                
                if self.LNodes[lmer_seq] == lnode:
                    del self.LNodes[lmer_seq]
                    break
            del lnode
        except KeyError:
            print 'KeyError: RNode class, remove method'
            pass
        if not self.LNodes:
            self.pdbn.removeRnode(self)
        pass
        
        
if __name__ == '__main__':
    x = 1
        