'''
Created on Nov 15, 2016

@author: s3cha
'''
import sys

from Distance import Distance


class LNode(object):
    '''
    classdocs
    '''


    def __init__(self, _lmer,_rmer):
        '''
        Constructor
        '''
        self.Lmer_seq = _lmer
        self.Rmer_seq = _rmer
        self.prev_nodes = {}
        self.next_nodes = {}
        self.seq = ['NA']
        self.visited = False
        
    def get_visit(self):
        return self.visited
    
    def set_visit(self):
        self.visited = True
        pass
    
    def reset_visit(self):
        self.visited = False
    
    #less than method for sorting
    def __lt__(self,other):
        if len(self.GetLmer()) < len(other.GetLmer()):
            return True
        else:
            return max(sum(self.GetPrev().values()),sum(self.GetNext().values())) < max(sum(other.GetPrev().values()),sum(other.GetNext().values()))
    
    def replacePrev(self,origin,alter,depth_value):
        if self.prev_nodes.has_key(origin) and origin.prev_nodes.has_key(alter):
            self.prev_nodes[alter] = depth_value
        else:
            print 'ERR in LNode class, replacePrev method: there is no edges between self to origin or origin to alter'
        pass
     
    def replaceNext(self,origin,alter,depth_value):
        if self.next_nodes.has_key(origin) and origin.next_nodes.has_key(alter):
            self.next_nodes[alter] = depth_value
        else:
            print 'ERR in LNode class, replaceNext method: there is no edges between self to origin or origin to alter'
        pass
    
    #connect previous edge to current lnode.
    def ConnectPrev(self,_prev):
        if _prev not in self.prev_nodes:
            self.prev_nodes[_prev] = 1
        else:
            self.prev_nodes[_prev] += 1
        pass
    
    #connect next edge to current lnode.
    def ConnectNext(self,_next):
        if _next not in self.next_nodes:
            self.next_nodes[_next] = 1
        else:
            self.next_nodes[_next] += 1
        pass
    
    def SplitLnode(self,other,prev):
        self.prev_nodes[prev] = other.GetPrev().get(prev)
        prev.GetNext()[self] = self.prev_nodes[prev]
        for next_lnode in other.GetNext():
            self.next_nodes[next_lnode] = other.GetNext().get(next_lnode)
            next_lnode.GetPrev()[self] = self.next_nodes[next_lnode]
        pass
    
    def addNext(self,nextnode,depth):
        self.next_nodes[nextnode] = depth
        pass
    
    def addPrev(self,prevnode,depth):
        self.prev_nodes[prevnode] = depth
        pass
    
    def resetSeq(self):
        self.seq = []
        pass
    
    def setSeq(self,_seq):
        self.seq = [_seq]
        pass
    
    def addSeqHead(self,_seq):
        self.seq.insert(0,_seq)
        pass
    
    def addSeqTail(self,_seq):
        self.seq.append(_seq)
        pass
    
    def getSeq(self):
        return ''.join(self.seq)
    
    def GetNext(self):
        return self.next_nodes
    
    def GetPrev(self):
        return self.prev_nodes
    
    def GetRmer(self):
        return self.Rmer_seq
    
    def GetLmer(self):
        return self.Lmer_seq
    
    def distance(self,_lnode):
        return Distance.HammingDistance(self.GetLmer(), _lnode.GetLmer())
    
    def remove(self):
        for prev in self.prev_nodes:
            try:
                del prev.GetNext()[self]
            except KeyError:
                print 'Err: LNode class, remove method. Prev l node have no edge to current node',self.Lmer_seq,self.Rmer_seq
                sys.exit()
                pass
        for next in self.next_nodes:
            try:
                del next.GetPrev()[self]
            except KeyError:
                print 'Err: LNode class, remove method. Next l node have no edge to current node',self.Lmer_seq,self.Rmer_seq
        pass
    
    def printnode(self):
        current = self
        prev_path = []
        while len(current.prev_nodes) == 1:
            prev_path.insert(0,current)
            current = current.prev_nodes.keys()[0]
        next_path = []
        current = self
        while len(current.next_nodes) == 1:
            next_path.append(current)
            current = current.next_nodes.keys()[0]
        for node in prev_path:
            print node.GetRmer(),node.GetLmer(),node.getSeq()
        print 'current'
        for index,node in enumerate(next_path):
            print node.GetRmer(),node.GetLmer(),node.getSeq()
            
    