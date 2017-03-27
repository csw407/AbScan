'''
Created on Nov 21, 2016

@author: s3cha
'''

class LinkedNode(object):
    def __init__(self,value):
        self.value = value
        self.next_node = None 
        pass
    
    def next(self):
        return self.next_node
    
    def set_next(self,node):
        self.next_node = node
        pass
    
    def get_value(self):
        return self.value
    

class Stack(object):
    def __init__(self):
        self.head = None
        pass
    
    def push(self,value):
        if self.head == None:
            self.head = LinkedNode(value)
        else:
            new_head = LinkedNode(value)
            new_head.set_next(self.head)
            self.head = new_head
        pass
    
    def pop(self):
        value = self.head.get_value()
        self.head = self.head.next()
        return value
    
    def peek(self):
        if self.head == None:
            return None
        else:
            return self.head.get_value()

class Queue(object):
    def __init__(self):
        self.head = None
        self.tail = None
        pass
    
    def enqueue(self,value):
        if self.head == None:
            node = LinkedNode(value)
            self.head = node
            self.tail = node
        else:
            node = LinkedNode(value)
            self.tail.set_next(node)
            self.tail = node
    
    def dequeue(self):
        if self.head == None:
            return None
        node = self.head
        if node.next() == None:
            self.head = None
            self.tail = None
        else:
            self.head = self.head.next()
        return node.get_value()
    
    def peek(self):
        if self.head == None:
            return None
        else:
            return self.head.get_value()
        
    def printout(self):
        node = self.head
        while True:
            if node == None:
                break
            print node.get_value()
            node = node.next()
        pass