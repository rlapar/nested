import sys
import pprint

class IntervalTree:
    def __init__(self, intervalList, transposonList):
        self.intervalList = intervalList
        self.transposonList = transposonList
        self.__buildTree()

    def __buildTree(self):
        self.root = TreeNode([-3000000000,3000000000], None)
        for i in reversed(range(len(self.intervalList))):
            self.root.addNode(self.intervalList[i], self.transposonList[i])

    def toDict(self):
        return self.root.toDict()

class TreeNode:
    def __init__(self, interval, transposon):
        self.children = []
        self.interval = interval
        self.parent = None
        self.transposon = transposon

    def contains(self, value):
        return value >= self.interval[0] and value <= self.interval[1]

    def setParent(self, node):
        self.parent = node

    def addNode(self, interval, transposon):
        for child in self.children:
            if child.contains(interval[0]):
                child.addNode(interval, transposon)
                return

        node = TreeNode(interval, transposon)
        node.setParent(self)
        self.children.append(node) 

    def toDict(self):
        node = {
            'interval': self.interval,
            'transposon': self.transposon,
            'children': []
        }
        for child in self.children:
            node['children'].append(child.toDict())
        return node
