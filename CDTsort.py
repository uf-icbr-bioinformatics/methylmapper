#!/usr/bin/env python

## (c) 2017, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import sys
import csv
import numpy as np

class Node():
    name = ""
    left = ""
    right = ""
    simil = ""
    values = []
    mean = 0.0

    def __init__(self, n, l, r, s):
        self.name = n
        self.left = l
        self.right = r
        self.simil = s
        self.values = []
        self.mean = 0.0

class GTree():
    nodes = {}
    root = ""

    def __init__(self, gtrfile):
        self.nodes = {}
        node = None
        with open(gtrfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                node = Node(line[0], line[1], line[2], line[3])
                self.nodes[line[0]] = node
        self.root = node.name
        sys.stderr.write("{} nodes read, root={}\n".format(len(self.nodes), self.root))

    def fillTree(self, CDT):
        self.fillTreeAux(self.nodes[self.root], CDT)

    def fillTreeAux(self, node, CDT):
        l = node.left
        r = node.right
        if l.startswith("GENE"):
            leftavg = CDT.data[l]
            node.values.append(leftavg)
        else:
            leftnode = self.nodes[l]
            self.fillTreeAux(leftnode, CDT)
            for x in leftnode.values:
                node.values.append(x)
        if r.startswith("GENE"):
            rightavg = CDT.data[r]
            node.values.append(rightavg)
        else:
            rightnode = self.nodes[r]
            self.fillTreeAux(rightnode, CDT)
            for x in rightnode.values:
                node.values.append(x)
        node.mean = np.mean(node.values)

    def sortTree(self, CDT):
        self.sortTreeAux(self.nodes[self.root], CDT)
        self.writeTree()

    def sortTreeAux(self, node, CDT):
        l = node.left
        r = node.right

        if l.startswith("NODE"):
            leftnode = self.nodes[l]
            self.sortTreeAux(leftnode, CDT)
            leftavg = leftnode.mean
        else:
            leftavg = CDT.data[l]

        if r.startswith("NODE"):
            rightnode = self.nodes[r]
            self.sortTreeAux(rightnode, CDT)
            rightavg = rightnode.mean
        else:
            rightavg = CDT.data[r]

        if leftavg > rightavg:
            tmp = node.left
            node.left = node.right
            node.right = tmp

    def writeTree(self):
        self.writeTreeAux(self.nodes[self.root])
        
    def writeTreeAux(self, node):
        l = node.left
        r = node.right

        if l.startswith("NODE"):
            self.writeTreeAux(self.nodes[l])
        if r.startswith("NODE"):
            self.writeTreeAux(self.nodes[r])
        sys.stdout.write("{}\t{}\t{}\t{}\n".format(node.name, l, r, node.simil))

    def dumpTree(self, CDT):
        self.dumpTreeAux(self.nodes[self.root], CDT, 0)

    def dumpTreeAux(self, node, CDT, indent):
        l = node.left
        r = node.right
        if l.startswith("GENE"):
            sys.stdout.write(" "*indent + l + " avg={}\n".format(CDT.data[l]))
        else:
            leftnode = self.nodes[l]
            self.dumpTreeAux(leftnode, CDT, indent+1)
        if r.startswith("GENE"):
            sys.stdout.write(" "*indent + r + " avg={}\n".format(CDT.data[r]))
        else:
            rightnode = self.nodes[r]
            self.dumpTreeAux(rightnode, CDT, indent+1)
        sys.stdout.write(" "*indent + "{} avg={} values={}\n".format(node.name, node.mean, node.values))
        

class CDTFile():
    data = {}

    def __init__(self, cdtfile):
        self.data = {}
        with open(cdtfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            c.next()
            c.next()            # Skip two header lines
            for line in c:
                gene = line[0]
                avg = np.mean([float(x) for x in line[4:]])
                self.data[gene] = avg
        sys.stderr.write("Averages read for {} genes.\n".format(len(self.data)))

if __name__ == "__main__":
#    sys.setrecursionlimit(5000) # Use this if running out of stack
    G = GTree(sys.argv[1])
    C = CDTFile(sys.argv[2])
    G.fillTree(C)
    sys.stderr.write("Tree filled.\n")
    #G.dumpTree(C)
    G.sortTree(C)
