#!/usr/bin/env python

## (c) 2017, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import sys
import colormaps
from Creator import Creator
from Drawer import Drawer

class Drawable():
    width = 0
    height = 0
    xoffset = 0
    yoffset = 0
    margin = 5

    def drawLine(self, d, x1, y1, x2, y2, color):
        d.drawLine(x1 + self.xoffset, y1 + self.yoffset, x2 + self.xoffset, y2 + self.yoffset, color)

    def drawRectangle(self, d, x1, y1, x2, y2, color):
        d.drawRectangle(x1 + self.xoffset, y1 + self.yoffset, x2 + self.xoffset, y2 + self.yoffset, color)

    def drawFilledRectangle(self, d, x1, y1, x2, y2, color):
        d.drawFilledRectangle(x1 + self.xoffset, y1 + self.yoffset, x2 + self.xoffset, y2 + self.yoffset, color)

    def drawPolygon(self, d, points, color):
        newpoints = []
        for i in range(0, len(points), 2):
            newpoints.append(points[i] + self.xoffset)
            newpoints.append(points[i+1] + self.yoffset)
        d.drawPolygon(newpoints, color)

    def drawFilledPolygon(self, d, points, color):
        newpoints = []
        for i in range(0, len(points), 2):
            newpoints.append(points[i] + self.xoffset)
            newpoints.append(points[i+1] + self.yoffset)
        d.drawFilledPolygon(newpoints, color)

class ClustPanel(Drawable):
    rowh = 10
    cellw = 10
    nrow = 0                    # Rows in CDT file
    ncol = 0                    # Columns in CDT file
    data = []

    def __init__(self, cdtfile, rowh=10, cellw=10):
        self.rowh = rowh
        self.margin = rowh / 2
        self.cellw = cellw
        self.data = []

        with open(cdtfile, "r") as f:
            hdr = f.readline().split("\t")
            self.ncol = len(hdr) - 4
            f.readline()
            for line in f:
                self.nrow += 1
                self.data.append( [ (x if x == "." else int(float(x))) for x in line.rstrip("\r\n").split("\t")[4:] ] )
        self.width = self.cellw * self.ncol + self.margin * 2
        self.height = self.rowh * self.nrow + self.margin * 2

    def findLimits(self, line):
        """Returns the index of the first and last non-zero position in `line'."""
        a = 0
        b = len(line) - 1
        for x in range(a, b):
            if line[x] != 0:
                a = x
                break
        for x in range(b, a - 1, -1):
            if line[x] != 0:
                b = x
                break
        return (a, b)

    def draw(self, d, cmap):
        row = 0
        for line in self.data:
            col = 0
            (a, b) = self.findLimits(line)
            for v in line:
                if a <= col <= b:
                    x = col * self.cellw + self.margin
                    y = row * self.rowh + self.margin
                    # print v
                    color = cmap.getColor(v)
                    # print "{},{} {} => {}".format(x, y, v, color)
                    # raw_input()
                    self.drawFilledRectangle(d, x, y, x + self.cellw - 1, y + self.rowh - 2, color)
                col += 1
            row += 1

class SiteBar(Drawable):
    ncol = 0
    cellw = 10
    trianglesize = 5
    positions = []

    def __init__(self, ncol, positions, cellw=10, trianglesize=5):
        """ncol should be taken from the corresponding ClustPanel."""
        self.ncol = ncol
        self.positions = positions
        self.cellw = cellw
        self.trianglesize = trianglesize
        self.width = self.cellw * self.ncol + self.margin * 2
        self.height = 20

    def draw(self, d):
        y = self.height - self.margin
        self.drawLine(d, self.margin, y, self.width - self.margin, y, 1)
        for p in self.positions:
            x = self.margin + p * self.cellw + (self.cellw / 2)
            self.drawFilledPolygon(d, [x, y, x-self.trianglesize, y-(self.trianglesize*2), x+self.trianglesize+1, y-(self.trianglesize*2)], 1)

class ClustTree(Drawable):
    treewidth = 100
    rowh = 10
    spacing = 6
    genes = []
    names = []
    branches = []
    coords = {}
    maxnamelen = 0

    def __init__(self, cdtfile, gtrfile, treewidth=100, rowh=10):
        self.treewidth = treewidth
        self.rowh = rowh
        self.margin = rowh / 2
        self.genes = []
        self.names = []
        self.branches = []
        self.coords = {}
        
        row = 0
        ypos = 0
        self.writeNames = (self.rowh > 12)

        with open(cdtfile, "r") as f:
            f.readline()
            f.readline()
            for line in f:
                fields = line.split("\t")
                gname = fields[0]
                name = fields[1]
                self.genes.append(gname)
                self.names.append(name)
                self.maxnamelen = max(self.maxnamelen, len(name))
                ypos = self.margin + row * self.rowh
                self.coords[gname] = (self.margin + self.treewidth, ypos)
                row += 1
        self.height = ypos + self.margin
        if self.writeNames:
            self.width = self.treewidth + self.spacing + self.maxnamelen * 6 + self.margin * 2
        else:
            self.width = self.treewidth + self.margin * 2

        with open(gtrfile, "r") as f:
            for line in f:
                fields = line.rstrip("\r\n").split()
                fields[3] = float(fields[3])
                self.branches.append(fields)

    def draw(self, d):

        if self.writeNames:
            for gidx in range(len(self.genes)):
                gene = self.genes[gidx]
                name = self.names[gidx]
                d.drawString(name, self.xoffset + self.margin + self.treewidth + self.spacing, self.yoffset + self.coords[gene][1], 1, anchor=4)

        while True:
            found = None
            for gidx in range(len(self.genes) - 1):
                g1 = self.genes[gidx]
                g2 = self.genes[gidx+1]
                br = self.findBranch(g1, g2) # find branch containing g1 and g2
                if br:
                    found = True
                    break

            if not found:
                sys.stderr.write("Error: branch for {} and {} not found!\n".format(g1, g2))
                return
            top = br[0]
            left = br[1]
            leftcoord = self.coords[left]
            right = br[2]
            rightcoord = self.coords[right]
            # print leftcoord
            # print rightcoord
            blen = self.margin + int(br[3] * self.treewidth)
            self.coords[top] = (blen, (leftcoord[1] + rightcoord[1]) / 2)
            # print self.coords[top]
            self.drawLine(d, blen, leftcoord[1], blen, rightcoord[1], 1)
            self.drawLine(d, blen, leftcoord[1], leftcoord[0], leftcoord[1], 1)
            self.drawLine(d, blen, rightcoord[1], rightcoord[0], rightcoord[1], 1)

            self.genes[gidx:gidx+2] = [br[0]]
            # print self.genes
            if len(self.genes) == 1:
                break

    def findBranch(self, g1, g2):
        for br in self.branches:
            if ((br[1] == g1) and (br[2] == g2)) or ((br[1] == g2) and (br[2] == g1)):
                return br
        return None

class MultiClusterPlot():
    maps = []
    
    def __init__(self, maps):
        self.maps = maps

    def draw(self):
        m0 = self.maps[0]
        cdt = m0.cdtfile
        gtr = m0.gtrfile

class MapColors():
    name = ""
    rgbs = []
    colors = []

    def __init__(self, name, rgbs):
        self.name = name
        self.rgbs = rgbs
        self.colors = {}

    def allocate(self, cm, weights):
        n = 0
        idx = cm.allocate(".", 255, 255, 255)
        self.colors["."] = "."
        for rgb in self.rgbs:
            cname = '{}{}'.format(self.name, n)
            idx = cm.allocate(cname, *rgb)
            # print "{} => {}".format(n, idx)
            self.colors[weights[n]] = cname
            n += 1
        #print "Map {}: {}".format(self.name, self.colors)

    def getColor(self, idx):
        if idx in self.colors:
            return self.colors[idx]
        else:
            return [255, 255, 255]

def makeColormaps(cm, weights):
    weights.reverse()
    cmaps = {"rgb": MapColors("rgb", [[250, 250, 250],
                                      [  0,   0,   0],
                                      [128, 128, 128],
                                      [255,   0,   0],
                                      [128,   0,   0]]),
             "ygb": MapColors("ygb", [[250, 250, 250],
                                      [  0,   0,   0],
                                      [128, 128, 128],
                                      [254, 254,   0],
                                      [125, 125,   0]]) }
    for mc in cmaps.values():
        mc.allocate(cm, weights)
    #print cm.cmap
    #print cm.colors
    return cmaps

def plotMap(plotfile, methmaps, rowh=15, cellw=3):
    panels  = []
    bars    = []
    map0    = methmaps[0]

    tree  = ClustTree(map0.cdtfile, map0.gtrfile, rowh=rowh)
    totwidth  = tree.width
    for mmap in methmaps:
        pan = ClustPanel(mmap.cdtfile, rowh=rowh, cellw=cellw)
        panels.append(pan)
        totwidth += pan.width
        bars.append(SiteBar(pan.ncol, mmap.allPositions(), cellw=cellw))

    totheight = panels[0].height + bars[0].height

    c = Creator(pathname="gdcreate")
    d = Drawer(c)
    d.createImage(totwidth, totheight)
    cm = colormaps.StandardColorMap(d)
    cmaps = makeColormaps(cm, map0.weights)
    d.setColormap(cm)
    d.setFont(2)

    tree.yoffset = panels[0].margin + bars[0].height
    tree.draw(d)
    xo = tree.width
    yo = bars[0].height
    i = 0
    for mmap in methmaps:
        panels[i].xoffset = xo
        panels[i].yoffset = yo
        cmapname = 'ygb' if mmap.site.startswith("GC") else 'rgb'
        panels[i].draw(d, cmap=cmaps[cmapname])
        bars[i].xoffset = xo
        bars[i].draw(d)
        xo += panels[i].width
        i += 1

    d.saveImage(plotfile)
    d.close()
