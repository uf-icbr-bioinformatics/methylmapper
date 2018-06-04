#!/usr/bin/env python

## (c) 2017, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import sys
import shutil

from Utils import parseLine, OUTPUT, CLUSTER

### Utils

### Classes 

class Block():
    char = ""
    start = 0
    end = 0
    n = 0
    prev = None
    next = None

    def __init__(self, char, start, prev):
        self.char = char
        self.start = start
        self.end = start
        self.n = 1
        self.prev = prev

    def __str__(self):
        return "<Block {}{}{}-{}>".format(self.n, self.char, self.start, self.end)

class BaseFreq():
    pos = 0
    n = 0
    counts = {}

    def __init__(self, pos):
        self.pos = pos
        self.n = 0
        self.counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    def base(self, b):
        b = b.upper()
        if b in "ACGT":
            self.n += 1
            self.counts[b] += 1

    def writeRow(self, out):
        out.write("{}\t{}\t{}\t{}\t{}\n".format(self.pos, 
                                                1.0 * self.counts['A'] / self.n,
                                                1.0 * self.counts['C'] / self.n,
                                                1.0 * self.counts['G'] / self.n,
                                                1.0 * self.counts['T'] / self.n))

class MethMap():
    site = ""
    ref = None
    mapstrings = []
    mapvectors = {}
    sclvectors = {}
    sitefreqs  = {}    # frequencies for C positions in sites
    otherfreqs = {}    # frequencies for other C positions
    positions  = []
    weights    = [2.0, 1.0, 0.0, -1.0, -2.0]
    charvalues = {'*': 2.0, '+': 1.0, ' ': 0.0, '-': -1.0, '#': -2.0}
    scale      = True  # If true, generate scaled vectors

    # These are copied from the MethylMapper object
    openMin  = 2
    closeMin = 1
    top      = True
    bottom   = True

    # Files
    csvfile = None              # added by writeMapCSV
    cdtfile = None              # added by Clusterer
    gtrfile = None              # added by Clusterer

    def __init__(self, site, ref, weights=None):
        self.site = site 
        self.ref  = ref
        self.mapstrings = []
        self.mapvectors = {}
        self.sitefreqs  = {}
        self.otherfreqs = {}
        self.positions  = []
        for p in ref.cpositionsTop:
            self.sitefreqs[p] = BaseFreq(p)
            self.positions.append(p)
        for p in ref.cpositionsBot:
            self.sitefreqs[p] = BaseFreq(p)
            self.positions.append(p)
        self.positions.sort()

        # print "othercTop:", ref.othercTop
        for p in ref.othercTop:
            self.otherfreqs[p] = BaseFreq(p)
        for p in ref.othercBot:
            self.otherfreqs[p] = BaseFreq(p)
        if weights:
            self.weights = weights
            self.setCharvalues(weights)

    def setCharvalues(self, values):
        self.charvalues['*'] = values[0]
        self.charvalues['+'] = values[1]
        self.charvalues[' '] = values[2]
        self.charvalues['-'] = values[3]
        self.charvalues['#'] = values[4]

    def dump(self, s=sys.stdout):
        s.write("""Map for: {}
Sites: {}
""".format(self.site, self.positions))

    def allPositions(self):
        if self.top and self.bottom:
            return self.positions
        elif self.top:
            return self.ref.cpositionsTop
        else:
            return self.ref.cpositionsBot

    def sitesToString(self):
        s = [' ']*len(self.ref.sequence)
        for p in self.ref.cpositionsTop:
            s[p] = '+'
        for p in self.ref.cpositionsBot:
            s[p] = '-'
        return "".join(s)

    def makeBlocks(self, s):
        blocks = []
        current = None
        for i in range(len(s)):
            ch = s[i]
            if ch in ['*', '#']:
                if current and ch == current.char:
                    current.end = i
                    current.n += 1
                else:
                    b = Block(ch, i, current)
                    blocks.append(b)
                    current = b

        # Fix forward pointers
        for j in range(len(blocks) - 1):
            blocks[j].next = blocks[j+1]

        return blocks

    def findBlock(self, blocks, char, n):
        """Find a block in `blocks' with at least `n' occurrences of `char'."""
        for b in blocks:
            if b.char == char and b.n >= n:
                return b
        return None

    def findOpeningBlock(self, block, openChar, closeChar, best=None):
        """Follow the prev pointers in the blocks list, returning the left-most
openChar block with a number of sites greater than openMin. A closeChar block
with nsites >= closeMin terminates the search."""
        #print "*** Extend left: {}, best={}".format(block, best)
        if not block:
            return best
        if block.char == closeChar and block.n >= self.closeMin:
            return best
        if block.char == openChar and block.n >= self.openMin:
            best = block
        return self.findClosingBlock(block.prev, openChar, closeChar, best)

    def findClosingBlock(self, block, openChar, closeChar, best=None):
        """Follow the next pointers in the blocks list, returning the right-most 
openChar block with a number of sites greater than openMin. A closeChar block with
nsites >= closeMin terminates the search."""
        #print block
        if not block:
            #print "no more"
            return best
        if block.char == closeChar and block.n >= self.closeMin:
            #print "closing"
            return best
        if block.char == openChar and block.n >= self.openMin:
            #print "setting best to {}".format(block)
            best = block
        return self.findClosingBlock(block.next, openChar, closeChar, best)

    def fillMapStringSingle(self, s, blocks, openChar, closeChar):
        regStart = 0
        regEnd = 0
        regions = []
        fillChar = '+' if openChar == '*' else '-'
        result = list(s)

        for b in blocks:
            if b.char == openChar and b.n >= self.openMin:
                regStart = b.start
                regEnd   = b.end
                c = self.findOpeningBlock(b.prev, openChar, closeChar)
                if c:
                    # print "extended left to {}".format(c)
                    regStart = c.start
                c = self.findClosingBlock(b.next, openChar, closeChar)
                if c:
                    # print "extended right to {}".format(c)
                    regEnd = c.end
                for i in range(regStart, regEnd+1):
                    if result[i] == ' ':
                        result[i] = fillChar
        return "".join(result)

    def fillMapString(self, s, blocks):
        top = self.fillMapStringSingle(s, blocks, '*', '#')
        bot = self.fillMapStringSingle(s, blocks, '#', '*')

        nconflicts = 0
        result = []
        for i in range(len(top)):
            if top[i] == bot[i]:
                result.append(top[i])
            elif top[i] == ' ':
                result.append(bot[i])
            elif bot[i] == ' ':
                result.append(top[i])
            else:
                nconflicts += 1
                result.append(top[i])

        vector = [ self.charvalues[c] for c in result ]
        return (result, vector)
                
    def scaleVector(self, fmap, vect):
        patchStart = None
        patchType = "N"
        newVect = [x for x in vect]

        for i in range(len(fmap)):
            # print (i, fmap[i], patchType, patchStart)
            if patchType in '+-':
                if fmap[i] == patchType:
                    continue
                else:
                    plen = i - patchStart
                    for y in range(patchStart, i):
                        newVect[y] = float(newVect[y]) / plen
                    patchType = "N"
            else:
                if fmap[i] in '+-':
                    patchType = fmap[i]
                    patchStart = i
        return newVect

    def makeAllMaps(self, sequences):
        for seq in sequences:
            basemap = self.ref.makeMapString(str(seq.seq), top=self.top, bottom=self.bottom)
            # print basemap
            blocks = self.makeBlocks(basemap)
            # print [str(b) for b in blocks]
            (fmap, vect) = self.fillMapString(basemap, blocks)
            # print fmap
            # raw_input()
            # print vect
            # raw_input()
            # print "".join(fmap)
            # raw_input()
            self.mapstrings.append((seq.name, fmap))
            self.mapvectors[seq.name] = vect
            self.sclvectors[seq.name] = self.scaleVector(fmap, vect)

    def calcAllFrequencies(self, sequences, freqfile):
        for seq in sequences:
            for p in self.ref.cpositionsTop:
                self.sitefreqs[p].base(seq[p])
            for p in self.ref.cpositionsBot:
                self.sitefreqs[p].base(seq[p])
            # print sorted(self.otherfreqs.keys())
            for p in self.ref.othercTop:
                self.otherfreqs[p].base(seq[p])
            for p in self.ref.othercBot:
                self.otherfreqs[p].base(seq[p])

        with open(freqfile, "w") as out:
            out.write("# Site Cs, Top\n")
            out.write("Pos\tA\tC\tG\tT\n")
            for p in self.ref.cpositionsTop:
                self.sitefreqs[p].writeRow(out)
            out.write("\n# Site Cs, Bot\n")
            out.write("Pos\tA\tC\tG\tT\n")
            for p in self.ref.cpositionsBot:
                self.sitefreqs[p].writeRow(out)
            out.write("\n# Non-site Cs, Top\n")
            out.write("Pos\tA\tC\tG\tT\n")
            for p in self.ref.othercTop:
                self.otherfreqs[p].writeRow(out)
            out.write("\n# Non-site Cs, Bot\n")
            out.write("Pos\tA\tC\tG\tT\n")
            for p in self.ref.othercBot:
                self.otherfreqs[p].writeRow(out)

    def writeCSV(self, csvname, hdrline):
        self.csvfile = "{}-{}".format(self.site, csvname)
        sys.stderr.write(OUTPUT + "  " + self.csvfile + "\n")
        with open(self.csvfile, "w") as out:
            out.write(hdrline)
            for (name, data) in self.mapvectors.iteritems():
                out.write(name + "\t" + "\t".join(str(x) for x in data) + "\n")
        if self.scale:
            self.sclfile = "{}-scaled.csv".format(self.site)
            with open(self.sclfile, "w") as out:
                out.write(hdrline)
                for (name, data) in self.sclvectors.iteritems():
                    out.write(name + "\t" + "\t".join(str(x) for x in data) + "\n")

    def writeCDT(self, rownames, roworder, gtrfile):
        """Write a CDT file for this map, using the gene order specified in `gtrfile'."""
        self.cdtfile = self.site + "-map.cdt"
        self.gtrfile = self.site + "-map.gtr"
        sys.stderr.write(CLUSTER + "  {} ({})\n".format(self.cdtfile, self.gtrfile))
        shutil.copyfile(gtrfile, self.gtrfile)
        with open(self.cdtfile, "w") as out:
            with open(self.csvfile, "r") as f:
                hdr = parseLine(f.readline())
                newhdr = ["GID", hdr[0], "NAME", "GWEIGHT"] + hdr[1:]
                out.write("\t".join(newhdr) + "\n")
                newhdr = "EWEIGHT\t\t\t1.000000\t" + "\t".join(["1.000000"]*(len(hdr) - 1)) + "\n"
                out.write(newhdr)

                oldlines = {}   # read original file into dictionary
                for line in f:
                    parsed = parseLine(line)
                    name = parsed[0]
                    oldlines[parsed[0]] = parsed

                # Write rows out in new order
                for name in roworder:
                    parsed = oldlines[name]
                    out.write("\t".join([rownames[name], name, name, "1.000000"] + parsed[1:]) + "\n")
