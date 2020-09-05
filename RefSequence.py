## (c) 2017, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import sys

import Bio
import Bio.Seq
import Bio.SeqUtils

from Utils import INPUT, OUTPUT, WARNING, MAPS

class RefSequence():
    sequence = ""
    target = ""
    length = 0
    coffset = 0
    positionsTop =  []          # list of site positions on top strand
    cpositionsTop = []          # position of C nucleotide in each site (top)
    positionsBot =  []          # list of site position on bottom strand
    cpositionsBot = []          # position of C nucleotide in each site (bottom)
    othercTop     = []          # list of other C positions (top)
    othercBot     = []          # list of other C positions (bottom)

    def __init__(self, sequence, target):
        tg = str(target)
        self.sequence = str(sequence.seq)
        self.length = len(self.sequence)
        self.target = target
        self.othercTop = []
        self.othercBot = []
        if len(target) > 0:
            if "C" in tg:
                coffset = target.index("C")
                hits = Bio.SeqUtils.nt_search(self.sequence, tg)
                self.positionsTop = hits[1:]
                self.cpositionsTop = [h + coffset for h in self.positionsTop]

                # Now do bottom strand
                target = Bio.Seq.reverse_complement(target)
                coffset = len(target) - coffset - 1
                hits = Bio.SeqUtils.nt_search(self.sequence, str(target))
                self.positionsBot = hits[1:]
                self.cpositionsBot = [h + coffset for h in self.positionsBot]
                for i in range(len(self.sequence)):
                    b = self.sequence[i]
                    if b == 'C' and i not in self.cpositionsTop:
                        self.othercTop.append(i)
                    elif b == 'G' and i not in self.cpositionsBot:
                        self.othercBot.append(i)
                #print self.positionsTop
                #print self.cpositionsTop
                #print self.positionsBot
                #print self.cpositionsBot
                #print self.othercTop
                #print self.othercBot
                #raw_input()
                sys.stderr.write(MAPS + "{} map: {} sites ({} top, {} bot)\n".format(tg, len(self.positionsTop) + len(self.positionsBot),
                                                                                     len(self.positionsTop), len(self.positionsBot)))
                sys.stderr.write(MAPS + "        {} non-site C positions ({} top, {} bot)\n".format(len(self.othercTop) + len(self.othercBot),
                                                                                                    len(self.othercTop), len(self.othercBot)))
            else:
                sys.stderr.write(WARNING + "target `{}' does not contain a C.\n".format(tg))

    def makeMapString(self, read, top=True, bottom=True):
        if len(read) != self.length:
            sys.stderr.write("Warning: read length ({}) does not match reference sequence length ({}).\n".format(len(read), self.length))
            return None
        smap = [" "]*self.length
        pattern = []
        if top:
            for i in self.cpositionsTop:
                if read[i] == 'C':
                    smap[i] = '*'
                    pattern.append("0")
                else:
                    smap[i] = '#'
                    pattern.append("1")
        if bottom:
            for i in self.cpositionsBot:
                if read[i] == 'G':
                    smap[i] = '*'
                    pattern.append("0")
                else:
                    smap[i] = '#'
                    pattern.append("1")
        # print pattern
        # print (self.cpositionsTop, self.cpositionsBot)
        # print (len(self.cpositionsTop), len(self.cpositionsBot), len(pattern))
        return ("".join(smap), "".join(pattern))

    def methylStretch(self, read, maxunconv, top=True, bottom=True):
        """Returns True if `read' contains `maxunconv' or more non-converted Cs, otherwise False."""
        mu = 0
        for i in self.cpositionsTop:
            if read[i] == 'C':
                mu += 1
                if mu == maxunconv:
                    return True
            else:
                mu = 0
        return False
