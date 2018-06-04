#!/usr/bin/env python

## (c) 2017, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import os

### Some ANSI fun...

def color(s, c):
    """Return string `s' in color `c' (0..7)."""
    return u"\u001b[" + str(c+30) + "m" + s + u"\u001b[0m"

def bright(s, c):
    """Return string `s' in bright color `c' (0..7)."""
    return u"\u001b[" + str(c+30) + ";1m" + s + u"\u001b[0m"

def bold(s):
    """Return string `s' in bold."""
    return u"\u001b[1m" + s + u"\u001b[0m"

BANNER  = bold(bright("[#######] methylmapper.py\n[#######] (c) 2017, A. Riva, ICBR Bioinformatics Core, University of Florida\n", 5))
INPUT   =       bold("[input  ] ")
OUTPUT  =       bold("[output ] ")
MAPS    =       bold("[maps   ] ")
CLUSTER =       bold("[cluster] ")
WARNING = bold(color("[warning] ", 1))

def makeColHeaders(n):
    """Returns a list of n strings of the form C1, C2... Cn, to use as column headers."""
    return [ "C" + str(x) for x in range(1, n+1) ]

def parseLine(s):
    return s.strip("\r\n").split("\t")

def saferm(filename):
    try:
        os.remove(filename)
    except:
        pass
