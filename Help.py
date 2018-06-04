#!/usr/bin/env python

## (c) 2017, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import sys

from Utils import bold, bright

class Help():
    shortstrings = {}
    longstrings = {}
    longest = 0
    keynames = []

    def addHelp(self, keys, argp, shorthelp, longhelp):
        self.keynames.append(keys[0])
        for k in keys:
            self.shortstrings[k] = (keys, argp, shorthelp)
            self.longstrings[k] = longhelp
            if len(k) > self.longest:
                self.longest = len(k)

    def __init__(self):
        self.addHelp(["-i", "--fasta"], True, "Input file in FASTA format.", """
Name of the input file containing aligned sequences, in FASTA format. The first sequence
in the file is assumed to be the reference sequence, unless a different file is specified
with the -r option, in which case the reference sequence is read from that file. If this
option is not provided, sequences are read from standard input. All sequences in this file
should have the same length as the reference sequence.""")

        self.addHelp(["-r", "--ref", "--reference"], True, "File containing reference sequence in FASTA format.", """
If this option is specified, the reference sequence will be read from this file (in FASTA
format) rather than from the main input file.""")

        self.addHelp(["-s", "--site", "--sites"], True, "Sites to detect.", """
The value of this option should be one or more nucleotide strings describing the methylated
sites the program should detect. For example: "-s CG GCH". The first C in each site is 
assumed to be the potentially-converted base. Default: CG.""")

        self.addHelp(["-o", "--open"], True, "Number of methylated sites required to open a patch.", "")
        self.addHelp(["-c", "--close"], True, "Number of unmethylated sites required to close a patch.", "")
        self.addHelp(["-x", "--strand"], True, "Strand to be examined (one of t, b, tb, bt).", "")
        self.addHelp(["-w", "--weights"], True, "Weights for C positions and patches (a list of 5 numbers).", """
The value of this option should be a comma-separated list of 5 numeric values, representing the weight of 
open C position, open patch, undetermined, closed patch, closed C position respectively. Default: 2,1,0,-1,2.""")

        self.addHelp(["--map"], True, "Name of map output file.", "")
        self.addHelp(["--csv"], True, "Name of tab-delimited output file.", "")
        self.addHelp(["-C", "--cluster-on"], True, "Map(s) to perform clustering on.", "")
        self.addHelp(["-p", "--cluster-from"], True, "Start position of region for clustering.", "")
        self.addHelp(["-q", "--cluster-to"], True, "End position of region for clustering.", "")
        self.addHelp(["-g", "--cluster-dist"], True, "Distance metric to use for clustering.", "")
        self.addHelp(["-m", "--cluster-meth"], True, "Clustering method (see cluster3 docs).", "")
        self.addHelp(["--cluster-path"], True, "Path to the cluster3 executable.", "")

    def shortHelp(self):
        sys.stderr.write("metyhlmapper.py - Generate and plot methylation maps.\n\nOptions:\n\n")
        fstr = " {:" + str(self.longest + 4) + "}"
        for key in self.keynames:
            (keys, argp, shorthelp) = self.shortstrings[key]
            n = len(keys)
            nm1 = n - 1
            for i in range(n):
                k = keys[i]
                if i == nm1:
                    sys.stderr.write(fstr.format(k + (bright(" ___", 1) if argp else "")) + "| " + shorthelp + "\n")
                else:
                    sys.stderr.write(fstr.format(k + (bright(" ___", 1) if argp else "")) + "|\n")

    def shortHelp2(self):
        sys.stderr.write(bold("metyhlmapper.py") + " - Generate and plot methylation maps.\n\nOptions:\n\n")
        for key in self.keynames:
            (keys, argp, shorthelp) = self.shortstrings[key]
            sys.stderr.write(" " + ", ".join([k + (bright(" ___", 1) if argp else "") for k in keys]) + "\n   " + shorthelp + "\n\n")

    def longHelp(self, key):
        if key in self.shortstrings:
            (keys, argp, shorthelp) = self.shortstrings[key]
            sys.stderr.write("  " + key + bright(" ___", 1) if argp else "" + "\n")
            sys.stderr.write(self.longstrings[key] + "\n\n")
