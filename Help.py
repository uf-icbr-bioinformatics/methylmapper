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
        self.addHelp(["-C", "--cluster-on"], True, "Map(s) to perform clustering on.", """
The value of this option should be one or more nucleotide strings chosen from the ones listed in the 
-s option. For example, if the value of -s is `CG GC', possible values for this option are `CG', `GC', 
or `CG GC'. Data from the specified maps will be used when clustering.""")
        self.addHelp(["-W", "--cluster-weights"], True, "Weights of the site maps during clustering.", """
The value of this option should be a list of numbers (integers or floating point) describing the weight
to be assigned to each map used in clustering. For example, if -C is `CG GC' and -W is `2 1', then the
CG map will have a weight for clustering double that of GC. The provided values are automatically scaled
to 1.0.""")
        self.addHelp(["-p", "--cluster-from"], True, "Start position of region for clustering.", "")
        self.addHelp(["-q", "--cluster-to"], True, "End position of region for clustering.", "")
        self.addHelp(["-g", "--cluster-dist"], True, "Distance metric to use for clustering.", "")
        self.addHelp(["-m", "--cluster-meth"], True, "Clustering method (see cluster3 docs).", "")
        self.addHelp(["--cluster-path"], True, "Path to the cluster3 executable.", "")
        self.addHelp(["-d"], True, "Read only this number of reads (at random) from the input file.", "")
        self.addHelp(["-u"], False, "Remove duplicate input sequences.", """
If supplied, sequences from the input file that are identical to already seen ones will be discarded.
Useful to remove PCR artifacts.""")
        self.addHelp(["-U"], False, "Remove duplicate input sequences (by pattern).", """
If supplied, sequences showing a methylation pattern identical to an already seen one will be discarded.
Useful to display unique methylation patterns only.""")
        self.addHelp(["-n", "--unconv"], True, "Maximum number of consecutive unconverted Cs.", """
If supplied, sequences containing this number of consecutive unconverted Cs or more will be discarded
before starting the analysis.""")
        self.addHelp(["--plot"], True, "Name of heatmap output file.", "")
        self.addHelp(["-z"], False, "Display gaps and Ns as white in heatmap.", "")

    def shortHelp(self):
        sys.stdout.write("metyhlmapper.py - Generate and plot methylation maps.\n\nOptions:\n\n")
        fstr = " {:" + str(self.longest + 4) + "}"
        for key in self.keynames:
            (keys, argp, shorthelp) = self.shortstrings[key]
            n = len(keys)
            nm1 = n - 1
            for i in range(n):
                k = keys[i]
                if i == nm1:
                    sys.stdout.write(fstr.format(k + (bright(" ___", 1) if argp else "")) + "| " + shorthelp + "\n")
                else:
                    sys.stdout.write(fstr.format(k + (bright(" ___", 1) if argp else "")) + "|\n")

    def shortHelp2(self):
        sys.stdout.write(bold("metyhlmapper.py") + " - Generate and plot methylation maps.\n\nOptions:\n\n")
        for key in self.keynames:
            (keys, argp, shorthelp) = self.shortstrings[key]
            sys.stdout.write(" " + ", ".join([k + (bright(" ___", 1) if argp else "") for k in keys]) + "\n   " + shorthelp + "\n\n")

    def longHelp(self, key):
        if key in self.shortstrings:
            (keys, argp, shorthelp) = self.shortstrings[key]
            sys.stdout.write("  " + key + bright(" ___", 1) if argp else "" + "\n")
            sys.stdout.write(self.longstrings[key] + "\n\n")
