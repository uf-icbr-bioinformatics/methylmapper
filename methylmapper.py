#!/usr/bin/env python

## (c) 2017, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import sys
import random
import hashlib
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator

import Help
import MethMap
import Cluster
import RefSequence
from Utils import makeColHeaders, INPUT, OUTPUT, WARNING, BANNER, MAPS

# CG -> red black, GC -> yellow black

class MethylMapper():
    filename   = None
    reffile    = None
    refseq     = None
    sites      = []
    sequences  = []
    references = []
    maps       = []
    sampleseqs = None           # Number of sequences to retain from input using random sampling
    maxnamelen = 9              # Length of longest sequence name (at least as long as "Reference")
    remdups    = False          # If true, remove duplicate sequences (-u option)
    
    # Output files
    mapfile  = None
    csvfile  = None
    freqfile = None
    plotfile = None

    # Map parameters
    top = True                  # Look for sites on top strand?
    bottom = False              # Look for sites on bottom strand?
    openMin = 2                 # Number of sites to open a patch
    closeMin = 1                # Number of sites to close a patch
    weights = [2.0, 1.0, 0.0, -1.0, -2.0]

    # Clustering
    clust = Cluster.Clusterer()

    def getStrands(self):
        if self.top:
            if self.bottom:
                return "both"
            else:
                return "top"
        else:
            return "bottom"

    def readSequences(self, f):
        seen = set()
        ns = 0
        removed = 0
        for rec in FastaIterator(f):
            if self.refseq is None:
                self.refseq = rec
            else:
                good = True
                if self.remdups:
                    md5 = hashlib.md5(str(rec.seq)).hexdigest()
                    if md5 in seen:
                        good = False
                        removed += 1
                    else:
                        seen.add(md5)
                if good:
                    self.sequences.append(rec)
                    self.maxnamelen = max(self.maxnamelen, len(rec.name))
                    ns += 1
        if self.remdups:
            sys.stderr.write(INPUT + "{} duplicate sequence(s) removed.\n".format(removed))
        if self.sampleseqs and self.sampleseqs < ns:
            indices = range(ns)
            random.shuffle(indices)
            newseqs = [ self.sequences[w] for w in indices[:self.sampleseqs] ]
            self.sequences = newseqs

    def initialize(self):
        if self.reffile:
            with open(self.reffile, "r") as f:
                for rec in FastaIterator(f):
                    self.refseq = rec
                    break

        if self.filename:
            sys.stderr.write(INPUT + "Reading sequences from file `{}'.\n".format(self.filename))
            with open(self.filename, "r") as f:
                self.readSequences(f)
        else:
            sys.stderr.write(INPUT + "Reading sequences from standard input.\n")
            self.readSequences(sys.stdin)

        sys.stderr.write(INPUT + "Reference sequence: {}bp.\n".format(len(self.refseq)))
        sys.stderr.write(INPUT + "{} input sequences.\n".format(len(self.sequences)))
        sys.stderr.write(INPUT + "Detected sites: " + ", ".join(self.sites) + ".\n")
        sys.stderr.write(INPUT + "Detection strands: " + self.getStrands() + "\n")
        sys.stderr.write(INPUT + "Open/close: {}/{}\n".format(self.openMin, self.closeMin))
        sys.stderr.write(INPUT + "Weights: {}\n".format(self.weights))

        self.references = []
        self.maps       = []
        for site in self.sites:
            mref = RefSequence.RefSequence(self.refseq, site)
            mmap = MethMap.MethMap(site, mref, weights=self.weights)
            mmap.openMin = self.openMin
            mmap.closeMin = self.closeMin
            mmap.top = self.top
            mmap.bottom = self.bottom
            self.references.append(mref)
            self.maps.append(mmap)

    def generateMaps(self):
        for m in self.maps:
            m.makeAllMaps(self.sequences)
            if self.freqfile:
                outfile = m.site + "-" + self.freqfile
                sys.stderr.write(MAPS + "Saving {} frequencies to {}.\n".format(m.site, outfile))
                m.calcAllFrequencies(self.sequences, outfile)

    ### Output

    def writeMapsText(self):
        sys.stderr.write(OUTPUT + "Writing maps in text format to file {}\n".format(self.mapfile))
        fstr = "{:" + str(self.maxnamelen) + "} |{}\n"
        with open(self.mapfile, "w") as out:
            for m in self.maps:
                out.write("## " + m.site + "\n")
                out.write(fstr.format("Reference", m.ref.sequence))
                out.write(fstr.format("Sites", m.sitesToString()))
                for row in m.mapstrings:
                    name = row[0]
                    seq = "".join(row[1])
                    out.write(fstr.format(name, seq))
                out.write("\n")

    def writeMapsCSV(self):
        sys.stderr.write(OUTPUT + "Writing maps in CSV format to files:\n")
        ncols = len(self.refseq)
        hdr = makeColHeaders(ncols)
        hdrline = "#Seq\t" + "\t".join(hdr) + "\n"
        for m in self.maps:
            m.writeCSV(self.csvfile, hdrline)

    ### Top level

    def helpWanted(self, args):
        if len(args) == 0:
            return (True, [])
        for i in range(len(args)):
            if args[i] in ["-h", "--help"]:
                return (True, args[i+1:])
        return (False, [])

    def parseArgs(self, args):
        (helpwanted, what) = self.helpWanted(args)
        if helpwanted:
            h = Help.Help()
            if what:
                for w in what:
                    h.longHelp(w)
            else:
                h.shortHelp2()
            return False

        valuedArgs = ["-i", "--fasta", "-r", "--ref", "--reference", "-o", "--open", "-c", "--close", "-s", "--site", "--sites", "--map", "--csv",
                      "-f", "--freq", "-C", "--cluster-on", "-p", "--cluster-from", "-q", "--cluster-to", "-g", "--cluster-dist", "-m", "--cluster-meth",
                      "--cluster-path", "--plot", "-x", "--strand", "-w", "--weights", "-d"]
        next = ""
        for a in args:
            if next in ["-i", "--fasta"]:
                self.filename = a
                next = ""
            elif next in ["-r", "--ref", "--reference"]:
                self.filename = a
                next = ""
            elif next in ["-d"]:
                self.sampleseqs = int(a)
                next = ""
            elif next in ["-o", "--open"]:
                self.openMin = int(a)
                next = ""
            elif next in ["-c", "--close"]:
                self.closeMin = int(a)
                next = ""
            elif next in ["-s", "--site", "--sites"]:
                if a[0] == "-":
                    next = a
                else:
                    self.sites.append(a)
            elif next in ["-f", "--freq"]:
                self.freqfile = a
                next = ""
            elif next in ["--map"]:
                self.mapfile = a
                next = ""
            elif next in ["--csv"]:
                self.csvfile = a
                next = ""
            elif next in ["-C", "--cluster-on"]:
                if a[0] == "-":
                    next = a
                else:
                    self.clust.clusterOn.append(a)
            elif next in ["-p", "--cluster-from"]:
                self.clust.clusterFrom = int(a) - 1
                next = ""
            elif next in ["-q", "--cluster-to"]:
                self.clust.clusterTo = int(a) - 1
                next = ""
            elif next in ["-g", "--cluster-dist"]:
                self.clust.clusterDist = a
                next = ""
            elif next in ["-m", "--cluster-meth"]:
                self.clust.clusterMeth = a
                next = ""
            elif next == "--cluster-path":
                self.clust.clusterPath = a
                next = ""
            elif next == "--plot":
                self.plotfile = a
                next = ""
            elif next in ["-x", "--strand"]:
                if a == "t":
                    self.top = True
                    self.bottom = False
                elif a == "b":
                    self.top = False
                    self.bottom = True
                elif a == "tb" or a == "bt":
                    self.top = True
                    self.bottom = True
                next = ""
            elif next in ["-w", "--weights"]:
                newweights = self.parseWeights(a)
                self.weights = newweights
                next = ""
            elif a in valuedArgs:
                next = a
            elif a == '-u':
                self.remdups = True
            else:
                sys.stderr.write(WARNING + "Unknown command-line option `{}'.\n".format(a))

        if not self.sites:
            self.sites = ["CG"]
        if not self.clust.clusterOn:
            self.clust.clusterOn = [self.sites[0]]

        good = []
        bad = []
        for c in self.clust.clusterOn:
            if c in self.sites:
                good.append(c)
            else:
                bad.append(c)
        if bad:
            sys.stderr.write(WARNING + "Unknown sites for clustering: {}.\n".format(", ".join(bad)))
        self.clust.clusterOn = good
        return True

    def parseWeights(self, w):
        good = False
        try:
            new = [float(x) for x in w.split(",")]
            good = True
        except:
            pass
        if good:
            return new
        else:
            sys.stderr.write(WARNING + "Weights should be specified as a list of 5 numbers separated by commas, e.g. 2,1,0,-1,-2. Argument ignored.\n")
            return False

    def main(self):
        self.generateMaps()
        if self.mapfile:
            self.writeMapsText()
        if self.csvfile:
            self.writeMapsCSV()
        if self.clust.clusterOn:
            self.clust.run(self.maps, plotfile=self.plotfile)

### Main

if __name__ == "__main__":
    M = MethylMapper()
    if M.parseArgs(sys.argv[1:]):
        sys.stderr.write(BANNER)
        M.initialize()
        M.main()
