#!/usr/bin/env python

## (c) 2017, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import os
import sys
import os.path
import tempfile
import subprocess

import Draw
from Utils import saferm, makeColHeaders, INPUT, OUTPUT, WARNING, CLUSTER

class Clusterer():
    clusterOn   = []
    clusterFrom = 0             # Start position of region used for clustering
    clusterTo   = None          # End position of region used for clustering
    clusterDist = "7"
    clusterMeth = "m"
    clusterPath = "cluster3"    # Path to the cluster3 executable
    
    def run(self, maps, plotfile=None):
        wantedMaps = []
        for site in self.clusterOn:
            for m in maps:
                if m.site == site:
                    wantedMaps.append(m)
                    break
        # print [m.site for m in wantedMaps]
        m0 = wantedMaps[0]
        if not self.clusterTo:
            self.clusterTo = m0.ref.length
        ncols = self.clusterTo - self.clusterFrom
        if ncols < 0:
            sys.stderr.write(WARNING + "Clustering region [{}, {}] is empty!\n".format(self.clusterFrom + 1, self.clusterTo))
            return False
        sys.stderr.write(CLUSTER + "Clustering on {} maps, region=[{}, {}]\n".format(len(wantedMaps), self.clusterFrom + 1, self.clusterTo))
        totcols = len(wantedMaps) * ncols
        hdr = makeColHeaders(totcols)
        tmpfile = tempfile.mkstemp(dir=".")[1]
        saferm(tmpfile)
        csvfile = tmpfile + ".csv"
        cdtfile = tmpfile + ".cdt"
        gtrfile = tmpfile + ".gtr"
        try:
            with open(csvfile, "w") as out:
                out.write("#Sequence\t" + "\t".join(hdr) + "\n")
                for (name, fmap) in m0.mapstrings:
                    out.write(name)
                    for m in wantedMaps:
                        #vect = m.mapvectors[name]     # *** THIS SHOULD BE DECIDED BY THE scale FLAG!
                        vect = m.sclvectors[name]
                        for i in range(self.clusterFrom, self.clusterTo):
                            out.write("\t" + str(vect[i]))
                    out.write("\n")

            cmd = [self.clusterPath, "-f", csvfile, "-g", self.clusterDist, "-m", self.clusterMeth]
            sys.stderr.write(CLUSTER + "Executing: " + " ".join(cmd) + "\n")
            retcode = subprocess.call(cmd)
            if retcode != 0:
                sys.stderr.write(WARNING + "cluster3 command returned exit code {}!\n".format(retcode))
                return False
            if not(os.path.isfile(cdtfile) and os.path.isfile(gtrfile)):
                sys.stderr.write(WARNING + "Clustering failed - check cluster3 command line.\n")
                return False

            sys.stderr.write(CLUSTER + "Clustering successful.\n")
            sys.stderr.write(CLUSTER + "Writing CDT files:\n")

            rownames = {}
            roworder = []
            with open(cdtfile, "r") as f:
                f.readline()
                f.readline()
                for line in f:
                    fields = line.split("\t")
                    rownames[fields[1]] = fields[0]
                    roworder.append(fields[1])
            for m in maps:
                if m.csvfile:
                    m.writeCDT(rownames, roworder, gtrfile)
                else:
                    sys.stderr.write(WARNING + "Writing CDT file requires the --csv option.\n")
                # m.dump()
            if plotfile:
                sys.stderr.write(CLUSTER + "Saving heatmap to: {}.\n".format(plotfile))
                Draw.plotMap(plotfile, maps)
            return True
        finally:
            saferm(csvfile)
            saferm(cdtfile)
            saferm(gtrfile)
