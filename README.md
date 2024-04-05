# Methylmapper

Methylmapper is a program to produce and display methylation maps from FASTA files.

## Requirements

This program requires:

* The `gdcreate` program from the [gdprogs](https://github.com/albertoriva/gdprogs) repository. Please ensure that `gdcreate` is in PATH, otherwise use the GDCREATE_PATH variable in `bin/methylmapper` to specify its location.
* The `cluster3` program. If it is not in PATH, please use the CLUSTER3_PATH variable in `bin/methylmapper` to specify its location.

## Usage

Call the bin/methylmapper script with the appropriate options. The following table describes all available command-line options. `___` indicates that the option takes an argument.

Option | Description
--- | ---
*Input options*
 -i ___, --fasta ___ |   Input file in FASTA format (required). This will typically contain short reads in which unmethylated Cs are converted to Ts.
 -r ___, --ref ___, --reference ___ | File containing reference sequence in FASTA format (required).
 -s ___, --site ___, --sites ___ | Sites to detect (default: CG). Can be specified more than once.
 -d ___ |    Subsample: read only this number of reads (at random) from the input file.
 -u |    Remove duplicate input sequences.
 -U |    Remove duplicate input sequences (by pattern).
 *Map options*
 -o ___, --open ___  | Number of methylated sites required to open a patch (default: 2).
 -c ___, --close ___ | Number of unmethylated sites required to close a patch (default: 1)
 -x ___, --strand ___ | Strand to be examined (one of t, b, tb, bt) (default: t).
 -n ___, --unconv ___ |    Maximum number of consecutive unconverted Cs (default: no limit).
*Clustering options*
 -w ___, --weights ___ |    Weights for C positions and patches (a list of 5 numbers).
 -C ___, --cluster-on ___ |    Map(s) to perform clustering on.
 -W ___, --cluster-weights ___ |    Weights of the site maps during clustering.
 -p ___, --cluster-from ___ |    Start position of region for clustering.
 -q ___, --cluster-to ___ |    End position of region for clustering.
 -g ___, --cluster-dist ___ |    Distance metric to use for clustering (see cluster3 docs).
 -m ___, --cluster-meth ___ |    Clustering method (see cluster3 docs).
*Output options*
 --map ___ |    Name of map output file.
 --csv ___ |    Name of tab-delimited output file.
 --plot ___ |    Name of heatmap output file.
 -z |    Display gaps and Ns as white in heatmap.

## Algorithm
1. Each read is converted into one or more *maps* as follows:
  * The program locates all occurrences of the specified site (e.g. CG) and determines whether the are methylated (*) or unmethylated (#).
  * Methylated sites are grouped into *patches*. Two consecutive methylated sites generate a new patch, while a single unmethylated site closes it (these numbers can be changed with `-o` and `-c`).
