Scripts used in the manuscript "Evaluating and improving the representation of 
bacterial contents in long-read metagenome assemblies".

### Description

There is one ready-to-use set of scripts (bin merging), 
one demo script (k-mer spectra plot functions), 
and supporting scritps.
Python needs to be >=3.6 for 
formatted string literals. Scripts were used in a python=3.10.5 environment.

#### bin_merger.py: combine the results of circle rescue and traditional binner  

Given a hifiasm-meta assembly & circle rescue results and genome binning 
from a traditional binner (e.g. metaBAT2),
this script tries to merge the two without introducing too much duplication. 

Script has three required inputs: 1) a tab-delimited plain text file
containing contig lengths and names, or an assembly graph with sequences 
in gfa format, 2) circle rescue fasta file, and 3) traditional binner bins, 
either a directory containing the bins (fasta/fastq/fna/...), or 
a list of files (via bash expansion).
All files are read once.
See also `-h`. 

The output is a single file or lines written to STDOUT, containing 
circle resuce and bin names. Together, they form the merged binning result.

Example of a run from scratch, as done in the manuscript:

```
# assembly
hifiasm-meta -t32 -o asm hifi.fq.gz 2>log_asm

# binning plus post-processing of bins, if necessary.
# Here we assume binner_bin will be the output that 
# contain fasta files (suffix being ".fa"), each of which is a putative MAG. 

# merge circle rescue and binner bins
python bin_merger.py  -o result -x fa \
<(gfa2l.py asm.p_ctg.gfa) asm.rescue.fa binner_bin/
```

#### plotkmerspectrum.py: demo of k-mer spectrum plots

See [yam](https://github.com/xfengnefx/yam) for generating input file,
which is a 2D integer matrix of size `(1024, 1024)`. Value `x` at row `i` and column `j` 
means that there are `x` k-mers from the reads which appears exactly `x` times in reads 
and `y` times in the assembly.

`example/` provides a pair of plotting input and output. 

