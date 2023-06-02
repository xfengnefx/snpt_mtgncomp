import argparse
example_text = '''example:
   python {0} <(gfa2l.py asm.gfa) asm.rescue.fa /path/metabat2/
   python {0} contiglens.tsv asm.rescue.fa.gz /path/metabat2
   python {0} contiglens.tsv <(zcat asm*.fa.gz) /path/metabat2'''.format(os.path.basename(__file__))

parser = argparse.ArgumentParser(description="Merge rescued circles and genome binner bins"+\
" from traditional binner. Write a tab-delimited text file to stdout or specified file.",
                 epilog=example_text,
                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-x', type=str, default='fa', help='file suffix of binner bins, e.g. use fa for *.fa files. Ignored if files instead of a directory are specified.')
parser.add_argument('-o', type=str, default='', help='write output to this file instead of stdout')
parser.add_argument('fn_contigl_tsv', type=str, help='a tab-delimited file. col1=contig names, col2=contig lengths.')
parser.add_argument('fn_rescue_fasta', type=str, help='fasta file of rescued circles.')
parser.add_argument('binner', nargs='+', type=str, help='folder contains binner bins or fasta/q files, can be mixed. Will assume plain text if file is a pipe.')
args = parser.parse_args()

import sys, argparse, os, stat, statistics
from fileutil import opener
from readfq import readfq
from time import time



accept_list = {}     # in each entry: [type, size, n_contigs, contignames]
blacklist = set()    # remembers what contigs have been used
bin2tig = {}         # key=bin name, value=a list of contig names belong to that bin
reject = []          # (for sancheck)
ll = {}              # tig2len
long_linear_contigs = set()  # (for sancheck)



# collect contig lengths
T = time()
with open(args.fn_contigl_tsv) as file:
    for line in file:
        qn, l = line.strip().split('\t')
        l = int(l)
        qnshort = qn.split('.')[1]
        ll[qn] = l
sys.stderr.write('[M] prepare done. %.1fs\n'%(time()-T))

# collect rescued circles:
#   all are accepted; also mark all of their contigs as used
T = time()
with open(args.fn_rescue_fasta) as file:
    for name, seq, qual in readfq(file, get_full_name=True):
        name = name.strip().split(' ')
        pfname = name[0]
        tigs = [_[:-1] for _ in name[1:]]
        bin2tig[pfname] = tigs
        tigs_ll = [ll[tig] for tig in tigs]
        accept_list[pfname] = ['rc', len(seq), 1, pfname]
        for tig in tigs:
            blacklist.add(tig)
sys.stderr.write(f'[M] accept rescue circles: {len(accept_list)} (%.1fs)\n'%(time()-T))

# parse binner directory or a list of filesfiles
T = time()
def yield_binner_streams(fns, suffix='fa'):
    n = 0
    nskip = 0
    for fn in fns:
        if os.path.isfile(fn):
            yield fn, opener(fn)
            n+=1
        elif os.path.isdir(fn):  # need to iterate over files in the directory
            for f in os.listdir(fn):
                if f.endswith(suffix):
                    ff = os.path.join(fn,f)
                    yield ff, opener(ff)
                    n+=1
        elif stat.S_ISFIFO(os.stat(fn).st_mode) or fn.startswith('/dev/fd/'): # assume pipe
            fh = open(fn)
            yield fn, fh
            fh.close()
            n+=1
        else:  # unknown
            sys.stderr.write(f'[W] {fn} is not a file or a directory? will continue\n')
            nskip+=1
    sys.stderr.write(f'[M] processed {n} files ({nskip} ignored) from binner\n')
for b, stream in yield_binner_streams(args.binner):
    tigs = []

    # collect names of the bin's contigs
    for qname, seq, qual in readfq(stream):
        tigs.append(qname)

    # skip empty bin, if any
    if len(tigs)>0:
        bin2tig[b] = tigs
    else:
        continue

    # ignore binner bins that are too small or huge
    # this is for reducing checkM runtime.
    binsize = sum([ll[tig] for tig in tigs])
    if binsize<500000 or binsize>10000000:
        reject.append(b)
        continue

    # check redundancy
    used = [_ for _ in tigs if _ in blacklist]
    used_bp = sum([ll[tig] for tig in used])
    if used_bp>1e6 or len(used)>10:  # conflicts with circle rescue. ignore binner bin
        reject.append(b)
    else:                            # take the binner bin
        for tig in tigs:
            blacklist.add(tig)
        accept_list[b] = ['binner', binsize, len(tigs), ','.join(tigs)]
sys.stderr.write('[M] accepted total: {0}; rejected from binner: {1}. ({2})\n'.format(
                len(accept_list),
                len(reject),
                '%.1fs'%(time()-T)))

# output
if args.o!='':
    file_out = open(args.o, 'w')
else:
    file_out = sys.stderr
file_out.write('\t'.join(['#binname', 'type', 'size', 'n_contigs', 'contigs'])+'\n')
for b in accept_list:
    tmp = [b]+[str(_) for _ in accept_list[b]]
    file_out.write('\t'.join(tmp)+'\n')
if args.o!='':
    file_out.close()
