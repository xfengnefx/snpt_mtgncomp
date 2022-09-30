import sys, argparse, os, stat, statistics
from fileutil import opener
from readfq import readfq
from time import time

example_text = '''example:
   python {0} <(gfa2l.py asm.gfa) asm.rescue.fa /path/metabat2/
   python {0} contiglens.tsv asm.rescue.fa.gz /path/metabat2
   python {0} contiglens.tsv <(zcat asm*.fa.gz) /path/metabat2'''.format(os.path.basename(__file__))

parser = argparse.ArgumentParser(description="Merge rescued circles and genome binner bins, write to stdout or specified file.",
                                epilog=example_text,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-l', type=str, default='', help='a file conting a list of contigs, will report any of them if not used by rescueing *and* binner.')
parser.add_argument('-x', type=str, default='fa', help='file suffix of binner bins. Ignored if files instead of a directory are specified. ')
parser.add_argument('-o', type=str, default='', help='write output to this file instead of stdout')
parser.add_argument('fn_contigl_tsv', type=str, help='a tab-delimited file. col1=contig names, col2=contig lengths')
parser.add_argument('fn_rescue_fasta', type=str, help='fasta file of rescued circles')
parser.add_argument('binner', nargs='+', type=str, help='folder contains binner bins or fasta/q files, can be mixed. Will assume plain text if file is a pipe.')
args = parser.parse_args()

accept_list = {}  # [type, size, n_contigs, contignames]
blacklist = set()
bin2tig = {}
maptignames = {}
reject = []
ll = {}
long_linear_contigs = set()
threshold = 10

# collect contig lengths
# also record truncated names - TODO update hifiasm-meta
T = time()
with open(args.fn_contigl_tsv) as file:
    for line in file:
        qn, l = line.strip().split('\t')
        l = int(l)
        qnshort = qn.split('.')[1]
        maptignames[qnshort] = qn
        ll[qn] = l
        
# anything to preserve?
if args.l!='':
    file = opener(args.l)
    for line in file:
        long_linear_contigs.add(line.strip())
sys.stderr.write('[M] prepare done. %.1fs\n'%(time()-T))

# record rescued circles
T = time()
with open(args.fn_rescue_fasta) as file:
    for name, seq, qual in readfq(file, get_full_name=True):
        name = name.strip().split(' ')
        pfname = name[0]
        tigs = [maptignames[_[:-1]+'l'] for _ in name[1:]]
        bin2tig[pfname] = tigs
        tigs_ll = [ll[tig] for tig in tigs]
        if (len(tigs)>=threshold and statistics.median(tigs_ll)<500000) or\
           len([_ for _ in tigs_ll if _>100000])>2:
            accept_list[pfname] = ['rc', len(seq), 1, pfname]
            for tig in tigs:
                blacklist.add(tig)
sys.stderr.write(f'[M] accept rescue circles: {len(accept_list)} (%.1fs)\n'%(time()-T))

# record binner bins
T = time()
def yield_binner_streams(fns, suffix='fa'):
    n = 0
    nskip = 0
    for fn in fns:
        if os.path.isfile(fn):
            yield fn, opener(fn)
            n+=1
        elif os.path.isdir(fn):
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
    for line in stream:
        for qname, seq, qual in readfq(stream):
            tigs.append(qname)
        if len(tigs)>0:
            bin2tig[b] = tigs
            
    binsize = sum([ll[tig] for tig in tigs])
    if binsize<500000 or binsize>10000000:
        reject.append(b)
        continue

    used = [_ for _ in tigs if _ in blacklist]
    if len(used)/len(tigs)>0.1 and len(used)>5:
        reject.append(b)
    else:
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
    
       
