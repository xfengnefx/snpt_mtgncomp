import sys, argparse
from readgfa import readgfa
from fileutil import opener

parser = argparse.ArgumentParser(description="Read gfa from stdint or given files, print to stdout a tsv with fields {seqn:str, seql:int}.")
parser.add_argument('-t', default='LN', type=str, help='Sequence length tag (default: LN)')
parser.add_argument('fn_gfa', nargs='+', type=str, help='List of files. Use - for stdint. stdint must be not gzipped.')
args = parser.parse_args()
tag = args.t

n = 0
for fn in args.fn_gfa:
    if fn=='-': 
        stream = sys.stdin
    else: 
        stream = opener(fn)
    for d in readgfa(stream):
        if d.t=='S':
            if d.q_hasseq:
                sys.stdout.write(f'{d.q_n}\t{len(d.q_seq)}\n')
            else:
                try: 
                    l = d.tags[tag] 
                except: 
                    sys.stderr.write(f'[E] length of {d.q_n} not found.\n')
                sys.stdout.write(f'{d.q_n}\t{l}\n')
            n+=1
sys.stderr.write(f'[M] had {n} sequences.\n')
