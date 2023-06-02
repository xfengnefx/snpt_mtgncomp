import os, sys, argparse, random
from time import time

runID = str(int(random.random()*1e6))
def add_forwardslash(d):
    if not d.endswith('/'):return d+'/'
    return d

parser = argparse.ArgumentParser(description="run metabate(including depth estimatin via mapping)")
parser.add_argument('--dry_run', action='store_true', default=False, help='dry run')
parser.add_argument('-t', type=int, default=48, help='number of threads')
parser.add_argument('-n', type=str, default='placeholder', help='prefix of output files')
parser.add_argument('-f', action='store_true', default=True, help='remove contents in the target folder if there is any')
parser.add_argument('--min_length', type=int, action='store', default=-1, help='-s option for metabat2; default is -1 (i.e. use mt2 default which should be 2000)')
parser.add_argument('--tmp', type=str, default='/hlilab/xfeng/prj_meta/tmp', help='tmp folder')
parser.add_argument('--list_contigs', type=str, default='', help="give a list of contig files (can NOT mix gz and non-gz); positional arg will be treated as a dummy")
parser.add_argument('--list_reads', type=str, default='', help='give a list of read files (can NOT mix gz and non-gz); positional arg will be treated as a dummy')
parser.add_argument('f_ctg', type=str, default='', help='contig fa/fq(.gz)')
parser.add_argument('f_read', nargs='+',type=str, default='', help='read fa/fq(.gz)')
parser.add_argument('workdir', type=str, default='./binned_raw', help='output directory')

args = parser.parse_args()
is_dry = args.dry_run
threads = args.t
f_ctg = args.f_ctg
fs_read = args.f_read
forced = args.f
dir_tmp = add_forwardslash(args.tmp)
dir_out = add_forwardslash(args.workdir)

mt2_contig_length_option = ''
if args.min_length>0:
    mt2_contig_length_option = '-s '+str(args.min_length)




if not os.path.exists(dir_tmp):
    cmd = 'mkdir -p '+dir_tmp
    if is_dry: print(cmd)
    else: os.system(cmd)


if args.n=='placeholder':
    prefix = 'mb2_'+f_ctg.split('/')[-1]
else:
    prefix = args.n
if not os.path.exists(dir_out):
    cmd = 'mkdir -p '+dir_out
    if is_dry: print(cmd)
    else: os.system(cmd)
else:
    assert dir_out!='/'
    san = os.listdir(dir_out)
    if len(san)>0:
        if not forced:
            print('output folder exists and not empty, specify -f to force.')
            exit(1)
        else:
            cmd= 'rm -rf {0}*'.format(dir_out)
            if is_dry: print(cmd)
            else: os.system(cmd)

f_ctglist = args.list_contigs
f_readlist = args.list_reads
fs_ctg = []
f_tmp_ctg = 'mt2tmp_{1}.{0}.contig.fa'.format(args.n, runID)  # could be gz file; just suffix doesn't matter
f_tmp_read = 'mt2tmp_{1}.{0}.read.fa'.format(args.n, runID)  # could be gz file; just suffix doesn't matter.
if f_ctglist!='':
    with open(f_ctglist) as file:
        for line in file:fs_ctg.append(os.path.abspath(line.strip()))
    cmd = 'cat {0} > {1}'.format(' '.join(fs_ctg), os.path.join(dir_tmp, f_tmp_ctg))
    if is_dry: print(cmd)
    else: os.system(cmd)
    f_ctg = os.path.join(os.path.abspath(dir_tmp), f_tmp_ctg)
if f_readlist!='':
    fs_read = []
    with open(f_readlist) as file:
        for line in file:fs_read.append(os.path.abspath(line.strip()))
    cmd = 'cat {0} > {1}'.format(' '.join(fs_read), os.path.join(dir_tmp, f_tmp_read))
    if is_dry: print(cmd)
    else: os.system(cmd)


bamname = dir_tmp+prefix+'.aln.bam'
readfiles = ' '.join(fs_read)
if os.path.isfile(bamname):
    print('BAM file exists, skip alignment')
    cmd_aln = []
else:
    print('will align readsets ({0}): '.format(len(fs_read)), readfiles)
    cmd_aln = ['minimap2 -a -k 19 -w 10 -I 10G -g 5000 -r 2000 '+\
            '--lj-min-ratio 0.5 -A 2 -B 5 -O 5,56 -E 4,1 -z 400,50 ' +\
            f'--sam-hit-only -t {threads} {f_ctg} {readfiles} '+\
            f'2> {dir_out}{prefix}.log_mm2 | '+\
            f'samtools sort -m4g -@8 -o {bamname}']
depthname = f'{dir_out}{prefix}.depth'
if (os.path.isfile(depthname)):
    print('depth file exists, skip jgi')
    cmd_jgi = []
else:
    cmd_jgi = ['jgi_summarize_bam_contig_depths --outputDepth '+\
        f'{dir_out}{prefix}.depth {bamname} '+\
        f'2> {dir_out}{prefix}.log_depth']


cmds = ['metabat2 --seed 1 {6} -i {0} -a {1}{2} -o {1}{3} -t {4} -v &> {1}{5}'.format(f_ctg, dir_out, prefix+'.depth',
                                                                     prefix+'.bin',
                                                                   threads,
                                                                   prefix+'.log_metabat2', mt2_contig_length_option),
        f'touch {dir_out}endmark'
]
cmds = cmd_aln + cmd_jgi + cmds

start_time = time()

for cmd in cmds:
    if is_dry:
        print(cmd)
        continue
    print('executing: ', cmd)
    ret = os.system(cmd)
    if ret!=0:
        if len(fs_ctg)>0:
            os.system('rm {0}'.format(f_tmp_ctg))
        if len(fs_read)>0:
            os.system('rm {0}'.format(f_tmp_read))
        print('failed on the followingn command, check log files:')
        print(cmd)
        exit(1)
#cmd = ';'.join(cmds)
#os.system(cmd)
with open(dir_out+prefix+'.time', 'w') as file:
    file.write('{0} s\n'.format(time()-start_time))
