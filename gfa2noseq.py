import argparse
from fileutil import opener

parser = argparse.ArgumentParser(description='Strip sequences from a gfa files.')
parser.add_argument('f_in', type=str, help='input file name')
parser.add_argument('f_out', type=str, help='output file name')
args = parser.parse_args()

def remove_seq(f_in, f_out):
    with open(f_out, 'w') as file_out:
        file_in = opener(f_in)
        for line in file_in:
            if line[0]!='S':
                file_out.write(line)
            else:
                has_len = 'LN:i' in line
                line = line.strip().split('\t')
                l = len(line[2])
                if not has_len:
                    line = line[:2]+['*', f'LN:i:{l}']+line[3:]
                else:
                    line = line[:2]+['*']+line[3:]
                line = '\t'.join(line)+'\n'
                file_out.write(line)

if __name__=='__main__':
    remove_seq(args.f_in, args.f_out)
