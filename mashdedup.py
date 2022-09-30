import argparse
from fileutil import opener

parser = argparse.ArgumentParser(description='Arbitrarily drop duplicates based on mash distance.')
parser.add_argument('-i', type=float, default=99, help='ANI threshold in percentage. e.g. set to 99 to deduplicate at 99% ANI.')
parser.add_argument('f_in', type=str, help='the output of `mash dist`. Plain text.')
parser.add_argument('f_out', type=str, help='a file to write to')
args = parser.parse_args()


if __name__=='__main__':
    blacklist = set()
    accepted = set()
    with open(args.f_in) as file_in:
        for line in file_in:
            qn1, qn2, s = line.split('\t')[:3]
            s = float(s)
            if s>args.i:
                accepted.add(qn1)
                accepted.add(qn2)
            else:
                if qn1 in accepted:
                    blacklist.add(qn2)
                elif qn2 in accepted:
                    blacklist.add(qn1)
                else:
                    accepted.add(qn1)
                    blacklist.add(qn2)
    with open(args.f_out, 'w') as file_out:
        for qn in accepted: 
            file_out.write(qn+'\n')
