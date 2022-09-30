from collections import namedtuple


def _parse_tags(l: list) -> dict:
    # parse SAM/PAF/GFA-like list of tags
    ret = {}
    for tag in l:
        n, t, v = tag.split(':')
        if t=='i': ret[n] = int(v)
        elif t=='f': ret[n] = float(v)
        elif t=='Z' or t=='A': ret[n] = v
        elif t=='H': ret[n] = v  # TODO: parse byte array
        elif t=='B':
            tmp = v.split(',')
            if tmp[0]=='f': ret[n] = [float(_) for _ in tmp[1:]]
            else: ret[n] = [int(_) for _ in tmp[1:]]
    return ret

def readgfa(file):
    attrs = 't q_n q_seq q_hasseq tags read_n read_s '\
            'read_e read_s_ontig read_strand q2_n strand1 '\
            'strand2 l_cigar'.split(' ')
    attrs_len = len(attrs)
    container = namedtuple('gfaline', ' '.join(attrs))
    for line in file:
        d = {key:None for key in attrs}
        line = line.strip().split('\t')
        d['t'] = line[0]
        if d['t']=='S':
            d['q_n'] = line[1]
            if line[2]=='*':
                d['q_hasseq']=False
                d['q_seq'] = ''
            else:
                d['q_hasseq'] = True
                d['q_seq'] = line[2]
            d['tags'] = _parse_tags(line[3:])
        elif d['t']=='A':
            d['q_n'] = line[1]
            d['read_s_ontig'] = int(line[2])
            d['read_strand'] = line[3]
            d['read_n'] = line[4]
            d['read_s'] = int(line[5])
            d['read_e'] = int(line[6])
            d['tags'] = _parse_tags(line[7:])
        elif d['t']=='L':
            d['q_n'] = line[1]
            d['strand1'] = line[2]
            d['q2_n'] = line[3]
            d['strand2'] = line[4]
            d['l_cigar'] = line[5]
            d['tags'] = _parse_tags(line[6:])
        else:
            print('[E] unknown gfa line type: {0}'.format(line[0]))
            raise ValueError
        yield container(*[d[key] for key in attrs])
