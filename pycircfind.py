import argparse
from time import time
from fileutil import opener

parser= argparse.ArgumentParser(description='A freestanding circle rescue script.')
parser.add_argument('f_in', type=str, help='assembly graph file name') 
parser.add_argument('f_out', type=str, help='output file name')  
args = parser.parse_args() 

trans = str.maketrans('ATCGatcg', 'TAGCtagc')
def revcmp(s):
    return s.translate(trans)[::-1]
def toggle_strand(c):
    if c=='+': return '-'
    elif c=='-': return '+'
    else:
        print('unknown strand:', c)
        raise ValueError
    
class stack:
    def __init__(self):
        self.a = []
    def push(self, d):
        self.a.append(d)
    def pop(self):
        x = self.a[-1]
        self.a = self.a[:-1]
        return x
    def is_empty(self):
        if len(self.a)==0:
            return True
        else:
            return False
    
class node:
    def __init__(self, ID:str, strand:str):
        self.name = ID+strand
        self.strand = strand
        self.arc = {}
    def add_arc(self, targetID:str, targetStrand:str, l:int):
        targetName = targetID+targetStrand
        if targetName in self.arc:
            print('[W] trying to add duplicated arc, ignore:', self.ID, targetName)
        else:
            self.arc[targetName] = l
    def has_child(self):
        if len(self.arc)>0: 
            return True
        else:
            return False
class graph:
    def __init__(self):
        self.nodes = {}
        self.seqs = {}
        self.size = 0
    def build_from_gfa(self, f_in:str):
        file = opener(f_in)
        for line in file:
            line = line.strip().split('\t')
            if line[0]=='S':
                if line[2]=='*':
                    print('[E] must use gfa with seqs')
                    raise ValueError
                ID = line[1]
                seq = line[2]
                self.add_node(ID, seq)
            elif line[0]=='L':
                _, ID1, strand1, ID2, strand2, l = line[:6]
                l = int(l[:-1])  # assuming it is [0-9]+M
                self.add_arc(ID1, strand1, ID2, strand2, l)
            else:
                pass
    def add_node(self, ID:str, seq:str):
        if ID in self.nodes:
            print('[W] trying to add duplicated node, ignore:', ID, len(seq))
        else:
            self.nodes[ID+'+'] = node(ID, '+')
            self.nodes[ID+'-'] = node(ID, '-')
            self.seqs[ID] = seq
            self.size+=1
    def add_arc(self, ID1:str, strand1:str, 
                    ID2:str, strand2:str, 
                   len_arc:int):
        self.nodes[ID1+strand1].add_arc(ID2, strand2, len_arc)
        self.nodes[ID2+toggle_strand(strand2)].add_arc(ID1, toggle_strand(strand1), len_arc)
    
    def get_targets(self, ID, strand):
        for name2 in self.nodes[ID+strand].arc:
            yield self.nodes[name2]
    def get_target_IDs_with_strands(self, ID, strand):
        for name2 in self.nodes[ID+strand].arc:
            yield name2[:-1], name2[-1]
    def get_targetNames(self, name):
        for name2 in self.nodes[name].arc:
            yield name2
    def get_node(self, name):
        return self.nodes[name]
    def get_node_len(self, name):
        return len(self.seqs[name[:-1]])
    def get_node_seq(self, name):
        strand = name[-1]
        if strand=='+':
            return self.seqs[name[:-1]]
        else:
            return revcmp(self.seqs[name[:-1]])
    def get_path_seq(self, namelist, is_circ):
        if len(namelist)==0: return ''
        name = namelist[0]
        seq = self.get_node_seq(name)
        for i, name in enumerate(namelist[1:]):
            arcl = self.get_node(namelist[i]).arc[name]
            seqtmp = self.get_node_seq(name)
            seq+=seqtmp[arcl:]
        if is_circ:
            arcl = self.get_node(name).arc[namelist[0]]
            if arcl>=len(seq):
                print('[W] circular, but ovlp is longer than seq. Do nothing.')
            else:
                seq = seq[:-arcl]
        return seq
    
    def do_circle_finding_core_DFS(self, namebase, min_len, max_len, max_used, used:set):
        verbose = False
        if not self.get_node(namebase).has_child():
            return 0, []  
        s = stack()
        color = {}
        tot_len = 0
        tot_used = 0
        tot_usedlen = 0
        
        #init
        s.push(namebase)
        color[namebase] = 1  # gray
        if verbose: print(f'[d] push: {namebase}')
        
        # go
        report = False
        while (not s.is_empty()) and (tot_len<max_len) and (tot_used<max_used):
            name = s.pop()
            if verbose: print(f'[d] pop: {name}')
            if color[name]==2: continue # black, node has been finished
                
            # special case: found the circle
            if namebase in self.get_targetNames(name) and\
                tot_len>=min_len:
                report = True
                s.push(name)  # put the last one back!
                break
                
            # step
            tmp = []
            for name2 in self.get_targetNames(name):
                tmp.append([int(name2 in used), name2])
            tmp.sort(reverse=True)
            if verbose: print(f'[d]     buffer len: {len(tmp)}')
            for _, name2 in tmp:   # step into unused nodes first
                if verbose: print(f'[d]      looking at: {name2}')
                if name2 not in color: # white node
                    s.push(name)
                    s.push(name2)
                    if verbose: print(f'[d]   put back {name}, push {name2}')
                    color[name2] = 1
                    if name2 in used: 
                        tot_used +=1
                        tot_usedlen += self.get_node_len(name2)
                    tot_len+=self.get_node_len(name2)
                    break
            else: #no white node
                if verbose: print(f'[d] drop: {name}')
                color[name] = 2
                tot_len-=self.get_node_len(name)
                if name in used: tot_used-=1
        
        if tot_usedlen/tot_len>0.7:
            report = False
        if report: 
            return 1, s.a
        return 0, []
                    
    
    def do_circle_finding_core(self, min_len, max_len, max_used):
        used = set()
        for namebase in self.nodes:
            if namebase[-1]=='-': continue  # only need to traverse one direction
            
            found, names = self.do_circle_finding_core_DFS(namebase, 
                                                       min_len, max_len, 
                                                       max_used, used)
            if found:
                for name in names:
                    used.add(name)
                    used.add(name[:-1]+toggle_strand(name[-1]))
                yield names
    def do_circle_finding_debug(self, namebase, min_len=1e6, max_len=1e7, max_used = 200):
        used = set()
            
        found, names = self.do_circle_finding_core_DFS(namebase, 
                                                   min_len, max_len, 
                                                   max_used, used)
        return found, names
    def do_circle_finding(self, f_out,min_len=1e6, max_len=1e7, max_used = 200):
        with open(f_out, 'w') as file:
            i = 0
            for names in self.do_circle_finding_core(min_len, max_len, max_used):
                tmp = ','.join(names)
                seq = self.get_path_seq(names, True)
                file.write(f'>rc.{i} {len(seq)} {tmp}\n{seq}\n')
                i+=1
        return i

if __name__=='__main__':
    T = time()
    g = graph()
    g.build_from_gfa(args.f_in)
    n = g.do_circle_finding(args.f_out)
    Td = '%.1fs'%(time()-T)
    print(f'[M] reported {n} circles from {args.f_in}. (Td)')
