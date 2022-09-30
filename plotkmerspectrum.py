import argparse
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser(description='Plot right-accumulated raio plot of a k-mer profile.')
parser.add_argument('f_in', type=str, help='The k-mer profile.')
parser.add_argument('f_out', type=str, help='Output file name.')
args = parser.parse_args()


def kat_typetall(f, ax,  fzero='', separator=',', ignore_lastline=True):
    # a very tall plot showing absolute counts 
    with open(f) as file:
        d = []
        for line in file:
            if line[0]=='#': continue
            line = line.strip().split(separator)
            d.append([int(_) for _ in line if _!=''])
    d = np.array(d)
    d = d.T
    print(d.shape)

    LIM = 15
    bottom = np.zeros(len(d[0]))
    L = len(bottom)
    for i, _ in enumerate(d[:LIM]):
        ax.bar(np.arange(L), d[i],bottom=bottom)
        bottom += d[i]
    #     plt.bar(np.arange(L), np.log10(d[i]+1),bottom=bottom)
    #     bottom += np.log10(d[i]+1)

    tmp = np.sum(d[LIM:], axis=0)
    ax.bar(np.arange(L), tmp, bottom = bottom, label='+')
    # plt.bar(np.arange(L), np.log10(tmp+1), bottom = bottom, label='+')

    ##### zero dump
    if fzero!='':
        with open(fzero, 'rb') as file:
            zerotot, zero_n_one_mismatch, zero_data, zero_data_one_mismatch = pickle.load(file)
        d_zero = []
        for i in range(d.shape[0]):
            if i in zero_data: d_zero.append(zero_data[i])
            else: d_zero.append(0)
        d_zero = np.array(d_zero)
        ax.bar(np.arange(L), d_zero, bottom =d[0]-d_zero, color='gray')
    ###########

    ax.set_ylim(bottom=0, top=1000000)


def kat_type1(f, ax,
              fzero='', separator=',', ignore_lastline=True):
    # plain ratio plot
    with open(f) as file:
        d = []
        for line in file:
            if line[0]=='#': continue
            line = line.strip().split(separator)
            d.append([int(_) for _ in line if _!=''])
    if ignore_lastline: d = np.array(d[:-1])
    else: d = np.array(d)

    d2 = []
    for _ in d:
        tmp = _/np.sum(_)
        d2.append(tmp)
    d2 = np.array(d2)
    d2 = d2.T

    LIM = 15
    bottom = np.zeros(len(d2[0]))
    L = len(bottom)
    for i, _ in enumerate(d2[:LIM]):
        ax.plot(d2[i]+bottom)
        ax.fill_between(np.arange(L), d2[i]+bottom,bottom)
        bottom += d2[i]

    tmp = np.sum(d2[LIM:], axis=0)
    ax.plot(tmp+bottom)
    ax.fill_between(np.arange(L), tmp+bottom, bottom)

    ##### zero dump
    if fzero!='':
        with open(fzero, 'rb') as file:
            zerotot, zero_n_one_mismatch, zero_data, zero_data_one_mismatch = pickle.load(file)
        d_zero = []
        for i in range(d.shape[0]):
            if i in zero_data: 
                d_zero.append(zero_data[i]/sum(d[i]))
            else: d_zero.append(0)
        d_zero = np.array(d_zero)
        ax.plot(d2[0], color='gray')
        ax.fill_between(np.arange(L), d2[0], d2[0]-d_zero, color='gray')
    ###########
    ax.set_title(f+' | type1')

def kat_type2(f, ax,
              write_prefix='', show=True,
              fzero='', separator=',', ignore_lastline=True,
              title=''):
    # right-accumulated ratio plot
    with open(f) as file:
        d = []
        for line in file:
            if line[0]=='#': continue
            line = line.strip().split(separator)
            d.append([int(_) for _ in line if _!=''])
    d = np.array(d)

    d2 = []
    for i in range(d.shape[0]-1):
        tmp = d[i:, :].sum(axis=0)
        tmp = tmp/(np.sum(tmp)+1)
        d2.append(tmp)
    d2.append(d[-1, :]/(np.sum(d[-1, :])+1))
    d2 = np.array(d2)
    d2 = d2.T


    LIM = 15
    bottom = np.zeros(len(d2[0]))
    L = len(bottom)
    ret = ''
    for i, _ in enumerate(d2[:LIM]):
        if i>1:
            alpha=0.3
        else:
            alpha=1
        if i==0:
            ret = d2[i]+bottom
        ax.plot(d2[i]+bottom, alpha=alpha)
        ax.fill_between(np.arange(L), d2[i]+bottom,bottom, alpha=alpha)
        bottom += d2[i]

    tmp = np.sum(d2[LIM:], axis=0)
    ax.plot(tmp+bottom, alpha=alpha, color='gray')
    ax.fill_between(np.arange(L), tmp+bottom, bottom, 
                    alpha=alpha, color='gray')
    
    
    ##### absolute count right-accumulated curve
    a = np.sum(d, axis=1)
    a = a[::-1]
    aa = [a[0]]
    for i in range(1, len(a)):
        aa.append(aa[i-1]+a[i])
    aa = np.array(aa[::-1])
    ax2 = ax.twinx()
    ax2.plot(aa, color='white', linewidth=1, linestyle='dashed')
    ax2.set_yscale('log')
    

    ##### zero dump
    if fzero!='':
        with open(fzero, 'rb') as file:
            zerotot, zero_n_one_mismatch, zero_data, zero_data_one_mismatch = pickle.load(file)
        d_zero = []
        for i in range(d.shape[0]):
            if i in zero_data: 
                d_zero.append(zero_data[i]/sum(d[i:, 0]))
            else: d_zero.append(0)
        d_zero = np.array(d_zero)
        ax.plot(d2[0], color='gray')
        ax.fill_between(np.arange(L), d2[0], d2[0]-d_zero, color='gray')
    ###########
    
#     ##### white out
#     for i in range(len(aa)):
#         if aa[i]<1e6:break
#     ax.fill_between(range(i, len(aa)), 
#                     [0]*(len(aa)-i),
#                     [1]*(len(aa)-i),
#                     color='white', 
#                     alpha=0.8
#                    )
    #### alternative, put a red vertical bar
    for i in range(len(aa)):
        if aa[i]<1e6:break
    ax.axvline(i, color='red', alpha=0.5, linestyle='dashed')
    
    
    if title=='':
        ax.set_title(f+' | type2')
    else:
        ax.set_title(title)
    return ret, aa

if __name__=='__main__':
    fig, axes = plt.subplots(1, 1)
    fig.set_size_inches(16, 4)
    ax = plt.gca()
    kat_type2(args.f_in, ax, separator=' ', ignore_lastline=True,
                title=f'{args.f_in}')
    ax.set_xlim(0, 1023)
    ax.set_ylim(0, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    plt.savefig(args.f_out, bbox_inches='tight', format='pdf')
    plt.close()
