import argparse
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser(description='Plot right-accumulated raio plot of a k-mer profile. Demo takes one file; see other functions for comparison between files (hatched regions; up to three files).')
parser.add_argument('f_in', type=str, help='The k-mer profile.')
parser.add_argument('f_out', type=str, help='Output file name.')
args = parser.parse_args()



def kat_type2(f, ax,
              write_prefix='', show=True,
              fzero='', separator=',', ignore_lastline=True,
              title='',
              no_1xabove=False  # do not plot the "colorful" bands
             ):
    # accumulated ratio plot
    with open(f) as file:
        d = []
        for line in file:
            if line[0]=='#': continue
            line = line.strip().split(separator)
            tmp = [int(_) for _ in line if _!='']
            tmp = [0 if _<0 else _ for _ in tmp ]  # monkey patch at first column - need to fix the overflow
            d.append(tmp)

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
        if not no_1xabove or (no_1xabove and i<=1):
            ax.plot(d2[i]+bottom, alpha=alpha)
            ax.fill_between(np.arange(L), d2[i]+bottom,bottom, alpha=alpha)
        bottom += d2[i]

    tmp = np.sum(d2[LIM:], axis=0)
    if not no_1xabove:
        ax.plot(tmp+bottom, alpha=alpha, color='gray')
        ax.fill_between(np.arange(L), tmp+bottom, bottom,
                        alpha=alpha, color='gray')


    # absolute count right accumulated
    a = np.sum(d, axis=1)
    a = a[::-1]
    aa = [a[0]]
    for i in range(1, len(a)):
        aa.append(aa[i-1]+a[i])
    aa = np.array(aa[::-1])
    for i_whiteout in range(len(aa)):
        if aa[i_whiteout]<1e6:break

    ##### absolute count right-accumulated curve
    ax2 = ax.twinx()
    ax2.plot(aa, color='white', linewidth=1, linestyle='dashed')
    ax2.set_yscale('log')


#     ##### soft white out
#     ax.fill_between(range(i_whiteout, len(aa)),
#                     [0]*(len(aa)-i_whiteout),
#                     [1]*(len(aa)-i_whiteout),
#                     color='white',
#                     alpha=0.8
#                    )
    ##### alternative: hard white out
    ax.fill_between(range(i_whiteout, len(aa)),
                    [0]*(len(aa)-i_whiteout),
                    [1]*(len(aa)-i_whiteout),
                    color='white',
                    alpha=1, zorder=2
                   )
#     #### alternative, put a red vertical bar
#     for i_whiteout in range(len(aa)):
#         if aa[i_whiteout]<1e6:break
#     ax.axvline(i_whiteout, color='red', alpha=0.5, linestyle='dashed')


    if title=='':
        ax.set_title(f+' | type2')
    else:
        ax.set_title(title)
    return ret, aa


def kat_type3(f_prev, f_new, ax,
              write_prefix='', show=True,
              fzero='', separator=',', ignore_lastline=True,
              title='',
             no_1xabove=False):
    data = []
    for f in [f_prev, f_new]:
        # accumulated ratio plot
        with open(f) as file:
            d = []
            for line in file:
                if line[0]=='#': continue
                line = line.strip().split(separator)
                tmp = [int(_) for _ in line if _!='']
                tmp = [0 if _<0 else _ for _ in tmp ]  # monkey patch - need to fix the overflow
                d.append(tmp)

        d = np.array(d)

        d2 = []
        for i in range(d.shape[0]-1):
            tmp = d[i:, :].sum(axis=0)
            tmp = tmp/(np.sum(tmp)+1)
            d2.append(tmp)
        d2.append(d[-1, :]/(np.sum(d[-1, :])+1))
        d2 = np.array(d2)
        d2 = d2.T
        data.append(d2)


    LIM = 15
    # get the watermark from old run
    watermark = data[0][0]
    # plot new run
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
        if not no_1xabove or (no_1xabove and i<=1):
            ax.plot(d2[i]+bottom, alpha=alpha)
            ax.fill_between(np.arange(L), d2[i]+bottom,bottom, alpha=alpha,
                           linewidth=0.2)
        bottom += d2[i]

    tmp = np.sum(d2[LIM:], axis=0)
    if not no_1xabove:
        ax.plot(tmp+bottom, alpha=alpha, color='gray')
        ax.fill_between(np.arange(L), tmp+bottom, bottom, linewidth=0.2,
                        alpha=alpha, color='gray')
    # plot watermark
    ax.fill_between(np.arange(L), data[1][0], data[0][0], hatch='|',
                          linewidth=0.2,facecolor="#ffa154", edgecolor='black')
#     ax_water.plot(watermark, color='white', linewidth=1, alpha=0.7)

    # absolute count right accumulated
    a = np.sum(d, axis=1)
    a = a[::-1]
    aa = [a[0]]
    for i in range(1, len(a)):
        aa.append(aa[i-1]+a[i])
    aa = np.array(aa[::-1])
    for i_whiteout in range(len(aa)):
        if aa[i_whiteout]<1e6:break

    ##### absolute count right-accumulated curve
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



#     ##### soft white out
#     ax.fill_between(range(i_whiteout, len(aa)),
#                     [0]*(len(aa)-i_whiteout),
#                     [1]*(len(aa)-i_whiteout),
#                     color='white',
#                     alpha=0.8
#                    )
    ##### alternative: hard white out
    ax.fill_between(range(i_whiteout, len(aa)),
                    [0]*(len(aa)-i_whiteout),
                    [1]*(len(aa)-i_whiteout),
                    color='white',
                    alpha=1, zorder=2
                   )


    if title=='':
        ax.set_title(f+' | type3')
    else:
        ax.set_title(title)
    return ret, aa


def kat_type4(f_prev, f_mid, f_now, ax,
              write_prefix='', show=True,
              fzero='', separator=',', ignore_lastline=True,
              title='',
             no_1xabove = False):
    data = []
    for f in [f_prev, f_mid, f_now]:
        # accumulated ratio plot
        with open(f) as file:
            d = []
            for line in file:
                if line[0]=='#': continue
                line = line.strip().split(separator)
                tmp = [int(_) for _ in line if _!='']
                tmp = [0 if _<0 else _ for _ in tmp ]  # monkey patch - need to fix the overflow
                d.append(tmp)
        d = np.array(d)
        d2 = []
        for i in range(d.shape[0]-1):
            tmp = d[i:, :].sum(axis=0)
            tmp = tmp/(np.sum(tmp)+1)
            d2.append(tmp)
        d2.append(d[-1, :]/(np.sum(d[-1, :])+1))
        d2 = np.array(d2)
        d2 = d2.T
        data.append(d2)


    LIM = 15
    # plot new run
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
        if not no_1xabove or (no_1xabove and i<=1):
            ax.plot(d2[i]+bottom, alpha=alpha)
            ax.fill_between(np.arange(L), d2[i]+bottom,bottom, alpha=alpha,
                           linewidth=0.2)
        bottom += d2[i]

    tmp = np.sum(d2[LIM:], axis=0)
    if not no_1xabove:
        ax.plot(tmp+bottom, alpha=alpha, color='gray')
        ax.fill_between(np.arange(L), tmp+bottom, bottom, linewidth=0.2,
                        alpha=alpha, color='gray')
    # plot watermark
    ax.fill_between(np.arange(L), data[1][0], data[0][0], hatch='/',
                          linewidth=0.2,facecolor="#ffa154", edgecolor='black')
    ax.fill_between(np.arange(L), data[2][0], data[1][0], hatch='\\',
                          linewidth=0.2,facecolor="#ffc653", edgecolor='black')
#     ax_water.plot(watermark, color='white', linewidth=1, alpha=0.7)

    # absolute count right accumulated
    a = np.sum(d, axis=1)
    a = a[::-1]
    aa = [a[0]]
    for i in range(1, len(a)):
        aa.append(aa[i-1]+a[i])
    aa = np.array(aa[::-1])
    for i_whiteout in range(len(aa)):
        if aa[i_whiteout]<1e6:break

    ##### absolute count right-accumulated curve
    ax2 = ax.twinx()
    ax2.plot(aa, color='white', linewidth=1, linestyle='dashed')
    ax2.set_yscale('log')


#     ##### soft white out
#     ax.fill_between(range(i_whiteout, len(aa)),
#                     [0]*(len(aa)-i_whiteout),
#                     [1]*(len(aa)-i_whiteout),
#                     color='white',
#                     alpha=0.8
#                    )
    ##### alternative: hard white out
    ax.fill_between(range(i_whiteout, len(aa)),
                    [0]*(len(aa)-i_whiteout),
                    [1]*(len(aa)-i_whiteout),
                    color='white',
                    alpha=1, zorder=2
                   )


    if title=='':
        ax.set_title(f+' | type3')
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
