import numpy as np
import matplotlib.pyplot as plt
import seaborn.apionly as sns
import matplotlib.gridspec as gridspec

import datreant as dtr
import pandas as pd

b = dtr.discover()

#tags for Treants
source = ['Lustre', 'SSD']
size = ['9000', '900']
analysis = ['RDF', 'RMS']
scheduler = ['distr', 'multi']
nodes = ['3nodes', '6nodes']
scheduler_full_name = {'distr': 'distributed', 'multi': 'multiprocessing'}


# for Lustre distributed
for ana in analysis:
    fig = plt.figure(figsize=(6, 5), tight_layout=True)
    gs = gridspec.GridSpec(5, 1)

    #900
    ax1 = fig.add_subplot(gs[0:2,0])
    ax1.set_xlabel('Number of cores')
    ax1.set_ylabel('Serial Fraction')
    ax1.set_title('pmda.{} on 900 frames'.format(ana.lower()))
    sns.despine(offset=10, ax=ax1)
    ax1.set_ylim([0,1])

    # 9000
    ax2 = fig.add_subplot(gs[2:4,0])
    ax2.set_xlabel('Number of cores')
    ax2.set_ylabel('Serial Fraction')
    ax2.set_title('pmda.{} on 9000 frames'.format(ana.lower()))
    sns.despine(offset=10, ax=ax2)
    ax2.set_ylim([0,1])
    for s in source:
        for sche in scheduler:
            if (sche == 'distr') & (s == 'Lustre'):
                for node in nodes:
                    print(node)
                    alp=0.5
                    ts = b[b.tags[['900', node, ana]]]
                    stat_csv = ts.leaves().globfilter('stat*')[0]
                    stat = pd.read_csv(stat_csv.abspath)
                    label = 'Lustre-distributed-{}'.format(node)
                    t0 = stat.total[0]
                    speedup = t0/stat.total
                    speedup_std = stat.total_std*t0/(stat.total)**2+stat.total_std[0]/stat.total
                    sp = speedup[1:]
                    sp_std = speedup_std[1:]
                    p = stat.n[1:]
                    f = (1.0/sp-1.0/p)/(1-1.0/p)
                    f_std = 1.0/(1-1.0/p)/sp**2*sp_std
                    markers, caps, bars = ax1.errorbar(p, f, yerr=f_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                    [bar.set_alpha(alp) for bar in bars]
                    [cap.set_alpha(alp) for cap in caps]
                    tl = b[b.tags[['9000', node, ana]]]
                    stat_csv = tl.leaves().globfilter('stat*')[0]
                    stat = pd.read_csv(stat_csv.abspath)
                    t0 = stat.total[0]
                    speedup = t0/stat.total
                    speedup_std = stat.total_std*t0/(stat.total)**2+stat.total_std[0]/stat.total
                    sp = speedup[1:]
                    sp_std = speedup_std[1:]
                    p = stat.n[1:]
                    f = (1.0/sp-1.0/p)/(1-1.0/p)
                    f_std = 1.0/(1-1.0/p)/sp**2*sp_std
                    markers, caps, bars = ax2.errorbar(p, f, yerr=f_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                    [bar.set_alpha(alp) for bar in bars]
                    [cap.set_alpha(alp) for cap in caps]
            else:
                ts = b[b.tags[['900', s, sche, ana]]]
                sche_name = scheduler_full_name[sche]
                label = '{0}-{1}'.format(s, sche_name)
                stat_csv = ts.leaves().globfilter('stat*')[0]
                stat = pd.read_csv(stat_csv.abspath)
                if sche == 'distr' and s == 'SSD':
                    alp = 0.3
                else:
                    alp = 0.5
                t0 = stat.total[0]
                speedup = t0/stat.total
                speedup_std = stat.total_std*t0/(stat.total)**2+stat.total_std[0]/stat.total
                sp = speedup[1:]
                sp_std = speedup_std[1:]
                p = stat.n[1:]
                f = (1.0/sp-1.0/p)/(1-1.0/p)
                f_std = 1.0/(1-1.0/p)/sp**2*sp_std
                markers, caps, bars = ax1.errorbar(p, f, yerr=f_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                [bar.set_alpha(alp) for bar in bars]
                [cap.set_alpha(alp) for cap in caps]
                tl = b[b.tags[['9000', s, sche, ana]]]
                stat_csv = tl.leaves().globfilter('stat*')[0]
                stat = pd.read_csv(stat_csv.abspath)
                t0 = stat.total[0]
                speedup = t0/stat.total
                speedup_std = stat.total_std*t0/(stat.total)**2+stat.total_std[0]/stat.total
                sp = speedup[1:]
                sp_std = speedup_std[1:]
                p = stat.n[1:]
                f = (1.0/sp-1.0/p)/(1-1.0/p)
                f_std = 1.0/(1-1.0/p)/sp**2*sp_std
                markers, caps, bars = ax2.errorbar(p, f, yerr=f_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                [bar.set_alpha(alp) for bar in bars]
                [cap.set_alpha(alp) for cap in caps]


    handles, labels = ax1.get_legend_handles_labels()
    lgd = fig.legend(handles, labels, ncol=3, edgecolor=None, fontsize = 8, frameon=False, loc = 'upper center', bbox_to_anchor=(0.5, 0.1))
    text1 = ax1.text(-0.25, 1.0, 'A', fontsize = 20, transform=ax1.transAxes)
    text2 = ax2.text(-0.25, 1.0, 'B', fontsize = 20, transform=ax2.transAxes)
    plt.savefig('serial_fraction_{}.pdf'.format(ana.lower()), bbox_extra_artists=(lgd, text1, text2))
    plt.close()
