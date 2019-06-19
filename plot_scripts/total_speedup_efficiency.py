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
    fig = plt.figure(figsize=(12, 8), tight_layout=True)
    gs = gridspec.GridSpec(5, 9)

    ax0 = fig.add_subplot(gs[0:2,0:3])
    ax0.set_xlabel('Number of cores')
    ax0.set_ylabel('total time (s)')
    sns.despine(offset=10, ax=ax0)

    ax0.set_title('pmda.{} on 900 frames'.format(ana.lower()))
    ax0.set_xscale("log")
    ax0.set_yscale("log")

    ax1 = fig.add_subplot(gs[0:2,3:6])
    ax1.set_xlabel('Number of cores')
    ax1.set_ylabel('Efficiency')
    ax1.set_ylim([0,1])
    sns.despine(offset=10, ax=ax1)

    ax1.set_title('pmda.{} on 900 frames'.format(ana.lower()))

    ax2 = fig.add_subplot(gs[0:2,6:9])

    ax2.set_xlabel('Number of cores')
    ax2.set_ylabel('Speedup')
    ax2.set_title('pmda.{} on 900 frames'.format(ana.lower()))
    sns.despine(offset=10, ax=ax2)

    ax5 = fig.add_subplot(gs[2:4,0:3])

    ax5.set_xlabel('Number of cores')
    ax5.set_ylabel('total time (s)')
    ax5.set_title('pmda.{} on 9000 frames'.format(ana.lower()))
    sns.despine(offset=10, ax=ax5)
    ax5.set_xscale("log")
    ax5.set_yscale("log")

    ax3 = fig.add_subplot(gs[2:4,3:6])
    ax3.set_xlabel('Number of cores')
    ax3.set_ylabel('Efficiency')
    ax3.set_ylim([0,1])
    sns.despine(offset=10, ax=ax3)

    ax3.set_title('pmda.{} on 9000 frames'.format(ana.lower()))


    ax4 = fig.add_subplot(gs[2:4,6:9])

    ax4.set_xlabel('Number of cores')
    ax4.set_ylabel('Speedup')

    ax4.set_title('pmda.{} on 9000 frames'.format(ana.lower()))
    sns.despine(offset=10, ax=ax4)
    for s in source:
        for sche in scheduler:
            if (sche == 'distr') & (s == 'Lustre'):
                for node in nodes:
                    print(node)
                    ts = b[b.tags[['900', node, ana]]]
                    stat_csv = ts.leaves().globfilter('stat*')[0]
                    stat = pd.read_csv(stat_csv.abspath)
                    label = 'Lustre-distributed-{}'.format(node)
                    markers0, caps0, bars0 = ax0.errorbar(stat.n, stat.total, yerr=stat.total_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                    [bar.set_alpha(0.5) for bar in bars0]
                    [cap.set_alpha(0.5) for cap in caps0]
                    t0 = stat.total[0]
                    speedup = t0/stat.total
                    speedup_std = stat.total_std*t0/(stat.total)**2+stat.total_std[0]/stat.total
                    eff = speedup/stat.n
                    eff_std = speedup_std/stat.n
                    markers, caps, bars = ax1.errorbar(stat.n, eff, yerr=eff_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                    [bar.set_alpha(0.5) for bar in bars]
                    [cap.set_alpha(0.5) for cap in caps]
                    markers, caps, bars = ax2.errorbar(stat.n, speedup, yerr=speedup_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                    [bar.set_alpha(0.5) for bar in bars]
                    [cap.set_alpha(0.5) for cap in caps]
                    ax2.plot(stat.n, stat.n, '--', color = 'black')
                    tl = b[b.tags[['9000', node, ana]]]
                    stat_csv = tl.leaves().globfilter('stat*')[0]
                    stat = pd.read_csv(stat_csv.abspath)
                    markers5, caps5, bars5 = ax5.errorbar(stat.n, stat.total, yerr=stat.total_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                    [bar.set_alpha(0.5) for bar in bars5]
                    [cap.set_alpha(0.5) for cap in caps5]
                    t0 = stat.total[0]
                    speedup = t0/stat.total
                    speedup_std = stat.total_std*t0/(stat.total)**2+stat.total_std[0]/stat.total
                    eff = speedup/stat.n
                    eff_std = speedup_std/stat.n
                    markers, caps, bars = ax3.errorbar(stat.n, eff, yerr=eff_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                    [bar.set_alpha(0.5) for bar in bars]
                    [cap.set_alpha(0.5) for cap in caps]
                    markers, caps, bars = ax4.errorbar(stat.n, speedup, yerr=speedup_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                    [bar.set_alpha(0.5) for bar in bars]
                    [cap.set_alpha(0.5) for cap in caps]
                    ax4.plot(stat.n, stat.n, '--', color = 'black')
            else:
                ts = b[b.tags[['900', s, sche, ana]]]
                sche_name = scheduler_full_name[sche]
                label = '{0}-{1}'.format(s, sche_name)
                stat_csv = ts.leaves().globfilter('stat*')[0]
                stat = pd.read_csv(stat_csv.abspath)
                markers0, caps0, bars0 = ax0.errorbar(stat.n, stat.total, yerr=stat.total_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                if sche == 'distr' and s == 'SSD':
                    alp = 0.3
                else:
                    alp = 0.5
                [bar.set_alpha(alp) for bar in bars0]
                [cap.set_alpha(alp) for cap in caps0]
                t0 = stat.total[0]
                speedup = t0/stat.total
                speedup_std = stat.total_std*t0/(stat.total)**2+stat.total_std[0]/stat.total
                eff = speedup/stat.n
                eff_std = speedup_std/stat.n
                markers, caps, bars = ax1.errorbar(stat.n, eff, yerr=eff_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                [bar.set_alpha(alp) for bar in bars]
                [cap.set_alpha(alp) for cap in caps]
                markers, caps, bars = ax2.errorbar(stat.n, speedup, yerr=speedup_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                [bar.set_alpha(alp) for bar in bars]
                [cap.set_alpha(alp) for cap in caps]
                ax2.plot(stat.n, stat.n, '--', color = 'black')
                tl = b[b.tags[['9000', s, sche, ana]]]
                stat_csv = tl.leaves().globfilter('stat*')[0]
                stat = pd.read_csv(stat_csv.abspath)
                markers5, caps5, bars5 = ax5.errorbar(stat.n, stat.total, yerr=stat.total_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                [bar.set_alpha(alp) for bar in bars5]
                [cap.set_alpha(alp) for cap in caps5]
                t0 = stat.total[0]
                speedup = t0/stat.total
                speedup_std = stat.total_std*t0/(stat.total)**2+stat.total_std[0]/stat.total
                eff = speedup/stat.n
                eff_std = speedup_std/stat.n
                markers, caps, bars = ax3.errorbar(stat.n, eff, yerr=eff_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                [bar.set_alpha(alp) for bar in bars]
                [cap.set_alpha(alp) for cap in caps]
                markers, caps, bars = ax4.errorbar(stat.n, speedup, yerr=speedup_std, label=label, markeredgewidth=1, capsize=2, markersize=4)
                [bar.set_alpha(alp) for bar in bars]
                [cap.set_alpha(alp) for cap in caps]
                ax4.plot(stat.n, stat.n, '--', color = 'black')

    handles, labels = ax1.get_legend_handles_labels()
    lgd = fig.legend(handles, labels, ncol=3, edgecolor=None, frameon=False, loc = 'upper center', bbox_to_anchor=(0.5, 0.1))
    text0 = ax0.text(-0.25, 1.0, 'A', fontsize = 20, transform=ax0.transAxes)
    text1 = ax1.text(-0.25, 1.0, 'B', fontsize = 20, transform=ax1.transAxes)
    text2 = ax2.text(-0.25, 1.0, 'C', fontsize = 20, transform=ax2.transAxes)
    text5 = ax5.text(-0.25, 1.0, 'D', fontsize = 20, transform=ax5.transAxes)
    text3 = ax3.text(-0.25, 1.0, 'E', fontsize = 20, transform=ax3.transAxes)
    text4 = ax4.text(-0.25, 1.0, 'F', fontsize = 20, transform=ax4.transAxes)
    plt.savefig('Total_Eff_SU_{}.pdf'.format(ana.lower()), bbox_extra_artists=(lgd,text0, text1, text2, text5, text3, text4,))
    plt.close()
