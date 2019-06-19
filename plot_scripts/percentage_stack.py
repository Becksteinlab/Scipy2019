import numpy as np
import matplotlib.pyplot as plt
import seaborn.apionly as sns
import matplotlib.gridspec as gridspec
import datreant as dtr
import pandas as pd


b = dtr.discover()
analysis = ['RDF', 'RMS']
size = ['900', '9000']

for ana in analysis:
    fig = plt.figure(figsize=(6, 5), tight_layout=True)
    gs = gridspec.GridSpec(5, 1)

    #900
    ax1 = fig.add_subplot(gs[0:2,0])
    ax1.set_xlabel('Number of cores')
    ax1.set_ylabel('Time fraction (%)')
    ax1.set_title('pmda.{} on 900 frames'.format(ana.lower()))
    sns.despine(offset=10, ax=ax1)
    ts = b[b.tags[[size[0], '3nodes', ana]]]
    stat_csv = ts.leaves().globfilter('stat*')[0]
    stat = pd.read_csv(stat_csv.abspath)
    t = stat.prepare+stat.universe+stat.compute+stat.IO+stat.wait+stat.conclude
    t = t*0.01
    x = stat.n
    y = [stat.prepare/t, stat.universe/t, stat.compute/t, stat.IO/t, stat.wait/t, stat.conclude/t]
    labels = ['prepare', 'universe', 'compute', 'I/O', 'wait', 'conclude']
    ax1.stackplot(x, y, labels=labels)

    # 9000
    ax2 = fig.add_subplot(gs[2:4,0])
    ax2.set_xlabel('Number of cores')
    ax2.set_ylabel('Time fraction (%)')
    ax2.set_title('pmda.{} on 9000 frames'.format(ana.lower()))
    sns.despine(offset=10, ax=ax2)
    ts = b[b.tags[[size[1], '3nodes', ana]]]
    stat_csv = ts.leaves().globfilter('stat*')[0]
    stat = pd.read_csv(stat_csv.abspath)
    t = stat.prepare+stat.universe+stat.compute+stat.IO+stat.wait+stat.conclude
    t = t*0.01
    x = stat.n
    y = [stat.prepare/t, stat.universe/t, stat.compute/t, stat.IO/t, stat.wait/t, stat.conclude/t]
    labels = ['prepare', 'universe', 'compute', 'I/O', 'wait', 'conclude']
    ax2.stackplot(x, y, labels=labels)
    handles, labels = ax1.get_legend_handles_labels()
    lgd = fig.legend(handles, labels, ncol=3, edgecolor=None, frameon=False, loc = 'upper center', bbox_to_anchor=(0.5, 0.12))
    text1 = ax1.text(-0.2, 1.0, 'A', fontsize = 16, transform=ax1.transAxes)
    text2 = ax2.text(-0.2, 1.0, 'B', fontsize = 16, transform=ax2.transAxes)
    plt.savefig('percentage_stack_{}.pdf'.format(ana.lower()), bbox_extra_artists=(lgd, text1, text2,))
    plt.close()
