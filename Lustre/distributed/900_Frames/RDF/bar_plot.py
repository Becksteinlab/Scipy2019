import numpy as np
import matplotlib.pyplot as plt
import seaborn.apionly as sns
import matplotlib.gridspec as gridspec

a = np.loadtxt('benchmark_rdf_distr.dat')
ndx = a[:,0]
total = a[:,1]
universe = a[:,2]
wait = a[:,3]
prepare = a[:,4]
IO = np.loadtxt('io_rdf_distr.dat')
compute = np.loadtxt('compute_rdf_distr.dat')
conclude = a[:,7]

ind = np.arange(len(total))+1
width = 0.35
n = 0
ns = [1,2,4,8,16,24,32,40,48,56,64,72]
Dom_com = np.zeros(12)
Dom_IO = np.zeros(12)
stop = len(IO[1,:])-1
for i in range(12):
    n_frames = len(IO[1,:])
    n_blocks = ns[i]
    bsize = int(np.ceil(n_frames / float(n_blocks)))
    print(bsize)
    ith_IO = np.zeros(i+1)
    ith_com = np.zeros(i+1)
    for j in range(i+1):
        jth_IO = IO[i, (bsize*j):min(stop, (j+1)*bsize)]
        jth_com = compute[i,(bsize*j):min(stop, (j+1)*bsize)]
        ith_IO[j] = np.sum(jth_IO)
        ith_com[j] = np.sum(jth_com)

    ith_Dom = ith_IO+ith_com
    n = np.argmax(ith_Dom)
    Dom_com[i] = ith_com[n]
    Dom_IO[i] = ith_IO[n]

fig = plt.figure(figsize=(10, 6), tight_layout=True)
gs = gridspec.GridSpec(3, 3)
ax1 = fig.add_subplot(gs[0:2,:])
sns.despine(offset=10, ax=ax1)
rects1 = ax1.bar(ind - width/2, total, width, color='SkyBlue', label='total')
rects2 = ax1.bar(ind + width/2, universe, width, color='m', label='universe')
rects0 = ax1.bar(ind + width/2, wait, width, color='g', bottom=universe, label='wait')
rects3 = ax1.bar(ind + width/2, prepare, width, color='b', bottom=universe+wait, label='prepare')
rects4 = ax1.bar(ind + width/2, Dom_IO, width, color='navy', bottom=universe+prepare+wait, label='I/O')
rects5 = ax1.bar(ind + width/2, Dom_com, width, color='r', bottom=universe+prepare+Dom_IO+wait, label='compute')
rects6 = ax1.bar(ind + width/2, conclude, width, color='darkorange', bottom=universe+prepare+Dom_IO+Dom_com+wait, label='conclude')
ax1.legend(bbox_to_anchor=(0.5, 1.0), loc='upper center', ncol=3, edgecolor=None)
ax1.set_xticks(ind)
ax1.set_xticklabels(ns)
ax1.set_xlabel('Number of cores')
ax1.set_ylabel('Time(s)')
sns.despine(offset=10, ax=ax1)

ax1 = fig.add_subplot(gs[2,0])
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_adjustable("datalim")
ax1.plot(ndx, total)
ax1.set_xlabel('Number of cores')
ax1.set_ylabel('Tatol Time(s)')

ax1 = fig.add_subplot(gs[2,1])
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_adjustable("datalim")
ax1.plot(ndx, Dom_com)
ax1.set_xlabel('Number of cores')
ax1.set_ylabel('Computation Time(s)')
plt.suptitle('Comet lustre', y=0.995, fontsize=16)

ax1 = fig.add_subplot(gs[2,2])
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_adjustable("datalim")
ax1.plot(ndx, Dom_IO)
ax1.set_xlabel('Number of cores')
ax1.set_ylabel('I/O Time(s)')

#plt.tight_layout()
plt.savefig('bcmk.pdf')
