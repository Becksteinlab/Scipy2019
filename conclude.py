import datreant as dtr
import pandas as pd
from pmda.util import make_balanced_slices
import numpy as np
import glob
import os

def conclude(total, io, compute, univ, w):
    '''Function to gather the total time, prepare time,
    universe time and maximum IO, computation, and wait time.

    total: path to the total benchmark data
    io: path to the io data
    compute: pate to the compute data
    '''
    data_t = np.loadtxt(total)
    all_io = np.loadtxt(io)
    all_compute = np.loadtxt(compute)
    all_univese = np.load(univ)
    all_wait = np.load(w)
    ndx = data_t[:,0]
    total = data_t[:,1]
    prepare = data_t[:,4]
    conclude = data_t[:,7]
    universe = []
    wait = []
    IO = []
    com = []
    stop = len(all_io[1,:])-1
    start = 0
    
    for i, j in enumerate(ndx):
        n_frames = len(all_io[1,:])
        n_blocks = int(j)
        slices = make_balanced_slices(n_frames, n_blocks,
                                  start=start, stop=stop, step=1)
        io_block = np.zeros(n_blocks)
        com_block = np.zeros(n_blocks)
        for k, bslice in enumerate(slices):
            k_io = all_io[i, bslice.start:bslice.stop]
            k_com = all_compute[i, bslice.start:bslice.stop]
            io_block[k] = np.sum(k_io)
            com_block[k] = np.sum(k_com)
        main = all_univese[i]+all_wait[i]+io_block+com_block
        n = np.argmax(main)
        IO.append(io_block[n])
        com.append(com_block[n])
        universe.append(all_univese[i][n])
        wait.append(all_wait[i][n])
    d = {'n': ndx, 'total': total, 'prepare': prepare, 'conclude': conclude,
         'universe': universe, 'wait': wait, 'IO': IO, 'compute': com}
    df = pd.DataFrame(data = d)
    return df

#Gathering Treants
b = dtr.discover()

#tags for Treants
source = ['Lustre', 'SSD']
size = ['9000', '900']
analysis = ['RDF', 'RMS']
scheduler = ['distr', 'multi']
nodes = ['3nodes', '6nodes']


# for Lustre distributed
for n in size:
    print n
    for node in nodes:
        print node
        for ana in analysis:
            print ana
            results = pd.DataFrame([])
            t = b[b.tags[[n, node, ana]]]
            for data in t.trees():
                lvs = data.leaves()
                total = lvs.globfilter('benchmark*')[0].abspath
                io = lvs.globfilter('io*')[0].abspath
                compute = lvs.globfilter('compute*')[0].abspath
                univ = lvs.globfilter('universe*')[0].abspath
                w = lvs.globfilter('wait*')[0].abspath
                df = conclude(total, io, compute, univ, w)
                results = results.append(df)
            path = os.path.join(t[0].abspath, 'conclusion.csv')
            results.to_csv(path)
            
