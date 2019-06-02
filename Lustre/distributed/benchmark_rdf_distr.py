from pmda import rdf
from pmda import rms
import MDAnalysis as mda
import time
import numpy as np
import dask
import sys
import shutil 
import os

PDB = "/oasis/scratch/comet/sfan19/temp_project/step5_znm.pdb"
XTC = "/oasis/scratch/comet/sfan19/temp_project/step7_90ns_center.xtc"

from dask.distributed import Client
Scheduler_IP = sys.argv[1]
print (Scheduler_IP)
print (type (Scheduler_IP))
c = Client(Scheduler_IP)

u = mda.Universe(PDB, XTC)
g1 = u.select_atoms('name OH2')
g2 = u.select_atoms('name OH2')
rdfs = rdf.InterRDF(g1, g2, range=(0,5))
t_time = []
t_io = []
t_compute = []
t_u = []
t_w = []
ns = [1,2,4,8,16,24,32,40,48,56,64,72]
for n in ns:
    rdfs.run(n_jobs=n, n_blocks=n)
    t_time.append([n, rdfs.timing.total, max(rdfs.timing.universe), 
                   np.max(rdfs.timing.wait),
                   rdfs.timing.prepare, np.sum(rdfs.timing.io)/(n),
                   np.sum(rdfs.timing.compute)/(n), rdfs.timing.conclude])
    t_io.append(rdfs.timing.io) 
    t_compute.append(rdfs.timing.compute)
    t_u.append(rdfs.timing.universe)
    t_w.append(rdfs.timing.wait)
    np.save('/oasis/scratch/comet/sfan19/temp_project/lustre/distributed/9000_Frames/RDF/3nodes/data2/universe_rdf_distr', t_u)
    np.save('/oasis/scratch/comet/sfan19/temp_project/lustre/distributed/9000_Frames/RDF/3nodes/data2/wait_rdf_distr', t_w)
    np.savetxt('/oasis/scratch/comet/sfan19/temp_project/lustre/distributed/9000_Frames/RDF/3nodes/data2/benchmark_rdf_distr.dat', np.array(t_time))
    np.savetxt('/oasis/scratch/comet/sfan19/temp_project/lustre/distributed/9000_Frames/RDF/3nodes/data2/io_rdf_distr.dat', np.array(t_io))
    np.savetxt('/oasis/scratch/comet/sfan19/temp_project/lustre/distributed/9000_Frames/RDF/3nodes/data2/compute_rdf_distr.dat', np.array(t_compute))
    print('Complete benchmarking for n={}'.format(n))
    time.sleep(5)

