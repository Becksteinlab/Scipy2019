from pmda import rdf
from pmda import rms
import MDAnalysis as mda
import time
import numpy as np
import dask
import sys
import shutil 
import os

path = '/oasis/scratch/comet/sfan19/temp_project/'
PDB = os.path.join(path, "YiiP_system.pdb")
XTC = os.path.join(path, "YiiP_system_9ns_center.xtc")

from dask.distributed import Client
Scheduler_IP = sys.argv[1]
print (Scheduler_IP)
print (type (Scheduler_IP))
c = Client(Scheduler_IP)

u = mda.Universe(PDB, XTC)
ca = u.select_atoms('name CA')
u.trajectory[0]
ref = u.select_atoms('name CA')
rmsd = rms.RMSD(ca, ref)
t_time = []
t_io = []
t_compute = []
t_u = []
t_w = []
ns = [1,2,4,8,16,24,32,40,48,56,64,72]
folder = '/oasis/scratch/comet/sfan19/temp_project/lustre/distributed/9000_Frames/RMS/6nodes/data4'

for n in ns:
    n = n-1
    rmsd.run(n_jobs=n+1, n_blocks=n+1)
    t_time.append([n+1, rmsd.timing.total, max(rmsd.timing.universe), 
                   np.max(rmsd.timing.wait),
                   rmsd.timing.prepare, np.sum(rmsd.timing.io)/(n+1),
                   np.sum(rmsd.timing.compute)/(n+1), rmsd.timing.conclude])
    t_io.append(rmsd.timing.io) 
    t_compute.append(rmsd.timing.compute)
    t_u.append(rmsd.timing.universe)
    t_w.append(rmsd.timing.wait)
    np.save(os.path.join(folder ,'universe_rms_distr'), t_u)
    np.save(os.path.join(folder, 'wait_rms_distr', t_w)
    np.savetxt(os.path.join(folder, 'benchmark_rms_distr.dat'), np.array(t_time))
    np.savetxt(os.path.join(folder, 'io_rms_distr.dat'), np.array(t_io))
    np.savetxt(os.path.join(folder, 'compute_rms_distr.dat'), np.array(t_compute))
    print('Complete benchmarking for n={}'.format(n+1))
    time.sleep(5)

