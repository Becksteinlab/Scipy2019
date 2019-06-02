from pmda import rdf
from pmda import rms
import MDAnalysis as mda
import time
import numpy as np
import dask
import sys
import shutil 
import os

from dask.distributed import Client
Scheduler_IP = sys.argv[2]
print (Scheduler_IP)
print (type (Scheduler_IP))
c = Client(Scheduler_IP)

folder = sys.argv[1]
PDB = os.path.join(folder, 'test.pdb')
XTC = os.path.join(folder, 'test.xtc')

shutil.copyfile("/oasis/scratch/comet/sfan19/temp_project/YiiP_system.pdb", PDB)
shutil.copyfile("/oasis/scratch/comet/sfan19/temp_project/YiiP_system_90ns_center.xtc", XTC)

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
for n in range(24):
    rmsd.run(n_jobs=n+1, n_blocks=n+1)
    t_time.append([n+1, rmsd.timing.total, max(rmsd.timing.universe), 
                   np.max(rmsd.timing.wait),
                   rmsd.timing.prepare, np.sum(rmsd.timing.io)/(n+1),
                   np.sum(rmsd.timing.compute)/(n+1), rmsd.timing.conclude])
    t_io.append(rmsd.timing.io) 
    t_compute.append(rmsd.timing.compute)
    t_u.append(rmsd.timing.universe)
    t_w.append(rmsd.timing.wait)
    np.save('/oasis/scratch/comet/sfan19/temp_project/SSD/distributed/9000_Frames/RMS/data3/universe_rms_distr', t_u)
    np.save('/oasis/scratch/comet/sfan19/temp_project/SSD/distributed/9000_Frames/RMS/data3/wait_rms_distr', t_w)
    np.savetxt('/oasis/scratch/comet/sfan19/temp_project/SSD/distributed/9000_Frames/RMS/data3/benchmark_rms_distr.dat', np.array(t_time))
    np.savetxt('/oasis/scratch/comet/sfan19/temp_project/SSD/distributed/9000_Frames/RMS/data3/io_rms_distr.dat', np.array(t_io))
    np.savetxt('/oasis/scratch/comet/sfan19/temp_project/SSD/distributed/9000_Frames/RMS/data3/compute_rms_distr.dat', np.array(t_compute))
    print('Complete benchmarking for n={}'.format(n+1))

os.remove(XTC)
os.remove(PDB)
