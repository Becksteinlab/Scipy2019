from pmda import rdf
from pmda import rms
import MDAnalysis as mda
import time
import numpy as np
import dask
import sys
import shutil 
import os

dask.config.set(scheduler='multiprocessing')
folder = sys.argv[1]
PDB = os.path.join(folder, 'test.pdb')
XTC = os.path.join(folder, 'test.xtc')

shutil.copyfile("/oasis/scratch/comet/sfan19/temp_project/YiiP_system.pdb", PDB)
shutil.copyfile("/oasis/scratch/comet/sfan19/temp_project/YiiP_system_9ns_center.xtc", XTC)

u = mda.Universe(PDB, XTC)
g1 = u.select_atoms('name OH2')
g2 = u.select_atoms('name OH2')
rdfs = rdf.InterRDF(g1, g2, range=(0,5))
t_time = []
t_io = []
t_compute = []
t_u = []
t_w = []
for n in range(24):
    rdfs.run(n_jobs=n+1, n_blocks=n+1)
    t_time.append([n+1, rdfs.timing.total, max(rdfs.timing.universe), 
                   np.max(rdfs.timing.wait),
                   rdfs.timing.prepare, np.sum(rdfs.timing.io)/(n+1),
                   np.sum(rdfs.timing.compute)/(n+1), rdfs.timing.conclude])
    t_io.append(rdfs.timing.io) 
    t_compute.append(rdfs.timing.compute)
    t_u.append(rdfs.timing.universe)
    t_w.append(rdfs.timing.wait)
    np.save('/oasis/scratch/comet/sfan19/temp_project/SSD/multiprocessing/900_Frames/RDF/data4/universe_rdf_multi.npy', t_u)
    np.save('/oasis/scratch/comet/sfan19/temp_project/SSD/multiprocessing/900_Frames/RDF/data4/wait_rdf_multi.npy', t_w)
    np.savetxt('/oasis/scratch/comet/sfan19/temp_project/SSD/multiprocessing/900_Frames/RDF/data4/benchmark_rdf_multi.dat', np.array(t_time))
    np.savetxt('/oasis/scratch/comet/sfan19/temp_project/SSD/multiprocessing/900_Frames/RDF/data4/io_rdf_multi.dat', np.array(t_io))
    np.savetxt('/oasis/scratch/comet/sfan19/temp_project/SSD/multiprocessing/900_Frames/RDF/data4/compute_rdf_multi.dat', np.array(t_compute))
    print('Complete benchmarking for n={}'.format(n+1))

os.remove(XTC)
os.remove(PDB)
