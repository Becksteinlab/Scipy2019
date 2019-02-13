.. -*- mode: rst; fill-column: 9999; coding: utf-8 -*-

:author: Shujie Fan
:email: sfan19@asu.edu
:institution: Arizona State University

:author: Max Linke
:email: 
:institution: 

:author: Ioannis Paraskevakos
:email: i.paraskev@rutgers.edu
:institution: Rutgers University

:author: Richard J. Gowers
:email: richardjgowers@gmail.com
:institution: University of New Hampshire

:author: Michael Gecht
:email: michael.gecht@biophys.mpg.de
:institution: Max Planck Institute of Biophysics

:author: Oliver Beckstein
:email: obeckste@asu.edu 
:institution: Arizona State University 
:corresponding:

-------------------------------------------------------------------------

-------------------------------------------------------------------------

Short Summary
--------------

PMDA is a Python library that provides parallel analysis algorithms based on MDAnalysis. With a simple map-reduce scheme from Dask library, PMDA provides pre-defined analysis tasks and a common interface to create user-defined analysis taks. Although still in alpha stage, it is already used on resources ranging from multi-core laptops to XSEDE supercomputers to speed up analysis of molecular dynamics trajectories.


Brief abstract
--------------

MDAnalysis_ is an object-oriented Python library to analyze trajectories from molecular dynamics (MD) simulations in many popular formats [Gowers2016]_. With the development of highly optimized molecular dynamics software (MD) packages on HPC resources, the size of simulation trajectories is growing to terabyte size. Thus efficient analysis of MD simulations becomes a challenge for MDAnalysis, which does not yet provide a standard interface for parallel analysis. To address the challenge, we developed PMDA(https://www.mdanalysis.org/pmda/), a Python library that provides parallel analysis algorithms based on MDAnalysis.  PMDA parallelizes common analysis algorithms in MDAnalysis through a task-based approach with the Dask_ library [Rocklin2015]_.  We implement a simple map-reduce scheme for parallel trajectory analysis [Khoshlessan2017; Paraskevakos2018]. The trajectory is partitioned into blocks and analysis is performed separately and in parallel on each block (“map”). The results from each block are gathered and combined (“reduce”).  PMDA allows one to perform parallel trajectory analysis with pre-defined analysis tasks. In addition, it provides a common interface that makes it easy to create user-defined parallel analysis modules. PMDA supports all schedulers in Dask, and one can run analysis in a distributed fashion on HPC or ad-hoc clusters or on a single machine. We tested the performance of PMDA on single node and multiple nodes on local supercomputing resources and workstations. The results show that parallelization improves the performance of trajectory analysis. Although still in alpha stage, it is already used on resources ranging from multi-core laptops to XSEDE supercomputers to speed up analysis of molecular dynamics trajectories. PMDA is available under the GNU General Public License, version 2.


Keywords
--------
PMDA, MDAnalysis, High Performance Computing, Dask, Map-Reduce


Tested libraries
----------------

- MDAnalysis_ 0.19.2
- Dask_ 1.1.1
- PMDA 0.2.0


References
----------

.. [Gowers2016] R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler, D. L. Dotson, J. Domański, S. Buchoux, I. M. Kenney, and O. Beckstein. MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. In S. Benthall and S. Rostrup, editors, Proceedings of the 15th Python in Science Conference, pages 102 – 109, Austin, TX, 2016. SciPy. URL http://mdanalysis.org.

.. [Rocklin2015] Rocklin, Matthew. Dask: Parallel Computation with Blocked algorithms and Task Scheduling. Presented at the Python in Science Conference, Austin, Texas, pp. 126–132. https://doi.org/10.25080/Majora-7b98e3ed-013

.. [Khoshlessan2017] Khoshlessan, Mahzad; Beckstein, Oliver (2017): Parallel analysis in the MDAnalysis Library: Benchmark of Trajectory File Formats. figshare. doi:`10.6084/m9.figshare.4695742`_

.. [Paraskevakos2018] Ioannis Paraskevakos, Andre Luckow, Mahzad Khoshlessan, George Chantzialexiou, Thomas E. Cheatham, Oliver Beckstein, Geoffrey C. Fox and Shantenu Jha. Task- parallel Analysis of Molecular Dynamics Trajectories. In 47th International Conference on Parallel Processing (ICPP 2018). doi: 10.1145/3225058.3225128.

.. _MDAnalysis: http://mdanalysis.org
.. _Dask: http://dask.pydata.org
.. _PMDA: https://www.mdanalysis.org/pmda/
