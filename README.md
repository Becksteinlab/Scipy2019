# Data files and scripts for "PMDA – Parallel Molecular Dynamics Analysis"
[![DOI](https://zenodo.org/badge/170560828.svg)](https://zenodo.org/badge/latestdoi/170560828)


This repository contains data files and scripts for our Scipy2019 paper on [PMDA](https://www.mdanalysis.org/pmda/).


## Summary

PMDA is a Python library that provides parallel analysis algorithms based on MDAnalysis. With a simple map-reduce scheme from Dask library, PMDA provides pre-defined analysis tasks and a common interface to create user-defined analysis taks. Although still in alpha stage, it is already used on resources ranging from multi-core laptops to XSEDE supercomputers to speed up analysis of molecular dynamics trajectories.


## Abstract

[MDAnalysis](https://mdanalysis.org) is an object-oriented Python library to analyze trajectories from molecular dynamics (MD) simulations in many popular formats [Gowers2016]. With the development of highly optimized molecular dynamics software (MD) packages on HPC resources, the size of simulation trajectories is growing to terabyte size. Thus efficient analysis of MD simulations becomes a challenge for MDAnalysis, which does not yet provide a standard interface for parallel analysis. To address the challenge, we developed PMDA (https://www.mdanalysis.org/pmda/), a Python library that provides parallel analysis algorithms based on MDAnalysis.  PMDA parallelizes common analysis algorithms in MDAnalysis through a task-based approach with the [Dask](https://dask.org) library [Rocklin2015].  We implement a simple map-reduce scheme for parallel trajectory analysis [Khoshlessan2017]  [Paraskevakos2018]. The trajectory is partitioned into blocks and analysis is performed separately and in parallel on each block (“map”). The results from each block are gathered and combined (“reduce”).  PMDA allows one to perform parallel trajectory analysis with pre-defined analysis tasks. In addition, it provides a common interface that makes it easy to create user-defined parallel analysis modules. PMDA supports all schedulers in Dask, and one can run analysis in a distributed fashion on HPC or ad-hoc clusters or on a single machine. We tested the performance of PMDA on single node and multiple nodes on local supercomputing resources and workstations. The results show that parallelization improves the performance of trajectory analysis. Although still in alpha stage, it is already used on resources ranging from multi-core laptops to XSEDE supercomputers to speed up analysis of molecular dynamics trajectories. PMDA is available under the GNU General Public License, version 2.


Keywords
--------
PMDA, MDAnalysis, High Performance Computing, Dask, Map-Reduce


References
----------


[Gowers2016] Gowers, Richard J.; Linke, Max; Barnoud, Jonathan; Reddy, Tyler J. E.; Melo, Manuel N.; Seyler, Sean L.; Dotson, David L.; Domański, Jan; Buchoux, Sébastien; Kenney, Ian M.; and Beckstein, Oliver. MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. In S. Benthall and S. Rostrup, editors, Proceedings of the 15th Python in Science Conference, pages 102 – 109, Austin, TX, 2016. SciPy. URL: https://www.mdanalysis.org/.

[Rocklin2015] Rocklin, Matthew. Dask: Parallel Computation with Blocked algorithms and Task Scheduling. Presented at the Python in Science Conference, Austin, Texas, pp. 126–132. DOI: https://doi.org/10.25080/Majora-7b98e3ed-013.

[Khoshlessan2017] Khoshlessan, Mahzad; Beckstein, Oliver (2017): Parallel analysis in the MDAnalysis Library: Benchmark of Trajectory File Formats. figshare. DOI: https://doi.org/10.6084/m9.figshare.4695742.

[Paraskevakos2018] Ioannis Paraskevakos, Andre Luckow, Mahzad Khoshlessan, George Chantzialexiou, Thomas E. Cheatham, Oliver Beckstein, Geoffrey C. Fox and Shantenu Jha. Task- parallel Analysis of Molecular Dynamics Trajectories. In 47th International Conference on Parallel Processing (ICPP 2018). DOI: https://doi.org/10.1145/3225058.3225128.

