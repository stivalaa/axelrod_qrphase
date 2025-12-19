# Another phase transition in the Axelrod model

## Software

Imported from http://munk.cis.unimelb.edu.au/~stivalaa/axelrod_qrphase/

Also available with Zenodo DOI: [![DOI](https://zenodo.org/badge/132081707.svg)](https://doi.org/10.5281/zenodo.17992801)

This software is free under the terms of the GNU General Public License.
Some is derived from code developed for an earlier publication
[Ultrametric distribution of culture vectors in an extended Axelrod model of cultural dissemination](http://munk.cis.unimelb.edu.au/~stivalaa/ultrametric_axelrod/).
The simulation code uses Python
and parallelization using MPI (with [mpi4py](http://mpi4py.scipy.org/)). It also requires the Python libraries [NumPy](http://www.numpy.org/) (part of the [SciPy package](http://www.scipy.org/)).
The code for numerical integration of the mean-field approximation 
is written in MATLAB, and was run in MATLAB R2014a.

The Python code was run with NumPy version 1.9.1, SciPy version 0.14.1, igraph version 0.6 and mpi4py version 1.3.1 under Python version 2.7.9 on a Lenovo x86 cluster (992 Intel Haswell compute cores running at 2.3GHz) running Linux (RHEL 6) with Open MPI version 1.10.0.
The C++ code was compiled with gcc version 4.9.2. 

### Running the models

The model can be run with a command line such as: `mpirun --mca mpi_warn_on_fork 0 python ./lattice-python-mpi/src/axelrod/geo/expphysicstimeline/multiruninitmain.py m:100 F:5 ./lattice-dyadic-noise-cpp-end/model reinit end 10000`

## Reference

If you use our software, data, or results in your research, please cite:

- A. Stivala and P. Keeler 2016. [Another phase transition in the Axelrod model](https://arxiv.org/abs/1612.02537) arXiv preprint [arXiv:1612.02537](https://arxiv.org/abs/1612.02537).

