# Topological susceptibility in Yang-Mills theories
and
# Sp(2N) Yang-Mills theories on the lattice: scale setting and topology

## Analyis code release

This repository contains the analysis code used to prepare the plots and tables
included in [Topological susceptibility in Yang-Mills theories][short-paper] and
[Sp(2N) Yang-Mills theories on the lattice: scale setting and topology][long-paper].

## Health warning

Much of the code in this repository was written as "disposable" code specifically to
get the specific outputs needed for the above papers, and is not intended to be reusable
or generalisable to other data sets.
It does not (and is not intended to) provide an example of "good software engineering practice".
In some cases the code was edited by hand during the analysis phase, and in most cases the 
individual tools were run separately by hand, with the management of the flow of data being done
manually. The automation of processes with `make` was added after the fact. 

## Requirements

The analysis is primarily written in Python, with the help of a few small Bash scripts.
The various components are joined together using Make. The code has been tested with:

* Python 3.10 (with requirements as documented in `environment.yml`)
* GNU Make 4.3
* LaTeX (to allow Matplotlib to put formatted equations into plots)
* (Optionally) Conda to manage dependencies

In particular, versions of Make before 4.3 will not work. This includes the
version of Make installed on macOS by default.

## Setup

### Installation

* Download the repository
* Run

      conda env create -f environment.yml

  to create a new Conda environment including the requisite packages.
  (Alternatively, create a Python environment with the listed packages, and
  ensure that you have GNU Make 4.3 or later.)

### Data

* Download and extract the data. Specifically, you need the `data` directory,
  provided as a ZIP file in the accompanying data release. By default it is
  expected that this is placed in the root directory of the repository,
  but this can be changed in the `Makefile`.

### Running the analysis

* With the software and data downloaded, you should be able to reproduce the
  full analysis by typing

      make

  Plots will be created in the `plots` directory, and tables in the `tables`
  directory. Intermediary files not included in the publication are held in
  the `processed_data` directory. Each of these will be created by `make`
  automatically.
* The analysis will run in parallel by using `make -j N`, where `N` is replaced
  by the number of parallel jobs to use.
* On four cores of an Apple M1 CPU, the full analysis takes around an hour.

## References

As this analysis makes use of data from other publications that has not been
previously made available in an open, machine-readable way, this repository
quotes results from a number of other publications. All quoted data is in the
`quoted_data` directory. Specifically, the files there borrow data from:

* `clim_SU_AA_MT.dat`: SU(N) gauge theories in 3+1 dimensions: glueball spectrum, string tensions and topology, Andreas Athenodorou, Michael Teper, [JHEP 12 (2021) 082][athenodorou-teper]. [arXiv:2106.00364](https://arxiv.org/2106.00364)
* `clim_SU_BL_MT.dat`: SU(N) gauge theories in four-dimensions: Exploring the approach to N = infinity, B. Lucini, M. Teper, [JHEP 06 (2001) 050](https://doi.org/10.1088/1126-6708/2001/06/050). [arXiv:hep-lat/0103027](https://arxiv.org/abs/hep-lat/0103027)
* `clim_SU_CB_CB_MD.dat`: Large-N SU(N) Yang-Mills theories with milder topological freezing, Claudio Bonanno, Claudio Bonati, Massimo D'Elia, [JHEP 03 (2021) 111](https://doi.org/10.1007/JHEP03(2021)111). [arXiv:2012.14000](https://arxiv.org/abs/2012.14000)
* `clim_SU_CB_MD.dat`: θ dependence of 4D SU(N) gauge theories in the large-N limit, Claudio Bonati, Massimo D'Elia, Paolo Rossi, Ettore Vicari, [Phys.Rev.D 94 (2016) 8, 085017](https://arxiv.org/abs/1607.06360). [arXiv:1607.06360](https://arxiv.org/abs/1607.06360)
* `clim_SU_DD_EV_HP.dat`: Theta dependence of SU(N) gauge theories, Luigi Del Debbio, Haralambos Panagopoulos, Ettore Vicari (Apr, 2002), [JHEP 08 (2002) 044](https://doi.org/10.1088/1126-6708/2002/08/044). [arXiv:hep-th/0204125](https://arxiv.org/abs/hep-th/0204125)
* `sqrts_vs_beta.dat`: Glueballs and strings in Sp(2N) Yang-Mills theories, Ed Bennett, Jack Holligan, Deog Ki Hong, Jong-Wan Lee, C.-J. David Lin, Biagio Lucini, Maurizio Piai, Davide Vadacchino, [Phys.Rev.D 103 (2021) 5, 054509](https://doi.org/10.1103/PhysRevD.103.054509). [arXiv:2010.15781](https://arxiv.org/abs/2010.15781)


[long-paper]: https://arxiv.org/abs/2205.09364
[short-paper]: https://arxiv.org/abs/2205.09254
[athenodorou-teper]: https://doi.org/10.1007/JHEP12(2021)082
