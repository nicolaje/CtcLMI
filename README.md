# Contractor for Linear Matrix Inequalities

This repository contains the source code of :

- the contractor proposed for the article "Contractors and Linear Matrix Inequalities" by Nicola & Jaulin published in [ASCE-ASME Journal of Risk and Uncertainty in Engineering Systems, Part B: Mechanical Engineering](http://ascelibrary.org/doi/abs/10.1115/1.4030781) available [here](https://www.ensta-bretagne.fr/jaulin/paper_nicola_lmi.pdf).
- the examples that illustrate the article

## Installation

### Requirements

The contractor relies on :

- the [IBEX](http://ibex-lib.org/) library for interval computations and constraints processing. It has been developped with IBEX 2.1.5, and might require some adjustments to work with the latest version of IBEX (particularly for specifying the dimension in the constructor).
- the [SDPA](http://sdpa.sourceforge.net/) library as an SDP solver

### Setup

This project assumes the IBEX, SOPLEX, SDPA, MUMPS, PORD, MPISEQ, BLAS libraries are available respectively in the ibex-dev/lib, soplex-1.7.2/lib, sdpa/lib, sdpa/share/sdpa/mumps/build/lib/, sdpa/share/sdpa/mumps/build/lib/, sdpa/share/sdpa/mumps/build/libseq/, sdpa/share/sdpa/blas folders. Check the CMakeLists.txt file for additional information.

The project has been developped in a Linux environment and might require some modifications (especially regarding building SDPA, although other solvers can be used) to work on other platforms.
