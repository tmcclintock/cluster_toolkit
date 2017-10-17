# Cluster_WL
This repository contains routines used in the Dark Energy Survey Year 1 stacked cluster weak lensing analysis. These routines form the basis of the model of cluster lensing profiles.

## Requirements
* Gnu Scientific Library [GSL](https://www.gnu.org/software/gsl/)

Linking to the GSL can either be done through the gsl-config, in which case the links are automatically identified, or manually by setting a path to ```gsl/lib``` with ```GSLL``` and a path to ```gsl/include``` with ```GSLI```. This code is confirmed to work with GSL version 2.4.

## Installation
The code is installed automatically in your Python path by running ```python setup.py install```.

## Functionality
As is, this package is only a Python module and does not contain a linkable C library. A re-implementation in the [DESC Core Cosmology Library](https://github.com/LSSTDESC/CCL) is underway, which will contain access to these funcitons in C. No gaurantee is made the the CCL implementation will be up to date with this module.

The functions in this module are in three groups: correlation functions in ```clusterwl.xi```, the linear bias calculated in ```clusterwl.bias```, surface mass density profiles and differential profiles are in ```clusterwl.deltasigma```, miscentering effects are in ```clusterwl.miscentering```, and finally radial bin averaged quantities are calculated in ```clusterwl.averaging```. Full documentation for all of these can be accessed in the docstrings, and examples of how to call these can be found in the notebooks in the ```examples``` directory.