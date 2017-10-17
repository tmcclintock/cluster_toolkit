# Cluster_WL
This repository contains routines used in the Dark Energy Survey Year 1 stacked cluster weak lensing analysis. These routines form the basis of the model of cluster lensing profiles.

## Requirements
* Gnu Scientific Library [GSL](https://www.gnu.org/software/gsl/)

Linking to the GSL can either be done through the gsl-config, in which case the links are automatically identified, or manually by setting a path to ```gsl/lib``` with ```GSLL``` and a path to ```gsl/include``` with ```GSLI```. This code is confirmed to work with GSL version 2.4.

## Installation
The code is installed automatically in your python path by running ```python setup.py install```.