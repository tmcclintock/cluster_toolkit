******************************
Cluster Toolkit Documentation
******************************

Cluster Toolkit is a Python package specifically built for calculating weak lensing signals from galaxy clusters and cluster cosmology.
The core functionality of the package includes:

 * 3D correlation functions :math:`\xi(r)`
 * 3D density functions :math:`\rho(r)`
 * Halo bias models :math:`b(M)`
 * Projected density and differential profiles :math:`\Sigma(R)` and :math:`\Delta\Sigma`
 * Boost factor models :math:`\mathcal{B} = (1-f_{\rm cl})^{-1}`
 * Miscentering affects on projected profiles :math:`f_{\rm mis}`
 * Halo mass functions :math:`\frac{dn}{dM}(M,z)` (in development)
 * Sunyaev-Zel'dovich (SZ) cluster signals :math:`Y_{SZ}` (in development)
 * Cluster magnification :math:`\kappa(\theta)` and shear profiles :math:`\gamma(\theta)` (in development)

The source code is publically available at https://github.com/tmcclintock/cluster_toolkit.

.. note::
   Unless stated otherwise, all distances are assumed to be :math:`{\rm Mpc}/h` comoving and masses :math:`{\rm M}_\odot/h`.

***************
Getting Started
***************

.. toctree::
   :maxdepth: 1

   installation
   density_profiles

