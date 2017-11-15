******************************
Cluster Toolkit Documentation
******************************

Cluster Toolkit is a Python package specifically built for calculating weak lensing signals from galaxy clusters and cluster cosmology.
It consists of a Python front end wrapped around a well optimized backend in C, merged with cffi.
The core functionality of the package includes:

 * 3D density functions :math:`\rho(r)`
 * 3D correlation functions :math:`\xi(r)`
 * Halo bias models :math:`b(M)`
 * Projected density and differential profiles :math:`\Sigma(R)` and :math:`\Delta\Sigma`
 * Radially averaged profiles :math:`\overline{\Delta\Sigma}`
 * Boost factor models :math:`\mathcal{B} = (1-f_{\rm cl})^{-1}`
 * Miscentering affects on projected profiles :math:`f_{\rm mis}`
 * Halo mass functions :math:`\frac{dn}{dM}(M,z)` (in development)
 * Sunyaev-Zel'dovich (SZ) cluster signals :math:`Y_{SZ}` (in development)
 * Cluster magnification :math:`\kappa(\theta)` and shear profiles :math:`\gamma(\theta)` (in development)

The source code is publically available at https://github.com/tmcclintock/cluster_toolkit.

.. note::
   Unless stated otherwise, all distances are assumed to be :math:`{\rm Mpc}/h` comoving and masses :math:`{\rm M}_\odot/h`. Furthermore, power spectra :math:`P(k)` must be in units of :math:`({\rm Mpc}/h)^3` with wavenumber :math:`k` in units of :math:`h/{\rm Mpc}`.

***************
Getting Started
***************

.. toctree::
   :maxdepth: 1

   installation
   density_profiles
   correlation_functions
   halo_bias
   projected_density_profiles
   averaged_profiles
