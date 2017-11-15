******************************
Correlation Functions
******************************

Cluster density profiles are closely related to the correlation function

.. math::
   
   \langle\rho(r)\rangle = \rho_0(1+\xi(r))

This module makes various correlation functions available. The correlation functions for the NFW and Einasto profiles can also be computed directly from the density profiles by inverting the above equation.

NFW Profile
===========

The NFW profile (`arxiv <https://arxiv.org/abs/astro-ph/9508025>`_) is a 3D correaltion function given by:

.. math::

   \xi_{nfw}(r) = \frac{\rho_{nfw}}{\Omega_m\rho_{crit}} - 1

The free parameters are the cluster mass :math:`M_\Delta` and concentration :math:`c_\Delta = r_\Delta/r_s`. In this module we choose to define the density with respect to the matter background density :math:`\Omega_m\rho_{crit}`. The scale radius :math:`r_s` is given in :math:`{\rm Mpc}/h`, however the code uses  as an argument instead. As written, because of the choice of units the only cosmological parameter that needs to be passed in is :math:`\Omega_m`. The arguments are identical to the density profile.

.. note::
   The correlation functions can use :math:`\Delta\neq 200`.

To use this, you would do:

.. code::

   from cluster_toolkit import xi
   import numpy as np
   radii = np.logspace(-2, 3, 100) #Mpc/h comoving
   mass = 1e14 #Msun/h
   concentration = 5 #arbitrary
   Omega_m = 0.3
   xi_nfw = xi.xi_nfw_at_R(radii, mass, concentration, Omega_m)


Einasto Profile
===============

The Einasto profile (`ADS <http://adsabs.harvard.edu/abs/1965TrAlm...5...87E>`_) is a 3D density profile given by:

.. math::

   \xi_{ein}(r) = \frac{\rho_{ein}}{\Omega_m\rho_{crit}} - 1
   

In this model, the free parameters are the scale radius :math:`r_s`, :math:`\alpha`, and the cluster mass :math:`M_\Delta`. The scale density :math:`\rho_s` is calculated internally, or can be passed in instead of mass. This is the same arguments as the density profile. To use this, you would do:

.. code::

   from cluster_toolkit import xi
   import numpy as np
   radii = np.logspace(-2, 3, 100) #Mpc/h comoving
   mass = 1e14 #Msun/h
   r_scale = 1.0 #Mpc/h comoving scale radius
   alpha = 0.19 #arbitrary; a typical value
   Omega_m = 0.3
   xi_ein = xi.xi_einasto_at_R(radii, mass, r_scale, alpha , Omega_m)

Matter-matter Correlation Function
=============================================

in progress

Halo-matter Correlation Function
=============================================

in progress

1-halo + 2-halo Terms
=============================================

in progress
