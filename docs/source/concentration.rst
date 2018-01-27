******************************
Mass-concentration Relations
******************************

The inner regions of clusters are called the 1-halo regime. The 1-halo regime is often modeled analytically using either an NFW or Einasto profile. Both of these profiles are functions of at least two variables. Numerous papers have examined the relationship between mass and NFW halo concentration. This module implements the Diemer-Kravtzov 2015 M-c relation. At present, only :math:`M_{200c}` and :math:`M_{200m}` mass definitions are supported. In the near future this will be expanded to :math:`\Delta\neq200` as well as :math:`M_{\rm vir}`.

To call this function you would do the following:

.. code::

   from cluster toolkit import concentration as conc
   M = 1e14 #Msun/h
   Omega_m = 0.3 #Matter fraction
   Omega_b = 0.05 #Baryon fraction
   ns = 0.96 #Power spectrum index
   h = 0.7 #Hubble constant
   #Assume that k and P come from somewhere, e.g. CAMB or CLASS
   #k are wavenumbers in h/Mpc and P is the linear power spectrum
   #in (Mpc/h)^3
   #The Mass_type argument can either be 'mean' or 'crit'
   Mass_type = 'mean'
   c = conc.concentration_at_M(M, k, P, ns, Omega_b, Omega_m, h, Mass_type)
