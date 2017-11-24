************************************************************************
Using CLASS
************************************************************************

[CLASS](http://class-code.net/) is a code used to compute the matter power spectrum. The power spectrum is a key input for the cluster-toolkit. This is meant to be a very short example of how you can call CLASS to get the linear and nonlinear power spectra.

.. note::
   The CLASS github page is [here](https://github.com/lesgourg/class_public). The CLASS documentation is [here](https://github.com/lesgourg/class_public/blob/master/explanatory.ini).

.. note::
   CLASS uses units of :math:`Mpc^{-1}` for math:`k` and :math:`Mpc^3` for math:`P`.

.. code-block:: python

   from classy import Class
   import numpy as np
   
   #Start by specifying the cosmology
   Omega_b = 0.05
   Omega_m = 0.3
   Omega_cdm = Omega_m - Omega_b
   h = 0.7 #H0/100
   A_s = 2.1e-9
   n_s = 0.96

   #Create a params dictionary
   #Need to specify the max wavenumber
   k_max = 10 #UNITS: 1/Mpc

   params = {
		'output':'mPk',
		'non linear':'halofit',
		'Omega_b':Omega_b,
		'Omega_cdm':Omega_cdm
		'h':h,
		'A_s':A_s,
		'n_s':n_s,
		'P_k_max_1/Mpc':k_max,
		'z_max_pk':10. #Default value is 10
   }

   #Initialize the cosmology andcompute everything
   cosmo = Class()
   cosmo.set(params)
   cosmo.compute()

   #Specify k and z
   k = np.logspace(-5, np.log10(k_max), num=1000) #Mpc^-1
   z = 1.

   #Call these for the nonlinear and linear matter power spectra
   Pnonlin = np.array([cosmo.pk(ki, z) for ki in k])
   Plin = np.array([cosmo.pk_lin(ki, z) for ki in k])
