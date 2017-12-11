************************************************************************
Using CAMB
************************************************************************

`CAMB <http://camb.readthedocs.io/en/latest/>`_ is similar to CLASS in that it is a Boltzmann code used to compute the matter power spectrum. Very recently, a Python wrapper was created for CAMB, however it is less documented and less modular compared to CLASS. However, for the sake of comparison you can follow this script to use CAMB to calculate the linear and nonlinear matter power spectra.

.. note::
   CAMB and CLASS have differences that can cause >1% level changes in things like the mass function and possibly lensing. In general, pick one and make it explicit that you are using it when you describe your work.

.. note::
   CAMB outputs tend to have different shapes than you would expect.

.. code-block:: python

   import camb
   from camb import model, initialpower

   #Set cosmological parameters
   pars = camb.CAMBparams()
   pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
   pars.set_dark_energy(w=-1.0)
   pars.InitPower.set_params(ns=0.965)

   #This sets the k limits and specifies redshifts
   pars.set_matter_power(redshifts=[0., 0.8], kmax=2.0)

   #Linear P(k)
   pars.NonLinear = model.NonLinear_none
   results = camb.get_results(pars)
   kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 1000)

   #Note: the above function has the maxkh argument for specifying a different
   #kmax than was used above.
   #Note: pk has the shape (N_z, N_k)
   
   #Non-Linear spectra (Halofit)
   pars.NonLinear = model.NonLinear_both
   results.calc_power_spectra(pars)
   khnl, znl, pknl = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 1000)

