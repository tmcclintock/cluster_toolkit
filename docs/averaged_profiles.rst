************************************************************************
Radially Averaged Projected Profiles
************************************************************************

Weak lensing measurements are bin-averaged quantities. That is, they are measurements of a quantity with a radial bin around the lens. This module allows for calculating radially averaged projected profiles from continuous profiles. Mathematically this is

.. math::

   \overline{\Delta\Sigma} = \frac{2}{R_2^2-R_1^2}\int_{R_1}^{R_2}{\rm d}R' R'\Delta\Sigma(R').

This can be computed in the code by using

.. code::

   from cluster_toolkit import averaging
   #Assume DeltaSigma at R_perp are computed here
   N_bins = 15
   bin_edges = np.logspace(np.log10(0.2), np.log10(30.), N_bins+1)
   #Bin edges are from 200 kpc/h to 30 Mpc/h
   averaged_DeltaSigma = np.zeros(N_bins)
   averaging.average_profile_in_bins(bin_edges, R_perp, DeltaSigma, averaged_DeltaSigma)

.. note::

   The :code:`average_profile_in_bins` function can work with any projected profile.

.. note::

   For now, :code:`averaged_profile` must passed in before calling this function. This will be updated soon so that this isn't needed.
