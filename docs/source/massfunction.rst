******************************
Halo Mass Function
******************************

Clusters and galaxies live inside of dark matter halos, which represent peaks in the matter density field. This means that it is possible to use the abundance of clusters to constrain cosmology if you have well understood mapping from clusters onto halos (which is outside the scope of this code).

The abundance of halos is also known as the *halo mass function*. A general mathematical form for the mass function is given in `Tinker et al. (2008) <https://arxiv.org/abs/0803.2706>`_ (where they cite `Press & Schechter (1974) <http://adsabs.harvard.edu/abs/1974ApJ...187..425P>`_, `Jenkins et al. (2000) <https://arxiv.org/abs/astro-ph/0005260>`_ and other papers that you can look up) is:

.. math::

   \frac{{\rm d}n}{{\rm d}M} = f(\sigma)\frac{\rho_m}{M}\frac{{\rm d}\ln\sigma^{-1}}{{\rm d}M}.

Where :math:`\sigma` is the RMS variance of a spherical top hat containing a mass :math:`M`, :math:`\rho_m=\rho_{\rm crit}\Omega_m` is the mean matter density and :math:`f(\sigma)` is known as the *halo multiplicity function*. Practically speaking, what sets one mass function model apart from another is how the multiplicity is written down.

At some point in the future the toolkit will have other options for the mass function, but for now it implements the mass function from `Tinker et al. (2008) <https://arxiv.org/abs/0803.2706>`_. Specifically, the version in Appendix C of that paper, which is usually the one people mean when they refer to this paper.

.. note::

   Implicitely, the mass function depends on the linear matter power spectrum :math:`P(k)` in order to map from :math:`M` to :math:`\sigma`. Since the toolkit doesn't have its own power spectrum implemented, the user must input one from, e.g. CLASS or CAMB.

Tinker Mass Function
====================

The Tinker mass function is defined by it's multiplicity function, which looks like

.. math::

   f(\sigma) = B\left[\left(\frac{\sigma}{e}\right)^{-d} + \sigma^{-f}\right]\exp(-g/\sigma^2).

In this mass function :math:`d`, :math:`e`, :math:`f`, and :math:`g` are free parameters that depend on halo definition. By default they are set to the values associated with :math:`M_{200m}`. If you want to switch to other values for other halo definitions, see `Tinker et al. (2008) <https://arxiv.org/abs/0803.2706>`_. The normalization :math:`B` is calculated internally so that you don't have to pass it in.

To use this in the code you would do:

.. code::

   from cluster_toolkit import massfunction
   import numpy as np
   #Assume that k and P come from somewhere, e.g. CAMB or CLASS
   #Units of k and P are h/Mpc and (Mpc/h)^3
   Mass = 1e14 #Msun/h
   Omega_m = 0.3 #example value
   dndM = massfunction.dndM_at_M(Mass, k, P, Omega_m)
   #Or could also use an array
   Masses = np.logspace(12, 16)
   dndM = massfunction.dndM_at_M(Masses, k, P, Omega_m)


Binned Mass Functions
=====================

In reality in a simulation or in cluster abundance the real observable is the number density of objects in some mass bin of finite width. Written out, this is

.. math::

   n = \int_{M_1}^{M_2}{\rm d}M\ \frac{{\rm d}n}{{\rm d}M}.

In the toolkit this is available by first calculating :math:`{\rm d}n/{\rm d}M` and then passing that back to the toolkit. This is available in the code by using

.. code::

   from cluster_toolkit import massfunction
   import numpy as np
   #Assume that k and P come from somewhere, e.g. CAMB or CLASS
   #Units of k and P are h/Mpc and (Mpc/h)^3
   Omega_m = 0.3 #example value
   M = np.logspace(12, 16) #Msun/h
   dndM = massfunction.dndM_at_M(M, k, P, Omega_m)
   M1 = 1e12 #Msun/h
   M2 = 1e13 #Msun/h
   n = n_in_bin(M1, M2, M, dndM)

You can also pass in many bin edges at once:

.. code::

   edges = np.array([1e12, 5e12, 1e13, 5e13])
   n = n_in_bins(edges, M, dndM)
