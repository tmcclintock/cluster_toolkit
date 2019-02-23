************************************************************************
Miscentering Effects
************************************************************************

If galaxy cluster centers are not properly identified on the sky, then quantities measured in annuli around that center will not match theoretical models. This effect is detailed in `Johnston et al. (2007) <http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:astro-ph/0507467>`_ and `Yang et al. (2006) <https://arxiv.org/abs/astro-ph/0607552>`_.

To summarize, if a cluster center is incorrectly identified on the sky by a distance :math:`R_{\rm mis}` then the surface mass density becomes:

.. math::

   \Sigma_{\rm mis}^{\rm single\ cluster}(R, R_{\rm mis}) = \int_0^{2\pi} \frac{{\rm d}\theta}{2\pi}\ \Sigma\left(\sqrt{R^2+R_{\rm mis}^2 + 2RR_{\rm mis}\cos\theta}\right).

That is, the average surface mass density at distance :math:`R` away from the incorrect center is the average of the circle drawn around that incorrect center. To get the miscentered profiles of a *single cluster* you would use

.. code::
   
   from cluster_toolkit import miscentering
   mass = 1e14 #Msun/h
   conc = 5 #arbitrary
   Omega_m = 0.3
   #Calculate Rp and Sigma here, where Sigma is centered
   Rmis = 0.25 #Mpc/h; typical value
   Sigma_mis_single = miscentering.Sigma_mis_single_at_R(Rp, Rp, Sigma, mass, conc, Omega_m, Rmis)

As you can see :code:`Rp` is passed in twice. It is first used as the location at which to evaluate :code:`Sigma_mis` and then as the locations at which :code:`Sigma` is known. So if you wanted those two radial arrays can be different.

The :math:`\Delta\Sigma` profile is defined the usual way

.. math::

   \Delta\Sigma(R,R_{\rm mis}) = \bar{\Sigma}_{\rm mis}(<R,R_{\rm mis}) - \Sigma_{\rm mis}(R,R_{\rm mis})

which can be calculated using this module using

.. code::

   DeltaSigma_mis_single = miscentering.DeltaSigma_mis_at_R(Rp, Rp, Sigma_mis_single)

Stacked Miscentering
==============================

In a stack of clusters, the amount of miscentering will follow a distribution :math:`P(R'|R_{\rm mis})` given some characteristic miscentering length :math:`R_{\rm mis}`. That is, some clusters will be miscentered more than others. `Simet et al. (2017) <https://arxiv.org/abs/1603.06953>`_ for SDSS and `Melchior et al. (2017) <https://arxiv.org/abs/1610.06890>`_ assume a Raleigh distribution for the amount of miscentering :math:`R'`:

.. math::

   P(R'|R_{\rm mis}) = \frac{R'}{R^2_{\rm mis}}\exp[-R'^2/2R_{\rm mis}^2]\,.

In `McClintock et al. (2019) <http://adsabs.harvard.edu/abs/2019MNRAS.482.1352M>`_ we used a Gamma profile for the mistnering:

.. math::

   P(R'|R_{\rm mis}) =  \frac{R'}{R^2_{\rm mis}}\exp[-R'/R_{\rm mis}]\,.

Both of these are available in the toolkit. We see that :math:`R_{\rm mis}` is a free parameter, giving rise to a miscentered projected stacked density profile:

.. math::

   \Sigma_{\rm mis}^{\rm stack}(R) = \int_0^\infty{\rm d}R'\ P(R'|R_{\rm mis})\Sigma_{\rm mis}^{\rm single\ cluster}(R, R')

which can then itself be integrated to get :math:`\Delta\Sigma_{\rm mis^{\rm stack}}`. To calculate these in the code you would use:

.. code::

   from cluster_toolkit import miscentering
   #Assume Sigma at R_perp are computed here
   Sigma_mis = miscentering.Sigma_mis_at_R(R_perp, R_perp, Sigma, mass, concentration, Omega_m, R_mis)
   DeltaSigma_mis = miscentering.DeltaSigma_mis_at_R(R_perp, R_perp, Sigma_mis)
   
