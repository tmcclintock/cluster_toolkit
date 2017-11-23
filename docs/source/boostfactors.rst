************************************************************************
Boost Factors
************************************************************************

In galaxy cluster weak lensing, a significant systematic issue is cluster member dilution also known as boost factors. The idea is that if some cluster galaxies are misidentified as source (background) galaxies, then your weak lensing signal is diluted due to the fact that the cluster member galaxy won't be sheared or magnified. Traditionally, one calculates or estimates this correction and "boosts" the data vector. This boost factor is radially dependent, since you will tend to misidentify cluster members close to the cluster more than those farther out. Mathematically this looks like

.. math::

   \Delta\Sigma_{corrected}(R) = (1-f_{\rm cl})^{-1}(R)\Delta\Sigma(R)

where :math:`f_{\rm cl}` is the fraction of cluster members misidentified as being source galaxies. For shorthand, we write :math:`\mathcal{B} = (1-f_{\rm cl})^{-1}`. This module provides multiple models for :math:`\mathcal{B}`.

NFW Boost Model
==================

In McClintock et al. (in prep.) we model the boost factor with an NFW model:

.. math::

   \mathcal{B}(R) = 1+B_0\frac{1-F(x)}{x^2-1}

where :math:`x=R/R_s` and

.. math::

   F(x) = \Biggl \lbrace
   \begin{eqnarray}
   \frac{\tan^{-1}\sqrt{x^2-1}}{\sqrt{x^2-1}} : x > 1\\
   1 : x = 1\\
   \frac{\tanh^{-1}\sqrt{1-x^2}}{\sqrt{1-x^2}} : x < 1.
   \end{eqnarray}
   
Parameters that need to be specified by the user are :math:`B_0` and the scale radius :math:`R_s`. To use this, you would do:

.. code::

   from cluster_toolkit import boostfactors
   import numpy as np
   R = np.logspace(-2, 3, 100) #Mpc/h comoving
   B0 = 0.1 #Typical value
   Rs = 1.0 #Mpc/h comoving; typical value
   B = boostfactors.boost_nfw_at_R(R, B0, Rs)

Powerlaw Boost Model
=========================

In `Melchior et al. <https://arxiv.org/abs/1610.06890>`_ we used a power law for the boost factor.

.. math::

   \mathcal{B} = 1 + B_0\left(\frac{R}{R_s}\right)^\alpha

Here, the input parameters are :math:`B_0`, the scale radius :math:`R_s`, and the exponent :math:`\alpha`. This is also available in this module:

.. code::

   from cluster_toolkit import boostfactors
   import numpy as np
   R = np.logspace(-2, 3, 100) #Mpc/h comoving
   B0 = 0.1 #Typical value
   Rs = 1.0 #Mpc/h comoving; typical value
   alpha = -1.0 #arbitrary
   B = boostfactors.boost_powerlaw_at_R(R, B0, Rs, alpha)
