************************************************************************
Boost Factors
************************************************************************

In galaxy cluster weak lensing, a significant systematic issue is cluster member dilution also known as boost factors. The idea is that if some cluster galaxies are misidentified as source (background) galaxies, then your weak lensing signal is diluted due to the fact that the cluster member galaxy won't be sheared or magnified. Traditionally, one calculates or estimates this correction and "boosts" the data vector. This boost factor is radially dependent, since you will tend to misidentify cluster members close to the cluster more than those farther out. Mathematically this looks like

.. math::

   \Delta\Sigma_{corrected}(R) = (1-f_{\rm cl})^{-1}(R)\Delta\Sigma(R)

where :math:`f_{\rm cl}` is the fraction of cluster members misidentified as being source galaxies. For shorthand, we write :math:`\mathcal{B} = (1-f_{\rm cl})^{-1}`. This module provides multiple models for :math:`\mathcal{B}`.

NFW Boost Model
==================

In McClintock et al. (in prep.) we model the boost factor with an NFW model.

in progress

