************************************************************************
Miscentering Effects
************************************************************************

If galaxy cluster centers are not properly identified on the sky, then quantities measured in annuli around that center will not match theoretical models. This effect is detailed in `Johnston et al. (2007) <http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:astro-ph/0507467>`_ and `Yang et al. (2006) <https://arxiv.org/abs/astro-ph/0607552>`_.

To summarize, if a cluster center is incorrectly identified on the sky by a distance :math:`R_{\rm mis}` then the surface mass density becomes:

.. math::

   \Sigma_{\rm mis}(R| R_{\rm mis}) = \int_0^{2\pi} \frac{{\rm d}\theta}{2\pi}\ \Sigma\left(\sqrt{R^2+R_{\rm mis}^2 + 2RR_{\rm mis}\cos\theta}\right).

That is, the average surface mass density at distance :math:`R` away from the incorrect center is the average of the circle drawn around that incorrect center.
