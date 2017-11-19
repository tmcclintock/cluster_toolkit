************************
Installation
************************

To install the cluster_toolkit you currently need to build it from source::
  
  git clone https://github.com/tmcclintock/cluster_toolkit.git
  cd cluster_toolkit
  python setup.py install

Requirements
============
This package has only ever been tested with Python 2.7.x and has some dependencies. The Python dependencies that you can get with pip are:
  
- `Numpy <http://www.numpy.org/>`_: 1.13 or later

- `cffi <https://cffi.readthedocs.io/en/latest/>`_: 1.10 or later

- `pytest <https://docs.pytest.org/en/latest/>`_: 3.x or later for testing

In addition, you must have the `GNU Science Library <https://www.gnu.org/software/gsl/>`_ (GSL) installed. If you follow the instructions in their INSTALL file you will be done in only a few lines in the terminal. There is a pip installable GSL, but I do not know if it will work with the cluster toolkit.

Furthermore, while it is not a dependency of this code, for some of the functions you will need a way to calculate the linear and nonlinear matter power spectra :math:`P(k)`. Two good options are `CAMB <http://camb.info/>`_ and `CLASS <http://class-code.net/>`_. Both are also available in the `Core Cosmology Library <https://github.com/LSSTDESC/CCL>`_.
