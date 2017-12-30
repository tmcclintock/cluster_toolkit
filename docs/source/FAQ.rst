************************************************************************
Frequently Asked Questions
************************************************************************

I'm getting an interpolation error
----------------------------------------

The backend of the toolkit is written in C. It is likely that you are passing in data that in Python is actually of floating point precision (32 bits), when everything in C is expected to be in double precision (64 bits). You must do a cast in Python to correct this.

I'm getting random GSL errors
--------------------------------

Python does not order arrays in memory the same way as in C. This can cause strange behavior as seemingly random memory addresses get used, and sometimes even segmentation faults. It is likely that your array in Python is not "C-ordered". Use `this numpy function <https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.ascontiguousarray.html>`_ to force your input arrays to have the correct ordering in memory.

These two points are outstanding issues on the `github issues page <https://github.com/tmcclintock/cluster_toolkit/issues>`_. One day they will be taken care of automatically without causing the code to slow significantly.
