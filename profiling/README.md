Profiling & Debugging
=====================

This directory is for routines used to profile and debug the C code without using the Python code at all.

If you are a user, don't bother looking at this directory. It's just used for optimization.

If you are a developer you can compile everything here with `make` and then run the code. To run it by itself you run `./profile`, and to run valgrind on it you use `valgrind --tool=memcheck ./profile`. This requires you to have [valgrind](http://valgrind.org/) installed.