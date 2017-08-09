from __future__ import print_function
import sys, os, glob
import setuptools
from setuptools import setup, Extension

os.system('ln -s ../include clusterwl/include')

sources = glob.glob(os.path.join('src','*.c'))
print('sources = ',sources)
headers = glob.glob(os.path.join('include','*.h'))
print('headers = ',headers)
ext=Extension("clusterwl._clusterwl", sources, depends=headers, include_dirs=['include'], extra_compile_args=[os.path.expandvars("-I${GSLI}")], extra_link_args=[os.path.expandvars("-L${GSLL}"),"-lgslcblas","-lgsl"])

dist = setup(name="clusterwl",
             author="Tom McClintock",
             author_email="tmcclintock89@gmail.com",
             description="Cluster weak lensing code.",
             license="MIT License",
             url="https://github.com/tmcclintock/Cluster_WL",
             packages=['clusterwl'],
             package_data={'clusterwl' : headers },
             ext_modules=[ext])

#setup.py doesn't put the .so file in the clusterwl directory, 
#so this bit makes it possible to
#import clusterwl from the root directory.  
#Not really advisable, but everyone does it at some
#point, so might as well facilitate it.
build_lib = glob.glob(os.path.join('build','*','clusterwl','_clusterwl*.so'))
if len(build_lib) >= 1:
    lib = os.path.join('clusterwl','_clusterwl.so')
    if os.path.lexists(lib): os.unlink(lib)
    os.link(build_lib[0], lib)
