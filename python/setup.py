from distutils.core import setup
from Cython.Build import cythonize

setup(name="mc_simulation", ext_modules=cythonize('MC_simulation.pyx'),)
