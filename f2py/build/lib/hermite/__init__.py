#!/usr/bin/env python
"""
Python wrapper to the fortran hermite polynomials library 

Module '_hermite' is auto-generated with f2py (version:2).
"""
from __version__ import __version__

import _hermite

##########################################################################################################
#
# functions
#
##########################################################################################################

def hermPol(p, x):
    return _hermite.hermite.hermpol(p, x)

def hermPol2(p, x):
    return _hermite.hermite.hermpol2(p, x)

