#!/usr/bin/env python
"""
Python wrapper to the fortran hermite polynomials library 

Module '_hermite' is auto-generated with f2py (version:2).
"""
from __version__ import __version__

import _PCE

##########################################################################################################
#
# functions
#
##########################################################################################################

def hermPol(p, x):
    return _PCE.hermite.hermpol(p, x)

def hermPol2(p, x):
    return _PCE.hermite.hermpol2(p, x)

def pce(p, seq, x):
    return _PCE.m_pce.pce(p, seq, x)

