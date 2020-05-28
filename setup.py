#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 11:52:50 2020

@author: filipe
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize('util.pyx'))