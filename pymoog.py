#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


# Why a single leading underscore? Because this is an internal function. See
# pep8 for more information here:
# http://legacy.python.org/dev/peps/pep-0008/#naming-conventions
def _run_moog():
    os.system('MOOGSILENT > tmp')
    return 1
