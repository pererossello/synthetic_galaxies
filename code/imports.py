import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sympy as sp
import PIL

# gala imports
import gala.coordinates as gc
import gala.dynamics as gd
import gala.potential as gp
from gala.units import galactic

# astropy imports
import astropy.units as u
from astropy.constants import G as G_val
G_val = G_val

import plot_utils as pu
import utils as ut
import gutils as gu

mpl.rcParams['text.usetex'] = True