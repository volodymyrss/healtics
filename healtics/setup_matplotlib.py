import os
import numpy as np # type: ignore
from matplotlib import rcParams, rc # type: ignore

# common setup for matplotlib
params = {'backend': 'pdf',
          'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 10,
          #'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 10,
          'text.usetex': os.environ.get("MATPLOTLIB_USETEX", 'no') == 'yes',
          'font.family':'sans-serif',
         # free font similar to Helvetica
          'font.sans-serif':'FreeSans'}

# use of Sans Serif also in math mode
if params['text.usetex']:
    rc('text.latex', preamble=r'\usepackage{sfmath}')

rcParams.update(params)

def cm2inch(cm):
    """Centimeters to inches"""
    return cm * 0.393701

