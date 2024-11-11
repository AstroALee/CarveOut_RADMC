"""
   main_CV.py

   Purpose:
        Driver for RadMC carve routines. Call this file when you want
        to run the code. Should not ever need to edit this file.

   Author:
        Aaron T. Lee, atl8@stmarys-ca.edu
        Spring 2020

   Written/Tested with Python 3.7.3, YT 3.4.1
"""

# python modules
import numpy as np


# inputs (also read in inside globals?)
from inputs_CV import inputs

# personal modules
from globals_CV import *
from yfields_CV import *

# Carving routine
from carver_CV import CarvingWriter


if __name__ == "__main__":
    # Welcome to the code!
    #print("The value of __name__ is:", repr(__name__))
    Messages.Welcome()

# Loads file into YT
ds = yt.load(inputs.O2gas_fname)
try:
    Messages.Print(inputs.verb,0,"Loaded file " + str(ds)) # using classic print() for try/except to work
except NameError:
    assert False, "Unable to properly load file into YT!"

# Unclear why you would want to carve out at a resolution higher than what is
#
assert(inputs.max_level <= ds.max_level)

# Find the cell dimensions that work for the desired domain ?
# Convert to CGS, if needbe
box_Lcgs = [ Convert(x,inputs.box_units,'cm','cm') for x in inputs.box_L ]
# Find what cell numbers will place us just shy of this boundary
box_Lcell = [ int( ( box_Lcgs[x] - ds.domain_left_edge[x].d  )/ds.index.dds_list[0][x] - 1 ) for x in range(3) ]


writer = CarvingWriter(ds, box_Lcell, ds.domain_left_edge, ds.domain_dimensions, ds.index.dds_list[0])
