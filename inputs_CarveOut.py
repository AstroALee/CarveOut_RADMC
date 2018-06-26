"""
   inputs_CarveOut.py

   Purpose:
        Input file for RadMC carving routines. This is the only file you should
        edit. Follow the comments below to see what each variable represents.

   Author:
        Aaron T. Lee, aaron.t.lee@utexas.edu
        Spring 2018

   Written/Tested with Python 3.6.4, YT 3.4.1
"""

# How are we setting up the domain?
#    0 = will read the left-side locations and use the number of cells on base
#        level to determine right-side locations
#    1 = will ignore Nx, Ny, Nz and use the left and right box sides to determine
#        box. Will round to nearest cell on base level.
setup_type = 0

# Allow for periodic boundary conditions (setting to 1 overwrites setup_type to 0)?
is_periodic = 1

# Number of cells on base level
# Only used if setup_type = 0
Ncells = [8,8,8]

# Units of the below box values ('pc','cm','AU','ly' accepted)
# Follows Orion2 convenions where 0,0,0 is the center of the box.
# Depending on the setup_type value, the 'right' values may be overwritten.
box_units = 'pc'
box_L = [0.5,-0.5,0]
box_R = [1.5,0,0.5]

# AMR levels included (0 = unigrid)
max_level = 0

# Information about the species XYZ you want to analyze (e.g., CO, NH3)
x_XYZ   = 1.0e-4 # mass fraction of species XYZ
aw_XYZ  = 1.4    # atomic weight of species XYZ per H atom (1/mu = x / aw*mp)

# HDF5 and sink particle file names
O2gas_fname = "data.0622.3d.hdf5"
O2sink_fname = "data.0622.3d.sink"

# Output file names for use in RADMC3D
out_afname = "amr_grid.inp"       # output file name for amr grid
out_nfname = "numberdens_xyz.inp" # output file name for species
out_vfname = "gas_velocity.inp"   # output file name for velocity

# Names for various fields in YT. Unlikely you need to change these,
# unless Orion2 changes names of outputs.
density_field   = "density" # name of density field when accessing HDF5 file in YT
velocityX_field = "velocity_x" # name of x-velocity field in YT
velocityY_field = "velocity_y" # name of y-velocity field in YT
velocityZ_field = "velocity_z" # name of z-velocity field in YT
