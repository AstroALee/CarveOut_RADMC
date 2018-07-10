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
is_periodic = 0

# Number of cells on base level
# Only used if setup_type = 0
Ncells = [32,32,32]
Ncells = [64,64,64]

# Units of the below box values ('pc','cm','AU','ly' accepted)
# Follows Orion2 convenions where 0,0,0 is the center of the box.
# Depending on the setup_type value, the 'right' values may be overwritten.
box_units = 'pc'
#box_L = [1.5,1.5,1.5]
#box_L = [-0.5,-0.5,-0.5]
box_L = [1.61218e18, -1.10413e18, 2.44835e18] #32
box_L = [0.428724, -0.451574, 0.699705] #64
#box_L = [1.212e18, -1.104e18, 1.848e18]   # issue: 32x32, level 1 boxes have indices too large
box_R = [0.75,0,0.5]

# AMR levels included (0 = unigrid)
max_level = 2

# Max Temperature allowed
max_temp = 299.0

# Information about the species XYZ you want to analyze (e.g., CO, NH3)
x_XYZ   = 1.0e-7 # number fraction of species XYZ
x_H2    = 1.0    # number fraction of molecular hydrogen = is most of the mass

# adiabatic index and mean molecular weight for gas (used to calc temperature)
gas_gamma = 1.0001
gas_mmw   = 2.33
dustgas_ratio = 1.0e-2

# HDF5 and sink particle file names
O2gas_fname = "data.0736.3d.hdf5"
O2sink_fname = "data.0736.3d.sink"

# Output file names for use in RADMC3D
out_afname = "amr_grid.inp"       # output file name for amr grid
out_nfname = "numberdens_nh3.inp" # output file name for species
out_h2fname = "numberdens_h2.inp" # output file name for H_2
out_vfname = "gas_velocity.inp"   # output file name for velocity
out_tfname = "gas_temperature.inp"    # output file name for temperature
out_ddfname = "dust_density.inp" # output file name for dust density
out_dtfname = "dust_temperature.dat" # output for dust temperature (requires .dat)

# Names for various fields in YT. Unlikely you need to change these,
# unless Orion2 changes names of outputs.
density_field   = "density" # name of density field when accessing HDF5 file in YT
velocityX_field = "velocity_x" # name of x-velocity field in YT
velocityY_field = "velocity_y" # name of y-velocity field in YT
velocityZ_field = "velocity_z" # name of z-velocity field in YT
ienergy_field  = "energy-density"
