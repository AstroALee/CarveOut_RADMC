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

# How verbose is the output to the screen? 0 = bare min is printed, 1 = some printed, >=2 = all
verbosity = 1

# Allow for periodic boundary conditions (setting to 1 overwrites setup_type to 0)?
is_periodic = 1

# Number of cells on base level
Ncells = [128,128,128]

# Units of the below box values ('pc','cm','AU','ly' accepted)
# Follows Orion2 convenions where 0,0,0 is the center of the box.
# Set box_L to the coordinates for the lower-left corner of the carve-out domain
# The values should match the unit given by box_units
box_units = 'pc'
box_L = [1.61218e18, -1.10413e18, 2.44835e18]

# AMR levels included (0 = unigrid)
max_level = 1

# Information about the species XYZ you want to analyze (e.g., CO, NH3)
x_XYZ   = 1.0e-7 # number fraction of species XYZ
x_H2    = 1.0    # number fraction of molecular hydrogen = is most of the mass

# Max Temperature allowed when calulating the local gas temperature. Sometimes
# sink particle regions create spurious high temps.
max_temp = 299.0

# Properties of the gas
# adiabatic index and mean molecular weight for gas (used to calc temperature)
gas_gamma = 1.0001  # Necessary this exactly matches what you used in orion2.ini; Extra 0's will throw off T calculation
gas_mmw   = 2.33    # 2.33 = contemp gas
dustgas_ratio = 1.0e-2  # dust to gas mass ratio

# HDF5 and sink particle file names to read in
# (as of now, the sink file isn't actually used)
O2gas_fname = "data.0736.3d.hdf5"
O2sink_fname = "data.0736.3d.sink"

# Output file names for use in RADMC3D
out_afname = "amr_grid.inp"       # output file name for amr grid
out_nfname = "numberdens_nh3.inp" # output file name for species XYZ above
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
ienergy_field  = "energy-density" # name of internal energy field in YT

# At this point, call python main_CarveOut.py and call it a day! 
