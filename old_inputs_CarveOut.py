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

# How verbose is the output to the screen?
#    0 = bare min is printed,
#    1 = some printed, >=2 = all
verbosity = 2

# Allow for periodic boundary conditions (setting to 1 overwrites setup_type to 0)?
#    In the future this will always be set to 1.
is_periodic = 0

# Number of cells on base level of the carved out region
Ncells = [32,32,32]

# Units of the below box values ('pc','cm','AU','ly' accepted)
# Follows Orion2 convenions where 0,0,0 is the center of the box.
# Set box_L to the coordinates for the lower-left corner of the carve-out domain
# The values should match the unit given by box_units
box_units = 'pc'
box_L = [0.518976, 0.480555, 0.804028]

# AMR levels included (0 = unigrid at base level from simulation output)
max_level = 2

# Do you want output unigrid at max_level?
#   0 = AMR structure is mapped if max_level > 0
#   1 = unigrid output at max_level when max_level > 0
is_unigrid = 1

# force nested?
force_nested = 0

# Information about the species XYZ you want to analyze (e.g., CO, NH3)
x_XYZ   = 1.0e-7 # number fraction of species XYZ
x_H2    = 1.0    # number fraction of molecular hydrogen = is most of the mass

# Max Temperature allowed when calulating the local gas temperature. Sometimes
# sink particle regions create spurious high temps.
max_temp = 299.0

# Freezeout / Dissociation thresholds
allow_freezeout = 0 # if = 1, abundance of species XYZ is set to 0 if H_2 number density is beyond these limits
freeze_minN = 1e2 # min number density of H_2 necessary, must be >=0, set to 0 for no min cutoff
freeze_maxN = 1e5 # max number density of H_2 necessary, must be >=freeze_minN, set to 1e30 for no max cutoff


# Properties of the gas
# adiabatic index and mean molecular weight for gas (used to calc temperature)
gas_gamma = 1.0001  # Necessary this exactly matches what you used in orion2.ini; Extra 0's will throw off T calculation
gas_mmw   = 2.33    # 2.33 = contemp gas
dustgas_ratio = 1.0e-2  # dust to gas mass ratio

# Microturbulence included?
has_microturb = 0 # =1, creates a microturbulence file that employs the MicroTurb function in main_CarveOut. Default is to set the
                  # microturbulence speed to the same constant everywhere.
microturbulence_speed = 1e5 # cgs, please

# HDF5 and sink particle file names to read in
# (as of now, the sink file isn't actually used)
O2gas_fname = "data.0700.3d.hdf5"
O2sink_fname = "data.0700.3d.sink"

# Output file names for use in RADMC3D
out_afname = "amr_grid.inp"       # output file name for amr grid
out_nfname = "numberdens_nh3.inp" # output file name for species XYZ above
out_h2fname = "numberdens_h2.inp" # output file name for H_2
out_vfname = "gas_velocity.inp"   # output file name for velocity
out_tfname = "gas_temperature.inp"    # output file name for temperature
out_ddfname = "dust_density.inp" # output file name for dust density
out_dtfname = "dust_temperature.dat" # output for dust temperature (requires .dat)
out_mtfname = "microturbulence.inp" # output for microturbulence

# Names for various fields in YT. Unlikely you need to change these,
# unless Orion2 changes names of outputs.
density_field   = "density" # name of density field when accessing HDF5 file in YT
velocityX_field = "velocity_x" # name of x-velocity field in YT
velocityY_field = "velocity_y" # name of y-velocity field in YT
velocityZ_field = "velocity_z" # name of z-velocity field in YT
ienergy_field  = "energy-density" # name of internal energy field in YT

# At this point, call python main_CarveOut.py and call it a day!
