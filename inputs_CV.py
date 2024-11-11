"""
   inputs_CV.py

   Purpose:
        Input file for RadMC carving routines. This is the only file you should
        edit. Follow the comments below to see what each variable represents.

   Author:
        Aaron T. Lee, atl8@stmarys-ca.edu
        Spring 2020

   Written/Tested with Python 3.7.3, YT 3.4.1
"""

class inputs:

    ''' This information relates to the grid '''
    # How verbose is the output to the screen?
    #    0 = bare min is printed,
    #    1 = some printed, >=2 = all
    verb = 2

    # Do you want output unigrid at a desired level?
    #   0 = AMR structure is mapped if max_level > 0
    #   1 = unigrid output at max_level when max_level > 0
    #       unigrid output at base grid when max_level = 0
    is_unigrid = 0

    # If is_unigrid = 0, AMR levels included up to max_level
    #     (0 = unigrid at base level from simulation output)
    # If is_unigrid = 1, data is mapped to max_level resolution
    #     (set to 0 for mapping to base grid level)
    max_level = 1

    # Number of cells on BASE LEVEL of the carved out region
    Ncells = [128,128,128]

    # Units of the below box values ('pc','cm','AU','ly' accepted)
    # Follows Orion2 convenions where 0,0,0 is the center of the box.
    # Set box_L to the coordinates for the lowerleft corner of the carved domain
    # The values should match the unit given by box_units
    # The code will find the cell that is closest to this location.
    box_units = 'pc'
    box_L = [0.5,-0.5,0.5]

    # Allow for periodic boundary conditions.
    # In the future this will always be set to 1.
    # Currently, this is not employed.
    is_periodic = 0


    ''' This information relates to physics '''
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
    has_microturb = 0 # =1, creates a microturbulence file that employs the
                      # MicroTurb function in main_CarveOut. Default is to set
                      # the microturbulence speed to the same constant everywhere.
    microturbulence_speed = 1e5 # cgs, please


    ''' This information relates to filenames '''

    # HDF5 and sink particle file names to read in
    # (as of now, the sink file isn't actually used)
    O2gas_fname = "data.0497.3d.hdf5"
    O2sink_fname = "data.0497.3d.sink"

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




''' ALL DONE, below will just do some checks on the values you inputted, no need to change this!
            ------=------=------=------=------=------=------=------=------=------=------=       '''

# check some values from inputs
assert (inputs.x_XYZ>0 and inputs.x_XYZ<1.0), "Number fraction outside physical limits!"
assert (inputs.freeze_minN >= 0), "Really....?"
assert (inputs.freeze_maxN >= inputs.freeze_minN), "Really....?"
assert (any([inputs.is_periodic==y for y in [0,1]])), "is_periodic needs to be 0 or 1. C'mon!"
assert (any([inputs.is_unigrid==y for y in [0,1]])), "is_unigrid needs to be 0 or 1. C'mon!"
assert (any([inputs.allow_freezeout==y for y in [0,1]])), "allow_freezeout needs to be 0 or 1. C'mon!"
assert (any([inputs.has_microturb==y for y in [0,1]])), "has_microturb needs to be 0 or 1. C'mon!"
assert(inputs.verb>=0), "Verbosity needs to be at least 0"
