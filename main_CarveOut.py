'''
    asldkjfadsf
'''

import yt
from globals_CarveOut import *
import inputs_CarveOut as inputs
from carver_CarveOut import CarvingWriter


# Definition of the XYZ species field. Uses info from inputs
def _NumberDensityXYZ(field, data):
    aw = yt.YTQuantity(inputs.aw_XYZ*Csts.mp, 'g')
    return (inputs.x_XYZ/aw)*data[inputs.density_field]
yt.add_field(("gas", "number_density_XYZ"), function=_NumberDensityXYZ, units="cm**-3")


# -=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=-
#                           PROGRAM STARTS HERE
# -=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=-
# Greetings
Messages.Welcome()

# Gets box lengths ready (may overwrite the right-side ones)
boxLeft  = [inputs.box_xL,inputs.box_yL,inputs.box_zL]
boxRight = [inputs.box_xR,inputs.box_yR,inputs.box_zR]
boxLeft =  [ Convert(x,inputs.box_units,'cm','cm') for x in boxLeft] # convert to CGS
boxRight = [ Convert(x,inputs.box_units,'cm','cm') for x in boxRight] # convert to CGS

if not any([inputs.is_periodic==y for y in [0,1]]):
    Messages.Calamity("is_periodic needs to be 0 or 1. C'mon!")
if not any([inputs.setup_type==y for y in [0,1]]):
    Messages.Calamity("setup_type needs to be 0 or 1. C'mon!")
assert(inputs.aw_XYZ>0)
assert(inputs.x_XYZ>0 and inputs.x_XYZ<1.0)

# Loads file into YT
ds = yt.load(inputs.O2gas_fname)
try:
    print("Loaded file " + str(ds))
except NameError:
    Messages.Calamity("Unable to properly load file into YT!")

# Cell size on coarsest level (will assume factors of 2 refinement)
dxArray = ds.index.dds_list[0]

# If setup_type=0, we'll reset the right side of the box
# else we'll use the box dimensions to set the cell dimensions
if(inputs.is_periodic):
    inputs.setup_type = 0

if(inputs.setup_type==0):
    for i in range(3):
        boxRight[i] = boxLeft[i] + (inputs.Ncells[i])*dxArray[i]
else:
    for i in range(3):
        if any([(y-x)<0 for y,x in zip(boxRight,boxLeft)]):
            Messages.Calamity("Box physical dimensions are weird")
        inputs.Ncells[i] = int(np.ceil((boxRight[i] - boxLeft[i])/dxArray[i]))

# Checks & Balances
assert(inputs.max_level <= ds.max_level)
assert(all(np.array(inputs.Ncells) >= 2))
if(inputs.is_periodic):
    if(any( np.array(boxRight)-np.array(boxLeft) >  np.array(ds.domain_right_edge.d)-np.array(ds.domain_left_edge.d))):
        Messages.Calamity("Carve out box is bigger than domain. Weird!")
else:
    if(any(boxLeft < ds.domain_left_edge.d)):
        Messages.Calamity("Left side of desired carve out outside actual domain. " + str(boxLeft) + " , " + str(ds.domain_left_edge.d))
    if(any(boxRight > ds.domain_right_edge.d)):
        Messages.Calamity("Right side of desired carve out outside actual domain. " + str(boxRight) + " , " + str(ds.domain_right_edge.d))


print("Carving between Left = " + str(boxLeft))
print("            to Right = " + str(boxRight))
print("           w/ Ncells = " + str(inputs.Ncells))

# Call the writer class
writer = CarvingWriter(ds, boxLeft, boxRight, inputs.Ncells, inputs.max_level, inputs.is_periodic)

# Write the amr grid file
print("Writing amr grid file")
writer.write_amr_grid()

# Write the number density file for species or dust
print("Writing number density file")
writer.write_line_file(("gas", "number_density_XYZ"), inputs.out_nfname)

# Write the gas velocity file
print("Writing velocity file")
velocity_fields = [inputs.velocityX_field,inputs.velocityY_field,inputs.velocityZ_field]
writer.write_line_file(velocity_fields, inputs.out_vfname)

# We're done!
Messages.Goodbye()
