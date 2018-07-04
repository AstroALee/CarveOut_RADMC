"""
   main_CarveOut.py

   Purpose:
        Driving function for RadMC carve routines. Call this file when you want
        to run the code. Should not ever need to edit this file.

   Author:
        Aaron T. Lee, aaron.t.lee@utexas.edu
        Spring 2018

   Written/Tested with Python 3.6.4, YT 3.4.1
"""

import yt
from globals_CarveOut import *
import inputs_CarveOut as inputs
from carver_CarveOut import CarvingWriter


# Definition of the gas temperature. Uses info from inputs
def _GasTemperature(field, data):
    ethbyden = data[inputs.ienergy_field]/data[inputs.density_field]
    gmmfac = (inputs.gas_gamma - 1.0)*yt.YTQuantity(inputs.gas_mmw*Csts.mp, 'g')/(yt.units.kb)
    return (gmmfac*ethbyden)
yt.add_field(("gas", "gas_temperature"), function=_GasTemperature, units="K")

# Definition of the XYZ species field. Uses info from inputs
def _NumberDensityXYZ(field, data):
    aw = yt.YTQuantity(inputs.aw_XYZ*Csts.mp, 'g')
    return (inputs.x_XYZ/aw)*data[inputs.density_field]
yt.add_field(("gas", "number_density_XYZ"), function=_NumberDensityXYZ, units="cm**-3")


def BreakIntoPatches(domainPatches,fulldomain, mincellsize):
    fdL = fulldomain[0]
    fdR = fulldomain[1]

    # As of right now, domainPatches is just the desired carved out region .
    # if for some unknown reason you inputted a left edge that's beyond the
    # actual domain's right edge, let's figure that out now
    leftAdjustment = np.array(len(domainPatches[0][0])*[0])
    for i in range(len(domainPatches[0][0])):
        if(domainPatches[0][0][i] > fulldomain[1][i]):
            leftAdjustment[i] = domainPatches[0][0][i] - fulldomain[1][i]
        else:
            leftAdjustment[i] = 0.
    print(leftAdjustment)

    while(True):
        print("Beginning of while loop, current domainPatches:")
        print(domainPatches)
        # look at current set of domainPatches to see if any exceed fulldomain
        # If all are within fulldomain, we're done. List comprehension was giving
        # me a headache, so we're doing this old-school.
        done = 1
        for boxL,boxR in domainPatches:
            if (np.array(boxR)>np.array(fdR)).any():
                done = 0
                break # no need to keep looping through...
        if(done):
            print("All patches inside domain")
            break
        else:
            # now we need to find the first example of where this is not true
            curIdx = 0
            for boxL,boxR in domainPatches:
                if (np.array(boxR)>np.array(fdR)).any():
                    curIdx = domainPatches.index([boxL,boxR])
                    break
            # print(curIdx) #debug
            # Using that index, pop the entry (removes from list)
            curBoxL, curBoxR = domainPatches.pop(curIdx)
            # In place of that entry, put a box truncated to fit inside the domain
            # as well as a periodic adjusted box for every dimension that exceeds
            #      truncated box
            tempBoxR=[]
            tempBoxL=[]
            tempBoxL.extend(curBoxL) # extend copies lists by value, not reference
            for i in range(len(tempBoxL)):
                if( curBoxR[i] > fdR[i] ):
                    tempBoxR.append(fdR[i])
                else:
                    tempBoxR.append(curBoxR[i])
            domainPatches.append( [tempBoxL,tempBoxR] )
            #      adjusted box(es)
            for i in range(len(curBoxL)):
                del tempBoxL, tempBoxR
                tempBoxL=[]
                tempBoxR=[]
                tempBoxL.extend(curBoxL)
                tempBoxR.extend(curBoxR)
                if(tempBoxR[i] > fdR[i]): # the i-th dimension needs adjustment only
                    tempBoxR[i] = tempBoxR[i] - (fdR[i]-fdL[i])
                    tempBoxL[i] = fdL[i] + leftAdjustment[i]
                    domainPatches.append( [tempBoxL,tempBoxR] )
        # clean up the patches in the case of duplicates or zero-volume boxes
        for box in domainPatches:
            for i in range(len(box[0])):
                if(box[1][i] - box[0][i] < mincellsize): # zero volume on grid (using mincell instead to avoid roundoff?)
                    print("Removing zero-volume patch")
                    domainPatches.pop(domainPatches.index(box))
                    break # this box is gone, next box please
        idxToDelete = []
        #print(domainPatches)
        for i in range(len(domainPatches)):
            for j in range(i+1,len(domainPatches)):
                if(idxToDelete.count(j)>0): continue
                box1 = domainPatches[i]
                box2 = domainPatches[j]
                if(box1==box2):
                    idxToDelete.append(j)
        #print(idxToDelete)
        idxToDelete.sort()
        idxToDelete.reverse()
        for i in idxToDelete:
            print("Removing duplicate patch")
            domainPatches.pop(i) # the reverse makes the changing size of the list not affect things
        # end of while loop
    print("End of breaking up!")
    print(domainPatches)


# -=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=-
#                           PROGRAM STARTS HERE
# -=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=-

# Greetings
Messages.Welcome()

# Gets box lengths ready (may overwrite the right-side ones)
boxLeft  = inputs.box_L #[inputs.box_xL,inputs.box_yL,inputs.box_zL]
boxRight = inputs.box_R #[inputs.box_xR,inputs.box_yR,inputs.box_zR]
boxLeft =  [ Convert(x,inputs.box_units,'cm','cm') for x in boxLeft] # convert to CGS
boxRight = [ Convert(x,inputs.box_units,'cm','cm') for x in boxRight] # convert to CGS

# Checks on input values before we read in file
if not any([inputs.is_periodic==y for y in [0,1]]):
    assert False, "is_periodic needs to be 0 or 1. C'mon!"
if not any([inputs.setup_type==y for y in [0,1]]):
    assert False, "setup_type needs to be 0 or 1. C'mon!"
assert inputs.aw_XYZ>0, "Atomic weight needs to be positive!"
assert (inputs.x_XYZ>0 and inputs.x_XYZ<1.0), "Mass fraction outside physical limits!"

# Loads file into YT
ds = yt.load(inputs.O2gas_fname)
try:
    Messages.Print("Loaded file " + str(ds))
except NameError:
    assert False, "Unable to properly load file into YT!"

# Domain left and right edges
fullDomain_L = ds.domain_left_edge
fullDomain_R = ds.domain_right_edge

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
            assert False, "Box physical dimensions are weird"
        inputs.Ncells[i] = int(np.ceil((boxRight[i] - boxLeft[i])/dxArray[i]))

# Checks & Balances
assert(inputs.max_level <= ds.max_level)
assert(all(np.array(inputs.Ncells) >= 2))
if(inputs.is_periodic):
    if(any( np.array(boxRight)-np.array(boxLeft) >  np.array(fullDomain_R.d)-np.array(fullDomain_L.d))):
        assert False, ("Carve out box is bigger than domain. Weird!")
else:
    if(any(boxLeft < fullDomain_L.d)):
        print((boxLeft-fullDomain_L.d)/fullDomain_L.d )
        if(any( (boxLeft-fullDomain_L.d)/fullDomain_L.d > 0.001 ) ):
            assert False, ("Left side of desired carve out outside actual domain. " + str(boxLeft) + " , " + str(fullDomain_L.d))
    if(any(boxRight > fullDomain_R.d)):
        print((boxRight-fullDomain_R.d)/fullDomain_R.d )
        if(any( (boxRight-fullDomain_R.d)/fullDomain_R.d > 0.001 ) ):
            assert False, ("Right side of desired carve out outside actual domain. " + str(boxRight) + " , " + str(fullDomain_R.d))


Messages.Print("Carving between Left = " + str(boxLeft))
Messages.Print("            to Right = " + str(boxRight))
Messages.Print("           w/ Ncells = " + str(inputs.Ncells))

# If periodic allowed, breaks the domain region into a set of patches we'll use to check overlap
domainPatches = [ [boxLeft,boxRight] ]
if(inputs.is_periodic):
    # boxLeft should always be within domain
    # boxRight could be beyond the physical domain
    # -- periodicity allows us to loop back around
    BreakIntoPatches(domainPatches,[fullDomain_L.d,fullDomain_R.d], dxArray.max() ) # domainPatches changes


# Call the writer class constructor (super fast)
writer = CarvingWriter(ds, boxLeft, boxRight, domainPatches, fullDomain_L, fullDomain_R, inputs.Ncells, inputs.max_level, inputs.is_periodic)


# Write the amr grid file (fast)
Messages.Print("Writing amr grid file")
writer.write_amr_grid()

# Write the number density file for species or dust (slow)
Messages.Print("Writing number density file")
writer.write_line_file(("gas", "number_density_XYZ"), inputs.out_nfname)

# Write the temperature file for species or dust (slow)
Messages.Print("Writing temperature file")
writer.write_line_file(("gas", "gas_temperature"), inputs.out_tfname)

# Write the gas velocity file (slow x 3)
Messages.Print("Writing velocity file")
velocity_fields = [inputs.velocityX_field,inputs.velocityY_field,inputs.velocityZ_field]
writer.write_line_file(velocity_fields, inputs.out_vfname)

# We're done!
Messages.Goodbye()
