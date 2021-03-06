"""
   main_CarveOut.py

   Purpose:
        Driver for RadMC carve routines. Call this file when you want
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

# Definition of the dust field. Uses dust to gas ration from inputs file
def _DustMassDensity(field, data):
    return inputs.dustgas_ratio*data[inputs.density_field]
yt.add_field(("gas", "dust_density"), function=_DustMassDensity, units="g*cm**-3")

# Definition of the gas temperature. Uses info from inputs
def _GasTemperature(field, data):
    ethbyden = data[inputs.ienergy_field]/data[inputs.density_field]
    gmmfac = (inputs.gas_gamma - 1.0)*yt.YTQuantity(inputs.gas_mmw*Csts.mp, 'g')/(yt.units.kb)
    Temp = gmmfac*ethbyden
    Temp[Temp > inputs.max_temp] = yt.YTQuantity(inputs.max_temp, 'K') # overwrites spurious high T
    return (Temp)
yt.add_field(("gas", "gas_temperature"), function=_GasTemperature, units="K")

# Definition of the dust temperature. Assumes same as gas for now
def _DustTemperature(field, data):
    return data["gas_temperature"]
yt.add_field(("gas", "dust_temperature"), function=_DustTemperature, units="K")

# Definition of the XYZ species field. Uses info from inputs
def _NumberDensityXYZ(field, data):
    aw = yt.YTQuantity(inputs.gas_mmw*Csts.mp, 'g')
    H2numDen = (inputs.x_H2/aw)*data[inputs.density_field]
    XYZnumDen= (inputs.x_XYZ/aw)*data[inputs.density_field]
    #print(H2numDen.min())
    #print(H2numDen.max())
    if(inputs.allow_freezeout):
        XYZnumDen[H2numDen < inputs.freeze_minN] = yt.YTQuantity(0.0, "cm**-3")
        XYZnumDen[H2numDen > inputs.freeze_maxN] = yt.YTQuantity(0.0, "cm**-3")
    return XYZnumDen
yt.add_field(("gas", "number_density_XYZ"), function=_NumberDensityXYZ, units="cm**-3")

# Definition of the H_2 species field. Uses info from inputs
def _NumberDensityH2(field, data):
    aw = yt.YTQuantity(inputs.gas_mmw*Csts.mp, 'g')
    return (inputs.x_H2/aw)*data[inputs.density_field]
yt.add_field(("gas", "number_density_H2"), function=_NumberDensityH2, units="cm**-3")

# Definition of the microturbulence at each point. Uses info from inputs
def _MicroTurb(field, data):
    turb = data[inputs.velocityX_field]
    turb[turb>=0] = yt.YTQuantity(inputs.microturbulence_speed, "cm/s")
    turb[turb<=0] = yt.YTQuantity(inputs.microturbulence_speed, "cm/s")
    #print(turb.min())
    #print(turb.max())
    return turb
yt.add_field(("gas", "microturbulence_speed"), function=_MicroTurb, units="cm/s")


# Routine that, when periodic boxes are allowed, breaks a domain that extends
# beyond the edges of the physical domain into a discrete set of patches equivalent
# to if periodic boundary conditions were employed.
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
        Messages.Print(2,"Beginning of while loop, current domainPatches:")
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
            Messages.Print(2,"All patches inside domain")
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
                    Messages.Print(1,"Removing zero-volume patch")
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
            Messages.Print(1,"Removing duplicate patch")
            domainPatches.pop(i) # the reverse makes the changing size of the list not affect things
        # end of while loop
    Messages.Print(0,"Domain broken into patches [ [LE,RE]_i ]:")
    Messages.Print(0,domainPatches)


# -=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=-
#                           PROGRAM STARTS HERE
# -=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=--=-=-

# Greetings
Messages.Welcome()

# Gets box lengths ready (may overwrite the right-side ones)
boxLeft  = np.array(inputs.box_L) #[inputs.box_xL,inputs.box_yL,inputs.box_zL]
boxRight = np.array( len(boxLeft)*[0] )
boxLeft =  [ Convert(x,inputs.box_units,'cm','cm') for x in boxLeft] # convert to CGS
boxRight = [ Convert(x,inputs.box_units,'cm','cm') for x in boxRight] # convert to CGS

# Checks on input values before we read in file
if not any([inputs.is_periodic==y for y in [0,1]]):
    assert False, "is_periodic needs to be 0 or 1. C'mon!"
if not any([inputs.is_unigrid==y for y in [0,1]]):
    assert False, "is_unigrid needs to be 0 or 1. C'mon!"
if not any([inputs.allow_freezeout==y for y in [0,1]]):
    assert False, "allow_freezeout needs to be 0 or 1. C'mon!"
if not any([inputs.has_microturb==y for y in [0,1]]):
    assert False, "has_microturb needs to be 0 or 1. C'mon!"
assert (inputs.x_XYZ>0 and inputs.x_XYZ<1.0), "Number fraction outside physical limits!"
assert (inputs.freeze_minN >= 0), "Really....?"
assert (inputs.freeze_maxN >= inputs.freeze_minN), "Really....?"


# Loads file into YT
ds = yt.load(inputs.O2gas_fname)
try:
    print("Loaded file " + str(ds)) # using classic print() for try/except to work
except NameError:
    assert False, "Unable to properly load file into YT!"

# Domain left and right edges
fullDomain_L = ds.domain_left_edge
fullDomain_R = ds.domain_right_edge

# Cell size on coarsest level (will assume factors of 2 refinement)
dxArray = ds.index.dds_list[0]

# Using cell size, sets the right edge of the box
boxRight = boxLeft + np.array(inputs.Ncells)*dxArray

# Checks & Balances
assert(inputs.verbosity>=0)
assert(inputs.max_level <= ds.max_level)
assert(all(np.array(inputs.Ncells) >= 2))
if(inputs.is_periodic):
    if(any( np.array(boxRight)-np.array(boxLeft) >  np.array(fullDomain_R.d)-np.array(fullDomain_L.d))):
        assert False, ("Without periodic turned on: Carve out box is bigger than domain. Weird!")
else:
    if(any(boxLeft < fullDomain_L.d)):
        print((boxLeft-fullDomain_L.d)/fullDomain_L.d )
        if(any( (boxLeft-fullDomain_L.d)/fullDomain_L.d > 0.001 ) ):
            assert False, ("Left side of desired carve out outside actual domain. " + str(boxLeft) + " , " + str(fullDomain_L.d))
    if(any(boxRight > fullDomain_R.d)):
        print((boxRight-fullDomain_R.d)/fullDomain_R.d )
        if(any( (boxRight-fullDomain_R.d)/fullDomain_R.d > 0.001 ) ):
            assert False, ("Right side of desired carve out outside actual domain. " + str(boxRight) + " , " + str(fullDomain_R.d))


Messages.Print(0,"Carving between Left = " + str(boxLeft))
Messages.Print(0,"            to Right = " + str(boxRight))
Messages.Print(0,"           w/ Ncells = " + str(inputs.Ncells))

# If periodic allowed, breaks the domain region into a set of patches we'll use to check overlap
domainPatches = [ [boxLeft,boxRight] ]
if(inputs.is_periodic):
    # boxLeft should always be within domain
    # boxRight could be beyond the physical domain
    # -- periodicity allows us to loop back around
    BreakIntoPatches(domainPatches,[fullDomain_L.d,fullDomain_R.d], dxArray.max() ) # domainPatches changes


# Call the writer class constructor (super fast)
writer = CarvingWriter(ds, boxLeft, boxRight, domainPatches, fullDomain_L, fullDomain_R, inputs.Ncells, inputs.max_level, inputs.is_periodic, inputs.is_unigrid)


# Write the amr grid file (fast)
Messages.Print(0,"1/7: Writing amr grid file (fast!)")
writer.write_amr_grid()

# Write the number density file for species (slow)
Messages.Print(0,"2/7: Writing number density file (slow.)")
writer.write_line_file(("gas", "number_density_XYZ"), inputs.out_nfname)

# Write the number density file for species (slow)
Messages.Print(0,"3/7: Writing number density file for H_2 (slow.)")
writer.write_line_file(("gas", "number_density_H2"), inputs.out_h2fname)

# Write the dust density file for dust (slow)
Messages.Print(0,"4/7: Writing dust density file (slow.)")
writer.write_dust_file(("gas", "dust_density"), inputs.out_ddfname)

if(inputs.has_microturb):
    Messages.Print(0,"4.5/7: Writing microturbulence file (slow.)")
    writer.write_line_file(("gas", "microturbulence_speed"), inputs.out_mtfname)

# Write the temperature file for species or dust (slow)
Messages.Print(0,"5/7: Writing temperature file (slower..)")
writer.write_line_file(("gas", "gas_temperature"), inputs.out_tfname)

# Assuming dust temperature is same as gas for now...
Messages.Print(0,"6/7: Writing dust temperature file (slower..)")
writer.write_dust_file(("gas", "dust_temperature"), inputs.out_dtfname)

# Write the gas velocity file (slow x 3)
Messages.Print(0,"7/7: Writing velocity file (slowest...)")
velocity_fields = [inputs.velocityX_field,inputs.velocityY_field,inputs.velocityZ_field]
writer.write_line_file(velocity_fields, inputs.out_vfname)

# We're done!
Messages.Goodbye()
