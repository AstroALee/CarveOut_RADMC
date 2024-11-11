"""
   yfields_CV.py

   Purpose:
        YT fields needed for the carving routines.

   Author:
        Aaron T. Lee, atl8@stmarys-ca.edu
        Spring 2020

   Written/Tested with Python 3.7.3, YT 3.4.1
"""

#import yt
import yt
from inputs_CV import inputs
from globals_CV import *

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
