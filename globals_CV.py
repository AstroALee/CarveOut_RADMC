"""
   globals_CarveOut.py

   Purpose:
        Constants and simple functions for use in the RadMC carving routines.
        Should not need to edit this file ever.

   Author:
        Aaron T. Lee, atl8@stmarys-ca.edu
        Spring 2020

   Written/Tested with Python 3.7.3, YT 3.4.1
"""

# CGS constants
class Csts:
    G  = 6.67408e-8 # gravity
    kB = 1.380658e-16 # Boltzmann constant
    mp = 1.6726231e-24 # proton mass
    c  = 2.9979245e10 # speed of light
    h  = 6.62607e-27 # Planck
    sbC = 5.6704e-5 # Stefan-Boltzmann constant

# Unit conversions
class UnitConv:
    cgsunits = ['cm','g','sec']
    lenunits = ['cm','au','ly','pc']
    massunits = ['sol','g']
    timeunits = ['sec','yr','Myr']
    au2cm = 1.496e13 # AU to cm
    pc2cm = 3.08567758128e18 # parsec to cm
    ly2cm = 9.463e18 # light-year to cm
    sol2g = 1.989e33 # solar mass to g
    sol2ergs = 3.83e33 # solar luminosity to erg/s
    sol2cm = 6.9550e10 # solar radius to cm
    yr2sec = 3.154e7   # year in seconds
    Myr2sec = 3.154e13 # Million years in seconds

# Unit conversion function
# Inputs: x = numerical value,
#         unit_in = unit of input x, unit_out = desired out unit
#         cgsunit = fundamental cgs unit (e.g., 'g' if mass, 'cm' if length)
def Convert(x,unit_in,unit_out,cgsunit):
    assert(cgsunit in UnitConv.cgsunits)

    if(cgsunit=='cm'):
        possible_units = UnitConv.lenunits
    elif(cgsunit=='g'):
        possible_units = UnitConv.massunits
    else:
        possible_units = UnitConv.timeunits

    units = [ unit_in.lower(), unit_out.lower()]
    if False in [x in possible_units for x in units]:
        print("Unit conversion failure for units: " + str(unit_in) + " and " + str(unit_out))

    # Convert first to CGS
    val1 = 1.0
    if units[0] not in UnitConv.cgsunits :
        val1 =eval( 'UnitConv.' + units[0] + '2' + cgsunit )

    # Convert to what you want
    val2 = 1.0
    if units[1] not in UnitConv.cgsunits :
        val2 = 1.0/eval( 'UnitConv.' + units[1] + '2' + cgsunit )

    return(x*val1*val2)


# Chatty Cathy
class Messages:
    # print function that only prints if your verbosity is high enough
    def Print(verb,v,s):
        if(verb >= v):
            print(s)
    # welcome
    def Welcome():
        print("\nWelcome to the carve out routine.")
        print("Depending on the file size,this could take a while.")
        print("Go refill your coffee.\n")
    # goodbye
    def Goodbye():
        print("\n\nAll done. Now get to work!")
