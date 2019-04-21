# /***********************************************************************
# Copyright (c) Fisheries and Oceans Canada, 2010
#
# This program may be freely redistributed for non-commercial purposes
# under the condition that the copyright notices are not removed. You may
# distribute modified versions of this program under the conditions that
# both source and object code are made available without charge and that
# clear notice is given of the modifications. Modified programs and source
# code are not to be represented as being endorsed by, or made in
# cooperation with, Fisheries and Oceans Canada.
#
# This program is distributed by the copyright holder and contributors in
# the hope that it will be useful. It is provided "AS IS" without
# warranties of any kind, whether express or implied. This includes, but
# is not limited to, merchantability or fitness for a particular purpose.
#
# In no event shall the copyright holder be liable for damages, including
# any general, special, incidental or consequential damages arising out of
# the use or inability to use the program (including, but not limited to,
# loss of use, data or profits).
#
# This program is not to be used for navigation purposes.
#
# This program is part of the WebTide software package.
# /**********************************************************************/
#
# /* PROGRAM: tidecor
#    VERSION: 2.42
#    DATE :   January 14, 2004
#    Author : Jason Chaffey
#    Modifications: Shawn Oakey, Jason Chaffey
# */
#
# /* Version 1.0 */
# /*******************************************************************************
#  * This program for tide correction of hydrographic data is derived from
#  * interp2d.c - Modifications and additions made by Patrick Roussel May 21,1999
#  ******************************************************************************/
#
# /* Version 2.0 */
# /* The FORTRAN components of Ver1.0 were translated to c for portability.
#    Original FORTRAN written by M. Foreman at IOS.
#    RAYBOUND bug fixed.
#    Jason Chaffey              December 15, 2000
# */
#
# /* Version 2.1 */
# /* Added configuration file (tidecor2.1.cfg) so that the following are not
#    hard-wired into the program:
#        - mesh filenames (.nod and .ele);
#        - the number of constituents to be used; and
#        - for each constituent:
#              - name of constituent; and
#              - filename for data of that constituent.
#    Jason Chaffey              January 15, 2001
# */
#
# /* Version 2.2 */
# /* Added support for complex input data (.v2c files)
#    Jason Chaffey              March 23, 2001
# */
#
# /* Version 2.21 */
# /* Made various changes for use in WebTide
#  *
#  * Web(Tide/Drogue)
#  *
#  * Ocean Science Division
#  * Bedford Institute of Oceanography
#  * Fisheries and Oceans Canada
#  *
#  * Shawn Oakey
# */
#
# /* Version 2.3 */
# /*
#  * Tested against Pawlawicz's T_Tide and
#  * Foreman's tide4_r2
#  *
# */
#
# /* Version 2.4 */
# /*
#  * Debugging (float - double and int - long int conflicts) done
#  * Switched to Manhattan distance vs "real" distance in closestnode
#  *          - speed increase consideration
#  * Thanks to Herman Varma for the above
#  * Jason Chaffey              August 7, 2003
#  *
# */
# /* Version 2.41 */
# /*
#  * Returned to using "real" distance calculations versus Manhattan distance.
#  * Problems were discovered with Manhattan distance with the NE Pacific mesh.
#  * Jason Chaffey              January 6, 2004
#  *
# */
# /* Version 2.42 */
# /*
#  * Bug in Manhattan distance calculation fixed. Returned to using it
#  * for optimization.
#  * Added fix for elements that cross International Dateline or
#  * Greenwich Meridian in raybound.
#  * Jason Chaffey              January 14, 2004
#  *
# */
# /* Version 2.5 */
# /*
#  * Fixed some problems with raybound for elements that cross International
#  * Dateline or Greenwich Meridian in raybound. Problems became evident in the
#  * global data set.
#  * Jason Chaffey              June 16, 2008
#  *
# */
# /* Version 2.5.1 */
# /*
#  * Fixed numbering of nodes in RAYBOUND after old bug snuck back in.
#  * Jason Chaffey              June 26, 2010
#  *
# */
# /* Version AGRG */
# /*******************************************************************************
#  * This script was translated into Python from the original webtidecor2.5.1.c
#  * Modifications and additions made by Andrea Zagar - April, 2019
#  ******************************************************************************/


#!/usr/bin/env python

import math
import sys
import datetime as dt
import numpy as np

#------------------------------------------------------------------------------------

# The following are structures used to hold the info of the constituents
# First is the structure to hold the satellites for each main constituent
class SatType(object):
    __slots__ = ('deld','phase','ratio','corr')

    def __init__(self, deld, phase, ratio, corr):
        self.deld = deld
        self.phase = phase
        self.ratio = ratio
        self.corr = corr

# This is for the main constituents
class MainType(object):
    __slots__ = ('name','dood','phase','nsats','sats','f','freq','vu')

    def __init__(self, name, dood, phase, nsats):
        self.name = name
        self.dood = dood
        self.phase = phase
        self.nsats = nsats
        self.sats = []

# This is for the main constituent factors for the shallow water constituents
class ShallowConstType(object):
    __slots__ = ('name','factor')

    def __init__(self, name, factor):
        self.name = name
        self.factor = factor

# This is for shallow water constituents
class ShallowType(object):
    __slots__ = ('name','numcon','cons', 'f', 'freq', 'vu')

    def __init__(self, name, numcon):
        self.name = name
        self.numcon = numcon
        self.cons = []

#------------------------------------------------------------------------------------

# Global constants:
R = 6.3675e-8           # Radius of the earth
nconst = 5              # Number of constituents

# Define the value for 2 * pi
twopi = 2.0 * math.acos( -1.0 )

# This defines the (Julian) date of the switch to the Gregorian calendar
# (Used in the julday subroutine)
IGREG = (15 + 31 * (10 + 12 * 1582 ))

# These define the indices for the 5 constituents used in this program
# (Perhaps someday users will be allowed to choose which and how many to use)
m2 = 0
s2 = 1
n2 = 2
k1 = 3
o1 = 4

# Table of number of days for each month
daytable = [
    [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
    [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]]

# Global variables:
# global nodesArr, elemsArr
nodesList, elemsList = [],[]
numNodes, numElems = 0,0
numSatAttr = 5
#nodeMinTidesList = []

# Dictionaries for reading iosfile (stores ALL possible constituents in those dictionaries)
mainConstDict, shallConstDict = {}, {}

# Dictionary of main constituents used in calculations (around 5)
constDict, constFilePathsDict = {}, {}

#------------------------------------------------------------------------------------
#
# FUNCTIONS
#
#------------------------------------------------------------------------------------

# Calculates the following ephermides of the sun and moon:
# h  = mean longitude of the sun;
# pp = mean longitude of the solar perigee;
# s  = mean longitude of the moon;
# p  = mean longitude of the lunar perigee; and
# np = negative of the longitude of the mean ascending node.
# Also calculates their rate of change.
# Units are cycles ( cycles / 365 days for rates of change ).
# These formulae were taken from pp 98 and 107 of the Explanatory
# Supplement to the Astronomical Ephermeris and the American
# Ephermis and Nautical Almanac (1961)
def astro_angles(d1):

    # Local variables to calculate
    d12 = d1 * d1
    d2 = d1 * 1.0E-04
    d22 = d2 * d2
    d23 = pow( d2, 3.0 )
    f = 360.0
    f2 = f / 365.0

    # Calculate global variables
    h = ( 2.79696678E+02 + d1 * 9.856473354E-01 + d22 * 2.267E-05 ) / f
    h = math.modf(h)[0]

    pp = ( 2.81220833E+02 + d1 * 4.70684E-05 + d22 * 3.39E-05 + d23 * 7.0E-08 ) / f
    pp = math.modf(pp)[0]

    s = ( 2.70434164E+02 + d1 * 1.31763965268E+01 - d22 * 8.5E-05 + d23 * 3.9E-08 ) / f
    s = math.modf(s)[0]

    p = ( 3.34329556E+02 + d1 * 1.114040803E-01 - d22 * 7.739E-04 - d23 * 2.6E-07 ) / f
    p = math.modf(p)[0]

    np = ( -2.59183275E+02 + d1 * 5.29539222E-02 - d22 * 1.557E-04 - d23 * 5.0E-08 ) / f
    np = math.modf(np)[0]

    dh = ( 9.856473354E-01 + d1 * 2.267E-05 * 2.0E-08 ) / f2
    dpp = ( 4.70684E-05 + d1 * 3.39E-05 * 2.0E-08 + d12 * 7.0E-08 * 3.0E-12 ) / f2
    ds = ( 1.31763965268E+01 - d1 * 8.5E-05 * 2.0E-08 + d12 * 3.9E-08 *3.0E-12 ) / f2
    dp = ( 1.114040803E-01 - d1 * 7.739E-04 * 2.0E-08 - d12 * 2.6E-07 * 3.0E-12 ) / f2
    dnp = ( 5.29539222E-02 - d1 * 1.557E-04 * 2.0E-08 - d12 * 5.0E-08 * 3.0E-12 ) / f2

    return h, pp, s, p, np, dh, dpp, ds, dp, dnp

#------------------------------------------------------------------------------------
# Gets satellite info from input line, puts the info into a new satellite structure and returns the structure
def add_sat(satin):

    return 'need to finish this'

#------------------------------------------------------------------------------------
# Finds the closest node to a point (ptx, pty) and the element containing
# that point, if one exists.
# Also gets the basis functions for interpolations to that point.
# Returns the element number containing the point (ptx, pty).
# Returns -999 if no element found containing the point (ptx, pty ).
def basis2d(ptx, pty):

    # Find node closest to point
    closest = closestNode( numNodes, ptx, pty )
    print("closest node: ", closest)

    # Try to find an element that contains the point and the closest node
    ele = -999
    xCoords, yCoords, basis = [],[],[]

    # Read in node vals for each element (subtract 1 because index of nodesDf starts at 0 and not 1)
    for index in range(len(elemsList)):
        row = elemsList[index]
        n11, n22, n33 = row[0]-1, row[1]-1, row[2]-1
        xCoords, yCoords = [],[]
        # If one or more of the nodes is in this element
        if (( closest == (n11) ) or ( closest == (n22) ) or ( closest == (n33) )):
            #print(n11, n22, n33)
            # Store lat and long of each node
            xCoords.extend((nodesList[n11][0], nodesList[n22][0], nodesList[n33][0]))
            yCoords.extend((nodesList[n11][1], nodesList[n22][1], nodesList[n33][1]))
            # See if the point is within this element
            flag, xCoords, yCoords = raybound(xCoords, yCoords, ptx, pty)
            #print("flag: ", flag)
            if flag == 1: # The point is within the element
                print("n11: ",n11,", n22: ", n22,", n33: ", n33)
                ele = index
                basis = phi2d(xCoords, yCoords, ptx, pty)
                #print("basis: ", basis)
                break
    # If the closest node's elements don't work, search through all elements
    if (ele < 0):
        print("searching through all elements")
        for index in range(len(elemsList)):
            row = elemsList[index]
            n11, n22, n33 = row[0]-1, row[1]-1, row[2]-1
            xCoords, yCoords = [],[]
            # Store lat and long of each node
            xCoords.extend((nodesList[n11][0], nodesList[n22][0], nodesList[n33][0]))
            yCoords.extend((nodesList[n11][1], nodesList[n22][1], nodesList[n33][1]))
            # See if the point is within this element
            flag, xCoords, yCoords = raybound(xCoords, yCoords, ptx, pty)
            if flag == 1: # The point is within the element
                print("n11: ",n11,", n22: ", n22,", n33: ", n33)
                ele = index
                basis = phi2d(xCoords, yCoords, ptx, pty)
                return ele

    return ele, basis

#------------------------------------------------------------------------------------
# Find the node that is closest to the point (ptx, pty).
def closestNode(numnodes, ptx, pty):

    # Keep track of closest index
    closeIdx = 0
    closeDist = 0

    # Go through node dataframe
    for index in range(len(nodesList)):
        row = nodesList[index]
        if index == 0:
            closeDist = math.fabs(pty - row[1]) + math.fabs(ptx - row[0])
        else:
            # Use Manhattan distance for calculation
            currentDist = math.fabs(pty - row[1]) + math.fabs(ptx - row[0])
            if currentDist < closeDist:
                closeDist = currentDist
                closeIdx = index

    return closeIdx

#------------------------------------------------------------------------------------
# Get the day and month from the day # of the year
def get_date(year, dayofyear):

    leap = (( year % 4 == 0 ) and ( year % 100 != 0 ) or ( year % 400 == 0 ))
    i = 1
    while (dayofyear > daytable[leap][i]):
        dayofyear -= daytable[leap][i]
        i += 1

    month = i
    day = dayofyear

    return month, day

#------------------------------------------------------------------------------------
# Calculate the Julian day number. Accounts for the change to the Gregorian calandar.
def julday(id, im, iy):

    jy = iy
    if( jy == 0 ):
        print("JULDAY: There is no year 0!\n" )
        sys.exit(2)

    if( jy < 0 ):
        jy += 1
    if( im > 2 ):
        jm = im + 1
    else:
        jy -= 1
        jm = im + 13

    jul = int(math.floor(365.25 * jy)) + math.floor(30.6001 * jm) + id + 1720995

    if(( id + 31 * ( im + 12 * iy )) >= IGREG ):
        ja = int(0.01 * jy)
        jul += 2 - ja + int(0.25 * ja)


    return jul

#------------------------------------------------------------------------------------
# Read in the constituent data from the IOS_tidetbl file
def openvuf(iosfilePath):

    # Global vars
    global mainConstDict, shallConstDict
    mainConstDict, shallConstDict = {}, {}

    # Open the ios table file
    iosfile = open(iosfilePath, "r+")

    # Get main constituents
    while (True):

        # Get line from file
        line = iosfile.readline().rstrip("/n").split()
        # Check for blank line to signal end of main constituents
        lineLen = len(line)
        if lineLen <= 1:
            break
        # Check to make sure the phase isn't accidentally joined with a dood number in the file
        try:
            int(line[6])
        except ValueError:
            checkPhase = line[6].split("-")
            #print("phase error: ", checkPhase)
            # Insert negative sign in front of both numbers
            if len(checkPhase) > 2:
                line[6] = "-" + checkPhase[2]
                line.insert(6, ("-" + checkPhase[1]))
            # Insert negative sign in front of second number only
            else:
                line[6] = "-" + checkPhase[1]
                line.insert(6, checkPhase[0])

        # Make new constituent object
        doods = list(map(int, line[1:7]))
        newMain = MainType(line[0], doods, float(line[7]), int(line[8]))

        # Get number of satellites, if above 0 then loop through satellites, assign to object
        if newMain.nsats > 0:
            satCount = newMain.nsats
            #print("initial satellites: ", satCount)
            while (satCount > 0):
                # Get next line and line length
                line = iosfile.readline().rstrip("/n").split()
                #print(line)
                lineLen = len(line)
                # Loop through sats and make new objects and add to sats
                numLineSats = int((lineLen-1) / numSatAttr)
                #print("current satellites: ", numLineSats)
                for x in range(numLineSats):
                    # Get deld array
                    deldIdx = 1 + x*numSatAttr
                    deld = list(map(int, line[deldIdx:deldIdx+3])) # 1-4, 6-9
                    # Phase value
                    phase = line[(4 + x*numSatAttr)]
                    # Ratio value
                    ratio = line[(5 + x*numSatAttr)]
                    # Get and assign corr value
                    corr = ratio[-2:]
                    if corr == "R1":
                        corrValue = 1
                        ratio = ratio[:-2]
                    elif corr == "R2":
                        corrValue = 2
                        ratio = ratio[:-2]
                    else:
                        corrValue = 0
                        ratio = float(ratio)

                    # Make new sat object and add to sats
                    newSat = SatType(deld, float(phase), float(ratio), corrValue) # *** convert phase to float?
                    #print("inserted new satellite into ", newMain.name)
                    newMain.sats.append(newSat)



                satCount -= numLineSats

        # Save main constituent into dictionary
        mainConstDict[newMain.name] = newMain

    print("number of main constituents: ", len(mainConstDict))

    # Get shallow constituents
    while (True):
        # Get line from file
        line = iosfile.readline().rstrip("/n").split()
        # Check for blank line to signal end of main constituents
        lineLen = len(line)
        if lineLen <= 1:
            break
        # Make new shallow constituent object
        numLineFactors = int(line[1])
        newShallow = ShallowType(line[0], numLineFactors)
        line = line[2:]
        # Load factors into objects and insert into shallow constituent object
        for x in range(numLineFactors):
            idx = x*2
            newShallowFactor = ShallowConstType(line[idx+1], float(line[idx]))
            newShallow.cons.append(newShallowFactor)

        # Save shallow constituent into dictionary
        shallConstDict[newShallow.name] = newShallow


    print("num shallow constituents: ", len(shallConstDict))
    # Close file
    iosfile.close()

#------------------------------------------------------------------------------------
#Calculates the basis functions for interpolating to a point inside an element.
def phi2d(xloc, yloc, ptx, pty):

    # Declare variables
    phi = []

    #print("phi xpoly: ", xloc)

    # Calculate area
    area = 0.5 * ( xloc[0] * ( yloc[1] - yloc[2] ) +
                xloc[1] * ( yloc[2] - yloc[0] ) +
                xloc[2] * ( yloc[0] - yloc[1] ))

    # Calculate the basis function
    for i in range(3):
        if i == 0:
            j = 1
            k = 2
        elif i == 1:
            j = 2
            k = 0
        else:
            j = 0
            k = 1

        a = ( xloc[j] * yloc[k] - xloc[k] * yloc[j] ) / ( area * 2 )
        b = ( yloc[j] - yloc[k] ) / ( area * 2 )
        c = -1 * ( xloc[j] - xloc[k] ) / ( area * 2 )

        #print("phi basis: ")
        #print(a, b, c)

        phi.append(a + b * ptx + c * pty)

    return phi

#------------------------------------------------------------------------------------
# Subroutine to check whether or not a point is inside a polygon.
# The process is as follows:
# 	Use an arbitrary ray (here, y = constant and x >= xref), starting from
#   the point and going off to infinity.
# 	Count the number of polygon boundaries it crosses.
# 	If an odd number, the point is inside the polygon, otherwise it is outside.
def raybound(xpoly, ypoly, ptx, pty):

    bcross = 0  # Number of boundary crossings

    #  Check to see if the element side crosses the International Dateline
    # (changes sign at +180/-180 degrees) and if so, change the longitudes
    # so that they all have the same sign.

    if ( ptx > 0.0 ):

        if(( xpoly[0] < -170.0 ) and ((xpoly[1] > 170.0 ) or ( xpoly[2] > 170.0 ))):
            xpoly[0] += 360.0
        if(( xpoly[1] < -170.0 ) and ((xpoly[0] > 170.0 ) or ( xpoly[2] > 170.0 ))):
            xpoly[1] += 360.0
        if(( xpoly[2] < -170.0 ) and ((xpoly[1] > 170.0 ) or ( xpoly[0] > 170.0 ))):
            xpoly[2] += 360.0

    else:
        if(( xpoly[0] > 170.0 ) and ((xpoly[1] < -170.0 ) or ( xpoly[2] < -170.0 ))):
            xpoly[0] -= 360.0
        if(( xpoly[1] > 170.0 ) and ((xpoly[0] < -170.0 ) or ( xpoly[2] < -170.0 ))):
            xpoly[1] -= 360.0
        if(( xpoly[2] > 170.0 ) and ((xpoly[1] < -170.0 ) or ( xpoly[0] < -170.0 ))):
            xpoly[2] -= 360.0

    # As above, except for the Greenwich meridian, for longitude coordinates
    # that range from 0 to 360 degrees.
    if( ptx > 350.0 ):
        if(( xpoly[0] < 10.0 ) and ((xpoly[1] > 350.0 ) or ( xpoly[2] > 350.0 ))):
            xpoly[0] += 360.0
        if(( xpoly[1] < 10.0 ) and ((xpoly[0] > 350.0 ) or ( xpoly[2] > 350.0 ))):
            xpoly[1] += 360.0
        if(( xpoly[2] < 10.0 ) and ((xpoly[1] > 350.0 ) or ( xpoly[0] > 350.0 ))):
            xpoly[2] += 360.0

    else:
        if(( xpoly[0] > 350.0 ) and ((xpoly[1] < 10.0 ) or ( xpoly[2] < 10.0 ))):
            xpoly[0] -= 360.0
        if(( xpoly[1] > 350.0 ) and ((xpoly[0] < 10.0 ) or ( xpoly[2] < 10.0 ))):
            xpoly[1] -= 360.0
        if(( xpoly[2] > 350.0 ) and ((xpoly[1] < 10.0 ) or ( xpoly[0] < 10.0 ))):
            xpoly[2] -= 360.0

    # Loop over each line segment around the element (3 total)
    for i in range(3):
        # If both endpoints of the line segment are on the same (vertical)
        # side of the ray, do nothing. Otherwise, count the number of times the ray intersects the segment.
        j = 0 if ( i == 2 ) else i + 1

        if(not( (( ypoly[i] < pty ) and ( ypoly[j] < pty )) or
                (( ypoly[i] >= pty ) and ( ypoly[j] >= pty )) ) ):

            if(xpoly[i] != xpoly[j]):
                m = ( ypoly[j] - ypoly[i] ) / ( xpoly[j] - xpoly[i] )
                b = ypoly[i] - m * xpoly[i]
                x = ( pty - b ) / m
                if( x > ptx ):
                    bcross += 1
            else:
                if( xpoly[j] > ptx ):
                    bcross += 1

    # Return the evenness/oddness of the boundary crossings i.e. the remainder from division by two.
    return bcross % 2, xpoly, ypoly

#------------------------------------------------------------------------------------
# Read in the mesh from a .nod/.ele set of files.
def readNodes(nodeFile, elementFile):

    # These frames have to have indexes corresponding to element file indexes (starting at 1)
    global nodesArr, elemsArr, numNodes, numElems
    global nodesList, elemsList

    # Open the node file
    nodefile = open(nodeFile, "r+")
    # Loop through file
    while (True):
        # Get line from file
        line = nodefile.readline().rstrip("/n").split()
        # Check for blank line to signal end of file
        if len(line) <= 1:
            break
        # Append long-lat pair to node list
        nodesList.append([float(line[1]), float(line[2])]) #long, lat
    nodefile.close()

    # Convert to numpy array
    #nodesArr = np.asarray(nodesList)

    # Open the element file
    elementfile = open(elementFile, "r+")
    # Loop through file
    while (True):
        # Get line from file
        line = elementfile.readline().rstrip("/n").split()
        # Check for blank line to signal end of main constituents
        if len(line) <= 1:
            break
        # Append nodes to element list
        elemsList.append([int(line[1]), int(line[2]), int(line[3])])

    elementfile.close()

    # Convert to numpy array
    #elemsArr = np.asarray(elemsList)

    # Number of nodes and elements
    numNodes = len(nodesList)
    numElems = len(elemsList)

    return

#------------------------------------------------------------------------------------
# Calculate the amplitudes, phases, etc. for each of the constituents
def setvuf(kh, xlat):

    # Global and local variables
    global mainConstDict, shallConstDict
    dd = []
    adj = 0
    kd = julday(31, 12, 1899) # this is a date
    d1 = (float)(kh - kd) - 0.5
    ktmp = kh * 24

    # Compute things
    h, pp, s, p, np, dh, dpp, ds, dp, dnp = astro_angles(d1)
    #print(h, pp, s, p, np, dh, dpp, ds, dp, dnp)
    hh = (float)(ktmp) - ( math.floor( (float)( ktmp / 24.0 )) * 24.0 )
    tau = hh / 24.0 + h - s
    dtau = 365.00 + dh - ds
    slat = math.sin(( twopi / 2.0 ) * ( xlat / 180.0 ))

    # The main constituents
    # Loop through main constituents
    for const in mainConstDict.keys():
        # Reset variables
        dd = []
        # Get object
        mainObj = mainConstDict[const]
        # Get dood numbers
        for i in range(6):
            dd.append(mainObj.dood[i])
        # Freq and phase
        mainObj.freq = (dd[0] * dtau + dd[1] * ds + dd[2] * dh + dd[3] *dp + dd[4] * dnp + dd[5] * dpp) / ( 24.0 * 365.0)
        vdbl = dd[0] * tau + dd[1] * s + dd[2] * h + dd[3] * p + dd[4] * np + dd[5] * pp + mainObj.phase
        v = vdbl - (math.floor( math.floor( vdbl ) / 2.0 ) * 2.0 )
        sumc = 1.0
        sums = 0.0

        #print(mainObj.name)

        # Loop through satellites
        for sat in mainObj.sats:
            if sat.corr == 0:
                adj = sat.ratio
            elif sat.corr == 1:
                adj = sat.ratio * 0.36309 * ( 1.0 - 5.0 * slat * slat ) / slat
            elif sat.corr == 2:
                adj = sat.ratio * 2.59808 * slat
            else:
                print("Unknown corr value!!!")
                sys.exit(2)
            # Calculations
            uudbl = float(sat.deld[0]) * p + float(sat.deld[1]) * np + float(sat.deld[2]) * pp + float(sat.phase)
            uu = math.modf(uudbl)[0]
            sumc += ( adj * math.cos( uu * twopi ))
            sums += ( adj * math.sin( uu * twopi ))

        # Insert values
        mainObj.f = math.sqrt(( sumc * sumc ) + ( sums * sums ))
        mainObj.vu = v + math.atan2( sums, sumc ) / twopi
        mainConstDict[const] = mainObj
        #print("f: ",mainObj.f, "      vu: ", mainObj.vu)

    # The shallow water constituents
    for const in shallConstDict.keys():
        # Temporary object, reassign at end of loops
        shallObj = shallConstDict[const]
        # Assign values
        shallObj.f = 1.0
        shallObj.vu = 0.0
        shallObj.freq = 0.0
        # Go through each factor
        for factor in shallObj.cons:
            # Go through all main constituents
            for mainCon in mainConstDict.keys():
                # Make temp object
                mainObj = mainConstDict[mainCon]
                if mainCon == factor.name:
                    shallObj.f *= math.pow(mainObj.f, math.fabs(factor.factor))
                    shallObj.vu += (factor.factor * mainObj.vu)
                    shallObj.freq += (factor.factor * mainObj.freq)
                    break
        # Reassign with new values
        #print(const, ":   F = ", shallObj.f) #, " VU = ", shallObj.vu, " FREQ = ", shallObj.freq)
        shallConstDict[const] = shallObj

    return

#------------------------------------------------------------------------------------
# double TideP( long int day, long int month, long int year, long int hour,
#        long int minute, double seconds, double latitude, double *ampl,
#        double *phase, char **innames, int numnames, mainptr *cons,
#        shallptr *shall, int nmain, int nshall )
# Calculates and returns the tidal correction
def TideP(day, month, year, hour, minute, seconds, latitude, nodeIdx, isMultiResult):

    # First, get the julday
    kd = julday(day, month, year)

    # Set which variables to use
    phaseRes, ampRes, mintideRes = 'phase','amp','min_tide'
    if isMultiResult:
        phaseRes = 'phase2'
        ampRes = 'amp2'
        mintideRes = 'min_tide2'

    # Calculate the amplitudes, phases, etc. for each of the constituents
    setvuf(kd, latitude)

    # More calculations
    dthr = (( float(hour) * 3600.0 ) + ( float(minute) * 60.0 ) + float(seconds) ) / 3600.0
    res = 0.0

    # For each of the desired constituents..(See top of program)
    for constName in constDict.keys():
        # Find the constituent from those loaded from IOS_tidetbl
        # Check main constituents first
        if constName in mainConstDict:
            # It's a main constituent
            constObj = mainConstDict[constName]
        else: # Check shallow water constituents next
            if constName in shallConstDict:
                # It's a shallow water constituent
                constObj = shallConstDict[constName]
            else:
                print("Couldn't find constituent! Exiting...")
                sys.exit(2)

        # Calculate phase and amp
        phase = constDict[constName][phaseRes][nodeIdx]
        amp = constDict[constName][ampRes][nodeIdx]
        # Calculations
        revgmt = constObj.freq * dthr + constObj.vu - phase / 360.0 #phase[i] / 360.0
        #print("node #: ", nodeIdx, " .. ", constName, ": PHASE = ", phase)
        radgmt = twopi * math.modf(revgmt)[0]
        res += constObj.f * amp * math.cos(radgmt)

    return res

#------------------------------------------------------------------------------------
# int vuf( char *inname, mainptr *cons, shallptr *shall, int nmain, int nshall )
# Finds constituent info corresponding to inname and returns the index to
# the node containing the info. Shallow water constituent indices are
# returned as their number greater than the max # of main constituents.
# e.g. if we want the 2nd shallow water constituent and there are 45
# main constituents, then the indice returned is 46, since constituents
# are counted from zero. ( 45 - 1 + 2 = 46 )
def vuf (inname):

    return 'need to finish this'

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

def main(cfgfilePath, inputfilePath, outputfilePath):
# def main(argv):

    # Variables
    global constDict, constFilePathsDict
    #constNames = []
    #constFiles = []
    #constFrames = []
    cplx = -9

    if cfgfilePath == '' or inputfilePath == '' or outputfilePath == '':
        print('To start tide correction, enter the config, input and output files.')
        sys.exit(2)

    # if len(argv) != 3:
    #     print('To start tide correction, enter the config, input and output files.')
    #     sys.exit(2)

    # cfgfilePath = argv[0]
    # inputfilePath = argv[1]
    # outputfilePath = argv[2]

    # Open the config file
    cfgfile = open(cfgfilePath, "r+")

    # Node, element and IOS files
    nodefile = cfgfile.readline().strip()
    elemfile = cfgfile.readline().strip()
    iosfile = cfgfile.readline().strip()
    #print(nodefile, elemfile, iosfile)

    # Num constituents, constituent names and their tidal files
    numConst = int(cfgfile.readline().strip())
    base_ext = ''
    for i in range(numConst):
        constName = cfgfile.readline().strip()
        filePath = cfgfile.readline().strip()
        constFilePathsDict[constName] = filePath
        print(constName, filePath)
        # Set base extension for comparison
        if i == 0:
            base_ext = filePath.split('.')[2].upper()
        else:
            # Check to make sure all extensions are the same for all tidal files
            ext = filePath.split('.')[2].upper()
            if ext != base_ext:
                print("Error. The extension (%s) for constituent data file (%s) is not the same as the first extension (%s).\nPlease check the filenames in the config file.", ext,filePath,base_ext)
                sys.exit(2)
            elif ext == "S2C":
                cplx = 0
            elif ext == "V2C":
                cplx = 1
            # TODO: add functionality for V2R extensions
            # elif ext == "V2R":
            #     continue
            else:
                print("Error. Unrecognized data file format (%s).\nPlease check the filenames in the config file.", ext)
                sys.exit(2)

    # Close opened file
    cfgfile.close()

    #----------------------------------------------------------
    # Read in the mesh (to numpy arrays)
    readNodes(nodefile, elemfile)

    # TODO: Append min tide to node in nodeMinTidesList - but do we need min_tide at all?
    #for node in nodesList: nodeMinTidesList.append(0.0)
    print("num nodes: ", numNodes)
    print("num elements: ", numElems)
    #----------------------------------------------------------
    # Load the model tidal data for each constituent into arrays
    for const in constFilePathsDict.keys():

        # Variables
        linesToSkip = 0
        infoHeader = ''
        indivConstDict = {}

        # Extract header info first
        tmp_file = open(constFilePathsDict[const],"r")
        tmp_line = ''
        while(True):
            tmp_line = tmp_file.readline()
            if (tmp_line.split(' ')[0] == '1'):
                break
            else:
                infoHeader += tmp_line.rstrip("\n")
            linesToSkip += 1
        #print('Header info: ', infoHeader)

        # This current line will be used in next stage
        tmp_line = tmp_line.rstrip("\n").split()

        # Make arrays for indivConstDict
        # TODO: add functionality for V2R extensions
        indivConstDict['amp'] = []
        indivConstDict['phase'] = []
        #indivConstDict['min_tide'] = []
        if cplx == 1: # V2C has extra columns
            indivConstDict['amp2'] = []
            indivConstDict['phase2'] = []
            #indivConstDict['min_tide2'] = []

        # Populate arrays with rest of data
        while(True):
            indivConstDict['amp'].append(float(tmp_line[1]))
            indivConstDict['phase'].append(float(tmp_line[2]))
            #indivConstDict['min_tide'].append(0.0)
            # Check if V2C file
            if cplx == 1:
                indivConstDict['amp2'].append(float(tmp_line[3]))
                indivConstDict['phase2'].append(float(tmp_line[4]))
                #indivConstDict['min_tide2'].append(0.0)
            # Get new line
            tmp_line = tmp_file.readline().rstrip("\n").split()
            if len(tmp_line) == 0:
                print("Done! length of arrays: ", len(indivConstDict['amp']))
                break

        tmp_file.close()

        # Populate and save new data frame ** ADD FOR V2R
        # TODO: Add for V2R
        constDict[const] = indivConstDict

    #----------------------------------------------------------
    # Min tides? doesn't seem to do anything in webtidecor
    #----------------------------------------------------------

    # Set initial values
    day = 0
    month = 0
    year = 0
    hour = 0
    minute = 0
    seconds = 0.0
    longitude = 0.0
    latitude = 0.0
    reslt = 0.0
    reslt2 = 0.0

    # Read in the constituent data from iosfile
    # TODO: change to arrays, not dataframes
    openvuf(iosfile)

    # Loop through the input file, for each line calculate the tidal correction
    # Fields: long, lat, year, dayofyear, hour, minute, seconds
    # Open the input/output files

    # open input and output files
    inputfile = open(inputfilePath,"r")
    outputfile = open(outputfilePath,"w")

    # Read in input values
    for line in inputfile:
        line = line.split()
        longitude = float(line[0])
        latitude = float(line[1])
        year = int(line[2])
        dayofyear = int(line[3])
        month, day = get_date(year, dayofyear)
        hour = int(line[4])
        minute = int(line[5])
        seconds = float(line[6])

        print(longitude, latitude, year, month, day, hour, minute, seconds)

        # Get element number containing the point
        ele, basis = basis2d( longitude, latitude )
        print("element: ", ele)
        print("basis: ", basis)

        if ele < 0:
            print("Error finding element! Some Markers not in the domain...\n")
            sys.exit(2)

        # If the point is inside and element, calculate the tidal correction for
        # each node of the element and interpolate to get the tidal correction
        # at the new position.
        elemRes, elemRes2 = [],[]
        # elemMin, elemMin2 = [],[]
        for i in range(3):
            nodeIdx = elemsList[ele][i] - 1
            # node numbers in nodesDf start at 1 (correct because arrays start at 0)
            #print("elem node#", i, " : ", node)
            elemRes.append(TideP( day, month, year, hour, minute, seconds, latitude, nodeIdx, 0))
            if (cplx > 0):
                elemRes2.append(TideP( day, month, year, hour, minute, seconds, latitude, nodeIdx, 1))
            #min_tide = constDict[constName]['amp'][nodeIdx]
            #elemMin.append(constDict)

        # Calculate results
        reslt = elemRes[0] * basis[0] + elemRes[1] * basis[1] + elemRes[2] * basis[2]
        print("result: ", reslt)
        # + elem_min[0] * basis[0] + elem_min[1] * basis[1] + elem_min[2] * basis[2]
        if cplx > 0:
            reslt2 = elemRes2[0] * basis[0] + elemRes2[1] * basis[1] + elemRes2[2] * basis[2]
            print("result2: ", reslt2)
            # + elem_min2[0] * basis[0] + elem_min2[1] * basis[1] + elem_min2[2] * basis[2]

        # Output the tidal correction into file
        if cplx == 0:
            str = '{0:.4f} {1:.8f} {2:.8f} {3} {4} {5} {6} {7:.2f}\n'.format(reslt, longitude, latitude, year, dayofyear, hour, minute, seconds)
            outputfile.write(str)
            print("res: ", str)
        else:
            str = '{0:.4f} {1:.4f} {2:.8f} {3:.8f} {4} {5} {6} {7} {8:.2f}\n'.format(reslt, reslt2, longitude, latitude, year, dayofyear, hour, minute, seconds)
            outputfile.write(str)
            print("res2: ", str)

    outputfile.close()
    inputfile.close()

    # FIN!

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

###########################################################################
# This will be run if the library is used "stand alone" on the command line
###########################################################################
# if __name__ == "__main__":
#
#     main(sys.argv[1:])



