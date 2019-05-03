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
#  * Modifications and additions made by Andrea Zagar @ COGS (NSCC) - April, 2019
#  ******************************************************************************/
#  WEBTIDE HELP: http://www.bio.gc.ca/science/research-recherche/ocean/webtide/help-aider-en.php

import math
import sys
import os
import argparse
import csv
import numpy as np
from scipy.spatial import distance
from AGRG_time import time_gps2datetime
from read_nav_file import readSol
from memory_profiler import profile

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
version = 0             # Are the constituent files S2C or V2C? (S2C is 0, V2C is 1)
solFreq = 200.00        # Default frequency for sol file data
osPathChar = "/"        # Separator for file path depending on operating system (change to \ for windows)

# Define the value for 2 * pi
twopi = 2.0 * math.acos( -1.0 )

# This defines the (Julian) date of the switch to the Gregorian calendar (Used in the julday subroutine)
IGREG = (15 + 31 * (10 + 12 * 1582 ))

# These define the indices for the 5 constituents used in this program
# (Perhaps someday users will be allowed to choose which and how many to use)
m2 = 0
s2 = 1
n2 = 2
k1 = 3
o1 = 4

# Indices for constituent arrays
ampIdx = 0
phaseIdx = 1
amp2Idx = 2
phase2Idx = 3

# Table of number of days for each month
daytable = [
    [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
    [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]]

# Global variables:
global nodesArr, elemsArr, datumCoordArr, datumZArr
numNodes, numElems = 0,0
numSatAttr = 5

# Dictionaries for reading iosfile (stores only constituents listed in config in those dictionaries)
mainConstDict = {}      # -> key: main constituent name, value: main constituent factors
shallConstDict = {}     # -> key: shallow constituent name, value: shallow constituent factors

# Dictionary of main constituents used in calculations (around 5)
constDict = {}          # -> key: constituent name, value: s2c or v2c data for constituent
constFilePathsDict = {} # -> key: constituent name, value: filename of s2c or v2c data

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
# @profile
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
# Finds the closest node to a point (ptx, pty) and the element containing
# that point, if one exists.
# Also gets the basis functions for interpolations to that point.
# Returns the element number containing the point (ptx, pty).
# Returns -999 if no element found containing the point (ptx, pty ).
# @profile
def basis2d(ptx, pty):

    # Find 3 node IDs closest to point (node IDs are 1 greater than their indexes because node file starts at 1, not 0
    closeNodes = closestNodesID(ptx, pty)

    # Try to find an element that contains the point and the closest node
    ele = -999
    elemsList = []
    outsideMeshFlag = 0
    closestDist = 0

    for closeNodeID in closeNodes:
        # Read in node vals for each element (subtract 1 because index of nodesArr starts at 0 and not 1)
        for index in range(elemsArr.shape[0]):
            row = elemsArr[index]
            id1, id2, id3 = row[0], row[1], row[2]
            xCoords, yCoords = [],[]
            # If one or more of the nodes is in this element
            if (( closeNodeID == (id1) ) or ( closeNodeID == (id2) ) or ( closeNodeID == (id3) )):
                # Check if element already added to list
                if index not in [row[1] for row in elemsList]:
                    elemsList.append([closeNodeID, index])
                # Store lat and long of each node (need to subtract 1 from ID to get index in nodesArr)
                xCoords.extend((nodesArr[id1-1][0], nodesArr[id2-1][0], nodesArr[id3-1][0]))
                yCoords.extend((nodesArr[id1-1][1], nodesArr[id2-1][1], nodesArr[id3-1][1]))
                # See if the point is within this element
                flag = raybound(xCoords, yCoords, ptx, pty)
                if flag == 1: # The point is within the element
                    ele = index
                    basis = phi2d(xCoords, yCoords, ptx, pty)
                    return ele, basis, closestDist, ptx, pty, outsideMeshFlag

    # If the closest node's elements don't work, search through all elements
    #print("searching through all elements")
    # for index in range(elemsArr.shape[0]):
    #     row = elemsArr[index]
    #     n11, n22, n33 = row[0]-1, row[1]-1, row[2]-1
    #     xCoordsFull, yCoordsFull = [],[]
    #     # Store lat and long of each node
    #     xCoordsFull.extend((nodesArr[n11][0], nodesArr[n22][0], nodesArr[n33][0]))
    #     yCoordsFull.extend((nodesArr[n11][1], nodesArr[n22][1], nodesArr[n33][1]))
    #     # See if the point is within this element
    #     flag = raybound(xCoordsFull, yCoordsFull, ptx, pty)
    #     if flag == 1: # The point is within the element
    #         # print("found index in full set: ", index)
    #         # print("n11: ",n11,", n22: ", n22,", n33: ", n33)
    #         ele = index
    #         # print("element-2: ", ele)
    #         basis = phi2d(xCoordsFull, yCoordsFull, ptx, pty)
    #         return closestNodeID, ele, basis, closestDist, ptx, pty, outsideMeshFlag

    # If not inside closest node's elements, find closest element to point using those elements (elemsList)
    closestDist = 9999
    closestCentX, closestCentY = 0,0
    point = np.array([[ptx,pty]])
    # For each element containing the closest node
    for closeNodeID, index in elemsList:
        row = elemsArr[index]
        # print("Nodes in element ", index, ": ", row, " and closest ID: ", closeNodeID)
        xCoordsClose, yCoordsClose = [],[]
        id1, id2, id3 = row[0], row[1], row[2]
        # For each element, find centroid of triangle and then find distance from point to centroid
        centX = (nodesArr[id1-1][0]+nodesArr[id2-1][0]+nodesArr[id3-1][0])/3
        centY = (nodesArr[id1-1][1]+nodesArr[id2-1][1]+nodesArr[id3-1][1])/3
        # Get Manhattan distance from point to centroid and compare
        centroid = np.array([[centX,centY]])
        dist = distance.cdist(point, centroid,'cityblock')
        # print("Dist from point: ", point, " .. to centroid: ", centroid, "----> ", dist)
        # If closest distance so far, update ele and basis
        if dist < closestDist:
            ele = index
            closestDist = dist
            outsideMeshFlag = 1
            closestCentX, closestCentY = centX, centY
            # Add node coordinates to list
            xCoordsClose.extend((nodesArr[id1-1][0], nodesArr[id2-1][0], nodesArr[id3-1][0]))
            yCoordsClose.extend((nodesArr[id1-1][1], nodesArr[id2-1][1], nodesArr[id3-1][1]))
            basis = phi2d(xCoordsClose, yCoordsClose, ptx, pty)

    # Closest dist, node is outside of mesh
    closestDist = closestDist[0][0]

    # This shouldn't ever happen, but a check in case
    if ele < 0:
        print("Error finding element! Some Markers not in the domain...\n")
        sys.exit(2)

    return ele, basis, closestDist, closestCentX, closestCentY, outsideMeshFlag

#------------------------------------------------------------------------------------
# Calculate the tidal correction for each node of the element and interpolate
# to get the tidal correction at the new position.
# @profile
def calculateTide(ele, basis, day, month, year, hour, minute, seconds, latitude):

    # Variables
    global version
    elemRes, elemRes2, result = [],[],[]

    for i in range(3):
        # Node numbers in nodesArr start at 1 (so need to correct because python arrays start at 0)
        nodeIdx = elemsArr[ele][i] - 1
        # Calculations
        elemRes.append(TideP( day, month, year, hour, minute, seconds, latitude, nodeIdx, 0))
        if (version > 0):
            elemRes2.append(TideP( day, month, year, hour, minute, seconds, latitude, nodeIdx, 1))

    # Calculate results
    result.append(elemRes[0] * basis[0] + elemRes[1] * basis[1] + elemRes[2] * basis[2])
    if version > 0: # Two results returned for V2C file
        result.append(elemRes2[0] * basis[0] + elemRes2[1] * basis[1] + elemRes2[2] * basis[2])

    return result

#------------------------------------------------------------------------------------
# Find the indexes of the 3 nodes that are closest to the point (ptx, pty).
# @profile
def closestNodesID(ptx, pty):

    # Add node IDs to this temporary array to make it easier to find smallest 3 distances
    # arrSize = len(nodesArrCopy)
    # nodeIds = np.arange(1,arrSize+1).reshape(arrSize, 1)
    # tempNodesArr = np.append(nodeIds, nodesArrCopy, axis=1)

    # Use Manhattan (cityblock) distance and find closest node to point, then return index of node
    point = np.array([[ptx,pty]])
    closeIdx1 = distance.cdist(point, nodesArr,'cityblock').argmin()
    closeRow1 = np.copy(nodesArr[closeIdx1])
    nodesArr[closeIdx1] = [999,999]
    # Find second closest point
    closeIdx2 = distance.cdist(point, nodesArr,'cityblock').argmin()
    closeRow2 = np.copy(nodesArr[closeIdx2])
    nodesArr[closeIdx2] = [999,999]
    # Find third closest point
    closeIdx3 = distance.cdist(point, nodesArr,'cityblock').argmin()

    # Replace rows with original values
    nodesArr[closeIdx1] = closeRow1
    nodesArr[closeIdx2] = closeRow2

    # Cdist.argmin returns index of closest node in nodesArr
    # Need to return closest node ID, so have to add 1 to index (nodesArr is 0 indexed, and node IDs are 1 indexed)
    return [closeIdx1+1, closeIdx2+1, closeIdx3+1]

#------------------------------------------------------------------------------------
# Find the z-value of the datum grid point closest to the point (ptx, pty).
def getDatumTransform(ptx, pty):

    # Find index of closest distance
    point = np.array([[ptx,pty]])
    closeIndex = distance.cdist(point, datumCoordArr,'cityblock').argmin()
    return datumZArr[closeIndex]

#------------------------------------------------------------------------------------
# Get the day and month from the day # of the year
# @profile
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
# @profile
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
# @profile
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
        if line[0] not in constDict.keys():
            continue
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
                    newMain.sats.append(newSat)

                satCount -= numLineSats

        # Save main constituent into dictionary
        mainConstDict[newMain.name] = newMain

    # Get shallow constituents
    while (True):
        # Get line from file
        line = iosfile.readline().rstrip("/n").split()
        # Check for blank line to signal end of main constituents
        lineLen = len(line)
        if lineLen <= 1:
            break
        if line[0] not in constDict.keys():
            continue
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

    # Close file
    iosfile.close()

#------------------------------------------------------------------------------------
#Calculates the basis functions for interpolating to a point inside an element.
# @profile
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
# @profile
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
    return bcross % 2

#------------------------------------------------------------------------------------
# Loop through the SOL input file, for each point calculate the tidal correction
# Write results to an output csv file
# @profile
def readInputSol(inputfilePath, outputfilePath, hertz):

    # Read in sol file
    rowsToSkip = int(solFreq/hertz)
    solArr = readSol(inputfilePath)[0::rowsToSkip]

    # Modify sol array before using
    colsToKeep = ["lat", "lon", "alt",
        "gps_week_number","time",
        "standard_deviation_latitude", "standard_deviation_longitude", "standard_deviation_height",
        "roll", "pitch", "heading",
        "standard_deviation_roll", "standard_deviation_pitch", "standard_deviation_true_heading",
        "nsspeed", "ewspeed", "vertspeed",
        "standard_deviation_north_velocity", "standard_deviation_east_velocity", "standard_deviation_up_velocity"]
    solArr['lon'] = np.rad2deg(solArr['lon'])
    solArr['lat'] = np.rad2deg(solArr['lat'])
    solArr = solArr[colsToKeep]

    # Make new 2d list to hold new values (convert to numpy array at end, then join with original solArr)
    # New values: result, result2, ellipsoidVal, flag, centroidX, centroidY, closestDistance
    newValsArr = []

    # Read in rows of input
    for line in solArr:

        # Variables needed for tidal calculation functions
        longitude = float(line['lon'])
        latitude = float(line['lat'])
        datetime = time_gps2datetime(line['gps_week_number'].item(),line['time'].item())
        year = datetime.date().year
        month = datetime.date().month
        day = datetime.date().day
        hour = datetime.time().hour
        minute = datetime.time().minute
        seconds = datetime.time().second

        # Get element containing the point, basis function, distance to centroid, centroid, flag (point outside of mesh = 1, inside = 0)
        ele, basis, closestDist, centroidX, centroidY, flag = basis2d(longitude, latitude)

        # If input point is inside an element in the mesh, then set centroid equal to that point
        if not flag:
            centroidX = longitude
            centroidY = latitude

        # Calculate results
        results = calculateTide(ele, basis, day, month, year, hour, minute, seconds, latitude)

        # Get datum transformation if necessary, then insert values into temporary new array
        if version == 0:
            # Get datum transformation based on centroid point, and difference between that and tide elevation
            ellipsoidVal = getDatumTransform(centroidX, centroidY)
            ellipsoidDiff = ellipsoidVal - results[0]
            newValsArr.append((results[0], ellipsoidVal, ellipsoidDiff, flag, centroidX, centroidY, closestDist))
        else:
            newValsArr.append((results[0], results[1], flag, centroidX, centroidY, closestDist))

    # Create numpy array based on 2d list
    if version == 0:
        newDtypes = [('tideElevMWL', 'float32'), ('WGS84Elev', 'float32'), ('MWL_WGS84_ElevDiff', 'float32'), ('flag', 'int8'), ('centroidLong', 'float32'), ('centroidLat', 'float32'), ('distToCentroid', 'float32')]
    else:
        newDtypes = [('tideCurrRes1', np.float32), ('tideCurrRes2', np.float32), ('flag', np.int8), ('centroidLong', np.float32), ('centroidLat', np.float32), ('distToCentroid', np.float32)]
    newValsNp = np.array(newValsArr, dtype=newDtypes)

    # Recreate solArr dtype list in order to join list with newDtypes
    solDtypes = []
    for colName in colsToKeep:
        solDtypes.append((colName, solArr.dtype[colName].name))

    # Extend original solArr to include new array of values
    new_dt = np.dtype(solDtypes + newDtypes)
    solArrFinal = np.empty(solArr.shape, dtype=new_dt)

    # Populate final array
    for name in colsToKeep:
        solArrFinal[name] = solArr[name]

    for name in [i[0] for i in newDtypes]:
        solArrFinal[name] = newValsNp[name]

    # Header and formatting
    header = solArrFinal.dtype.names
    header = ",".join(header)
    # "lat", "lon", "alt", "gps_week_number","time", "standard_deviation_latitude", "standard_deviation_longitude",
    # "standard_deviation_height", "roll", "pitch", "heading", "standard_deviation_roll", "standard_deviation_pitch",
    # "standard_deviation_true_heading", "nsspeed", "ewspeed", "vertspeed", "standard_deviation_north_velocity",
    # "standard_deviation_east_velocity", "standard_deviation_up_velocity",
    if version == 0:
        # "tideElevMWL", "WGS84Elev", "MWL_WGS84_ElevDiff", "flag", "centroidLong", "centroidLat", "distToCentroid"
        format = "%1.8f, %1.8f, %1.8f, %d, %1.8f, %1.8f, %1.8f," \
                 "%1.8f, %1.8f, %1.8f, %1.8f, %1.8f, %1.8f," \
                 "%1.8f, %1.8f, %1.8f, %1.8f, %1.8f," \
                 "%1.8f, %1.8f, " \
                 "%1.6f, %1.6f, %1.6f, %d, %1.8f, %1.8f, %1.8f"
    else:
        # "tideCurrRes1", "tideCurrRes2", "flag", "centroidLong", "centroidLat", "distToCentroid"
        format = "%1.8f, %1.8f, %1.8f, %d, %1.8f, %1.8f, %1.8f," \
                 "%1.8f, %1.8f, %1.8f, %1.8f, %1.8f, %1.8f," \
                 "%1.8f, %1.8f, %1.8f, %1.8f, %1.8f," \
                 "%1.8f, %1.8f, " \
                 "%1.6f, %1.6f, %d, %1.8f, %1.8f, %1.8f"
    # # Save array to CSV file
    np.savetxt(outputfilePath, solArrFinal, delimiter=",", header=header, fmt=format, comments='')

#------------------------------------------------------------------------------------
# Loop through the input file, for each point calculate the tidal correction
# Write results to an output csv file
# @profile
def readInputText(inputfilePath, outputfilePath):

    # Open input file
    inputfile = open(inputfilePath,"r")

    # Open output file for writing (csv)
    outputfile = open(outputfilePath, 'w', newline='')
    writer = csv.writer(outputfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    # Add headers to csv file
    if version == 0:
        writer.writerow(["tideElevMWL","WGS84Elev","MWL_WGS84_ElevDiff","longitude","latitude","year","month","day","hour","minute","second","flag","centroidLong","centroidLat","distToCentroid"])
    else:
        writer.writerow(["tideCurrent1","tideCurrent2","longitude","latitude","year","month","day","hour","minute","second","flag","centroidLong","centroidLat","distToCentroid"])

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

        # Get element number containing the point
        ele, basis, closestDist, centroidX, centroidY, flag = basis2d(longitude, latitude)

        if not flag:
            centroidX = longitude
            centroidY = latitude

        # Calculate results
        results = calculateTide(ele, basis, day, month, year, hour, minute, seconds, latitude)

        # Get datum transformation if necessary, then output the tidal correction into file
        if version == 0:
            # Get datum transformation based on centroid point, and difference between that and tide elevation
            ellipsoidVal = getDatumTransform(centroidX, centroidY)
            ellipsoidDiff = ellipsoidVal - results[0]
            str = '{0:.6f},{1:.6f},{2:.6f},{3:.8f},{4:.8f},{5},{6},{7},{8},{9},{10:.5f},{11},{12:.6f},{13:.8f},{14:.8f}'.format(results[0], ellipsoidVal, ellipsoidDiff, longitude, latitude, year, month, day, hour, minute, seconds, flag, centroidX, centroidY, closestDist)
        else:
            str = '{0:.6f},{1:.6f},{2:.8f},{3:.8f},{4},{5},{6},{7},{8},{9:.5f},{10},{11:.6f},{12:.8f},{13:.8f}'.format(results[0], results[1], longitude, latitude, year, month, day, hour, minute, seconds, flag, centroidX, centroidY, closestDist)
        arr = str.split(",")
        writer.writerow(arr)

    # Close files
    outputfile.close()
    inputfile.close()

#------------------------------------------------------------------------------------
# Read in the mesh from a .nod/.ele set of files.
# @profile
def readDatum(datumFile):

    # Reference (ID) starts at 0 in text file
    global datumCoordArr, datumZArr

    # Read datum file into numpy array (ONLY lat and long)
    datumCoordArr = np.loadtxt(datumFile, dtype=float, usecols=[1,2,3],delimiter=",",skiprows=1)
    datumZArr = np.copy(datumCoordArr[:,2])
    # Remove Z column from coordinate array
    datumCoordArr = datumCoordArr[:,[1,0]]

    return

#------------------------------------------------------------------------------------
# Read in the mesh from a .nod/.ele set of files.
# @profile
def readNodes(nodeFile, elementFile):

    # These frames have to have indexes corresponding to element file indexes (starting at 1)
    global nodesArr, elemsArr, numNodes, numElems

    # Read node file into numpy array
    nodesArr = np.loadtxt(nodeFile, dtype=float, usecols=[1,2])

    # Read element into numpy array
    elemsArr = np.loadtxt(elementFile, dtype=int, usecols=[1,2,3])

    # Number of nodes and elements
    numNodes = nodesArr.shape[0]
    numElems = elemsArr.shape[0]

    return

#------------------------------------------------------------------------------------
# Calculate the amplitudes, phases, etc. for each of the constituents
# @profile
def setvuf(kh, xlat):

    # Global and local variables
    global mainConstDict, shallConstDict
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
# @profile
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
#------------------------------------------------------------------------------------
# @profile
def main(cfgfilePath, inputfilePath, outputfilePath, datumFilePath, hertz=solFreq, baseHertz=solFreq):

    # Variables
    global constDict, constFilePathsDict
    version = -9

    if cfgfilePath == '' or inputfilePath == '' or outputfilePath == '' or datumFilePath == '':
        print('To start tide correction, enter the config, input and output files.')
        sys.exit(2)

    # Open the config file
    cfgfile = open(cfgfilePath, "r+")

    # Node, element and IOS files
    nodefile = cfgfile.readline().strip()
    elemfile = cfgfile.readline().strip()
    iosfile = cfgfile.readline().strip()

    # Num constituents, constituent names and their tidal files
    numConst = int(cfgfile.readline().strip())
    base_ext = ''
    for i in range(numConst):
        constName = cfgfile.readline().strip()
        filePath = cfgfile.readline().strip()
        constFilePathsDict[constName] = filePath
        # Set base extension for comparison
        if i == 0:
            base_ext = filePath.split('.')[2].upper()
        else:
            # Check to make sure all extensions are the same for all tidal files
            path = filePath.split(osPathChar)[-1]
            ext = path.split('.')[2].upper()
            if ext != base_ext:
                print("Error. The extension (%s) for constituent data file (%s) is not the same as the first extension (%s).\nPlease check the filenames in the config file.", ext,filePath,base_ext)
                sys.exit(2)
            elif ext == "S2C":
                version = 0
            elif ext == "V2C":
                version = 1
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

    # Read in the datum transformation file
    readDatum(datumFilePath)

    # TODO: Append min tide to node in nodeMinTidesList - but do we need min_tide at all?
    #for node in nodesList: nodeMinTidesList.append(0.0)
    #----------------------------------------------------------
    # Load the model tidal data for each constituent into arrays
    for const in constFilePathsDict.keys():

        # Variables
        linesToSkip = 0
        infoHeader = ''

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
        tmp_file.close()

        # TODO: add functionality for V2R extensions
        # TODO: add functionality to save min_tides in numpy arrays (right now values are all initialized to 0.0)
        # min_tide and min_tide2 (for V2C)

        # Read rest of file into numpy array
        # If V2C, has 2 more cols
        if version == 1:
            dts = [('id',int),('amp',float),('phase',float),('amp2',float),('phase2',float)]
        else:
            dts = [('id',int),('amp',float),('phase',float)]

        constArr = np.loadtxt(constFilePathsDict[const], dtype=dts, skiprows=linesToSkip)

        # Populate and save new array in dictionary of constituent data
        constDict[const] = constArr

    #----------------------------------------------------------
    # TODO: Min tides? doesn't seem to do anything in webtidecor
    #----------------------------------------------------------

    # Read in the constituent data from iosfile
    openvuf(iosfile)

    #----------------------------------------------------------
    # Read input file
    #----------------------------------------------------------

    # Figure out the input type based on file extension (.sol or .txt)
    fileExt = inputfilePath.split(".")[-1]

    if fileExt == 'sol': # Read sol file
        readInputSol(inputfilePath, outputfilePath, hertz)
    elif fileExt == 'txt': # Regular webtide input file
        readInputText(inputfilePath, outputfilePath)
    else:
        print("Invalid input file type. Exiting")
        sys.exit()

    # FIN!

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

###########################################################################
# This will be run if the script is used on the command line
###########################################################################

if __name__=='__main__':
    #Get the input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--config','-c',help ='Input webtide configuration (.cfg) file to read',metavar="<file>",required=True)
    parser.add_argument('--input','-i',help ='Input navigation (sol/txt) file to read',metavar="<file>",required=True)
    parser.add_argument('--output','-o',help ='Output CSV file to write',default="results.csv",metavar="<csvfile>")
    parser.add_argument('--datum','-d',help ='Input datum conversion (txt) file to read',metavar="<file>",required=True)
    parser.add_argument('--hertz','-hz',help ='Value in hertz to sample sol file',default=solFreq,metavar="<float>")
    parser.add_argument('--hertzdef','-hzdef',help ='Default value in hertz of sol file data',default=solFreq,metavar="<float>")

    commandline=parser.parse_args()

    # Check if the config file exists - exit if not
    if not os.path.exists(commandline.config):
        print("%s the config file cannot be found - are you sure it exists?"%commandline.config)
        sys.exit(1)

    # Check if the input navigation file exists - exit if not
    if not os.path.exists(commandline.input):
        print("%s the input file cannot be found - are you sure it exists?"%commandline.input)
        sys.exit(1)

    # Check if the datum conversion file exists - exit if not
    if not os.path.exists(commandline.datum):
        print("%s the datum conversion file cannot be found - are you sure it exists?"%commandline.datum)
        sys.exit(1)

    #read in the sol file - an arry with each row as a record
    try:
        print("Running program ...")
        main(commandline.config, commandline.input, commandline.output, commandline.datum, hertz=float(commandline.hertz), baseHertz=float(commandline.hertzdef))
    except Exception as e:
        print("Failed.\nSomething went wrong!")
        print(e)
