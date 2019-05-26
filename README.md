# AGRG Tide Project

This project is based on the existing application known as WebTide, which was originally created by the Bedford Institute of Oceanography (BIO) in Nova Scotia. The original code, written in C, was translated to Python 3.6 and then customized for AGRG's use.

**WebTide**: [http://www.bio.gc.ca/science/research-recherche/ocean/webtide/index-en.php](http://www.bio.gc.ca/science/research-recherche/ocean/webtide/index-en.php)
**Datasets**: Northwest Atlantic and Scotia-Fundy-Maine.

**Author**: Andrea Zagar (Centre of Geographic Sciences - NSCC, May 2019)

# Project Configuration


Consists of config files, input data, and core python files.

## Required External Python Libraries
- argparse
- csv
- numpy
- scipy.spatial
- memory_profiler (optional)


## Config Files
- tidecor2.4-nwatl-s2c.cfg
- config-nwatl
	- IOS_tidetbl
	- nwatl_ll.nod
	- nwatl.ele
	- K1.barotropic.s2c
	- K1.barotropic.v2c
	- { other constituent files ... }

|Filename  |Purpose |Contents  |
|----------|--------|----------|
|`IOS_tidetbl`|Contains overall constituent data|First section is *main constituents* in this format: name, dood numbers (6 ints in total), phase (float), number of satellites (int), satellite info for each satellite (if applicable). Second section is *shallow constituents* in this format: name, number of main constituent factors (int), main constituent factors (float and name)|
|`nwatl_ll.nod`|Contains node data|node ID, longitude, latitude|
|`nwatl.ele`|Contains element data|element ID, node 1, node 2, node 3 (elements are made up of 3 nodes, i.e. a triangle)|
|`K1.barotropic.s2c`|Contains data for specific constituent (for elevation calculations)|ID, amplitude, phase|
|`K1.barotropic.v2c`|Contains data for specific constituent (for current calculations)|ID, amplitude 1, phase 1, amplitude 2, phase 2|

### Master Configuration Text File Example

**File name**: tidecor2.4-nwatl-s2c.cfg

**Contents**:
config-nwatl/nwatl_ll.nod
config-nwatl/nwatl.ele
config-nwatl/IOS_tidetbl
5
M2
config-nwatl/M2.barotropic.s2c
S2
config-nwatl/S2.barotropic.s2c
N2
config-nwatl/N2.barotropic.s2c
K1
config-nwatl/K1.barotropic.s2c
O1
config-nwatl/O1.barotropic.s2c

## Input Data
The program can run using either a .sol file (for AGRG data) or a .txt file (general WebTide input format) for the input data. The datum transformation file is required to put the MWL measurements into ellipsoidal heights.

- inputData
	- MWLm_WGS84.txt
	- t200_TC_MULTI_RAPID_SHEL.sol
	- track.txt

|Filename  |Purpose |Contents  |
|----------|--------|----------|
|`MWLm_WGS84.txt`|Datum transformation file|Contains nodes for a mesh that represents elevation difference between two datums in metres (in this case, MWL and WGS84)|
|`.sol`|AGRG flight data file |see read_nav_file.py|
|`track.txt`|Example of a general WebTide input data file|Longitude, latitude, year, day of year, hour, minute, seconds|

## Core Python Files
The program can run using either a .sol file (for AGRG data) or a .txt file (general WebTide input format) for the input data. The datum transformation file is required to put the MWL measurements into WGS84 ellipsoidal heights.

- AGRG_time.py
- read_nav_file.py
- main.py
- webtide_calculator.py

|Filename  |Purpose |
|----------|--------|
|`AGRG_time.py`|Convert gps time to datetime (time_gps2datetime function)|
|`read_nav_file.py`|Read .sol files and convert results to numpy arrays (readSol function)|
|`main.py`|Runs program from within IDE (i.e. PyCharm)|
|`webtide_calculator.py`|Full WebTide script, can be run from command line|

## Script Output

Results are written to a .csv file that is generated in the same directory as the core python scripts.


# Usage Instructions
This script can be run in an IDE or by using command line.

## Command Line

|Argument  |Short form| Purpose |Required|
|----------|--------|--------|--------|
|`-config`|`-c`|Path to WebTide configuration (.cfg) file to be read|Yes|
|`-input`|`-i`|Path to input (.sol/.txt) file to be read|Yes|
|`-output`|`-o`|Name of output CSV file|No (default is results.csv)|
|`-datum`|`-d`|Path to datum conversion (.txt) file to read|Yes|
|`-hertz`|`-hz`|Value in hertz to sample .sol file|No (default specified in script)|
|`-hertzdef`|`-hzdef`|Default value in hertz of .sol file data|No (default specified in script)|

**Note**: The default value for `-hz` and `-hzdef` in the script is 200. It is only used for .sol input files and is ignored for .txt input files.

### Example (UNIX)
Navigate to project directory and then run the following commands.

**Sol file**:
$ python3 webtide_calculator.py -c tidecor2.4-nwatl-s2c.cfg -i inputData/t200_TC_MULTI_RAPID_SHEL.sol -d inputData/MWLm_WGS84.txt -hz 0.01666 -o resultsShel.csv

**Text file**:
$ python3 webtide_calculator.py -c tidecor2.4-nwatl-s2c.cfg -i inputData/track.txt -d inputData/MWLm_WGS84.txt -o resultsTxt.csv

## IDE
Change any settings in `main.py` and then run it to execute WebTide script. Settings include the paths to the input and datum transformation files, as well as the WebTide parameters (config file name, input file, name of output csv file, desired frequency, and default frequency of sol data file).


# WebTide Logic Summary
- **main**
- Check to make sure input, config and output file paths are valid
- Open config file
	- Read in constituent file paths
	- Make sure constituent file extensions are the same
- **readNodes** Read in nodes and elements from .nod/.ele files
- **readDatum** Read in datum mesh from datum transformation file
- Load the model tidal data for each constituent into arrays
- **openvuf** Read in the overall constituent data from IOS_tidetbl
- Read input file
	- **readInputSol** Read a .sol file as input
		- Modify numpy array columns before reading in data
		- Read in data line by line
			-  Get longitude, latitude, year, month, day, hour, minute, seconds (**time_gps2datetime**)
			- **basis2d** Get element containing the point, basis function, distance to centroid, centroid, and flag
				- **closestNodesID** Find 3 node IDs closest to point
				- **raybound, phi2d** Try to find an element that contains the point and the closest node
			- **calculateTide** Calculate the tidal correction for each node of the element and interpolate to get the tidal correction at the new position
				- **TideP** Calculates and returns the tidal correction
					-  **julday** Get day of year
					-  **setvuf, astro_angles** Calculate the amplitudes, phases, etc. for each of the constituents
					- Calculate phase and amplitude, and final results
			- **getDatumTransform** Get datum transformation based on centroid point
			- Create numpy array based on calculation results and original .sol file data
			- Save array to CSV file
	- **readInputText** Read a .txt file as input
		- Open input file and output file for reading/writing
		- Add headers to output csv before reading in data
		- Read in data line by line
			- Get longitude, latitude, year, month, day, hour, minute, seconds (**get_date**)
			- *Same steps as **readInputSol** from **basis2d** up to and including **getDatumTransform***
			- Write to output csv file

# Future Work
- Need to test more on Windows environment (was written and tested on OS X), especially command line
- Test other WebTide datasets
- Improve speed and efficiency of script
- Several things in original WebTide code that were not added, including min_tide in calculations (which did not do anything) and V2R file functionality (see TODOs in webtide_calculator.py)
