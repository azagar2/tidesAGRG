import webtide_calculator as wt
import os

# Directory
cwd = os.getcwd()
basePath = os.path.join(cwd, 'inputData')

# Sol files
solFile1 = 't200_TC_MULTI_URAPID_C10_C11_C12.sol'
solFile2 = 't200_TC_MULTI_RAPID_BART_LTWN.sol'
solFile3 = 't200_TC_MULTI_RAPID_SHEL.sol'
solPath = os.path.join(basePath, solFile3)

# Paths
txtPath = os.path.join(basePath, 'track.txt')
datumPath = os.path.join(basePath, 'MWLm_WGS84.txt')

# Webtide params
cfg = 'tidecor2.4-nwatl-s2c.cfg'  # Config file, can use ntwal or sshelf
inp = solPath           # Input file, needs to be a sol file or webtide input text file as of now
out = 'resultsTxtIDE.csv'     # Can be named anything
hertz = 0.01666         # One data point per minute
# baseHertz = 200       # Default frequency of the data itself

# Call webtide
wt.main(cfg, inp, out, datumPath, hertz)
# wt.main(cfg, inp, out, datumPath, hertz, baseHertz)
