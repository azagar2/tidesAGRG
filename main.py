import webtide_calculator as wt

# Sol file
basePath = '/Users/andrea/PycharmProjects/Tides3.6/inputData/'
# solPath = basePath + 't200_TC_MULTI_URAPID_C10_C11_C12.sol'
solPath = basePath + 't200_TC_MULTI_RAPID_BART_LTWN.sol'
# solPath = basePath + 't200_TC_MULTI_RAPID_SHEL.sol'
webPath = basePath + 'track.txt'
datumPath = basePath + 'MWLm_WGS84.txt'

# Webtide params
cfg = 'tidecor2.4.cfg'  # Config file
inp = solPath           # Input file, needs to be a sol file or webtide input text file as of now
out = 'results.csv'      # Can be named anything
hertz = 0.01

# Call webtide
wt.main(cfg, inp, out, datumPath, hertz)
