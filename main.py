from read_nav_file import readSol
import webtide_calculator as wt

wt.main('tidecor2.4.cfg','track.txt','output.txt')

# arr = readSol('/Users/andrea/PycharmProjects/Tides2.7/t200_TC_MULTI_URAPID_C10_C11_C12.sol')
# print(arr[0])
# print(len(arr))

# Call webtide cor
cfg = 'tidecor2.4.cfg'
inp = 'track.txt'
out = 'output.txt'

#wt.main(cfg, inp, out)
