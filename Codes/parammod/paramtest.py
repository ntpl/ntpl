## ## ## This script will provide examples of how to implement
## ## ## the ntpy.param module
ntplPATH = '/home/kevin/projects/ntpl/Codes'

## Import the module
import sys
sys.path.append(ntplPATH) # Needed to link the ntplpy module
import ntpy.param as pm

## ## We will use the param module to easily call common parameters.
## ## Simply call the system and pass a string.

##	First, create a key:value dict.
lat20 = (pm.lj('lat20') / pm.lj('sigma')) * 1e-10

print lat20


