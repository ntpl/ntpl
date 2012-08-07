## ## ## This script will provide examples of how to implement
## ## ## the ntpy.param module
ntplPATH = '/home/kevin/projects/ntpl/Codes'

## Import the module
import sys
sys.path.append(ntplPATH) # Needed to link the ntplpy module
import ntpy.param.lj as lj
import ntpy.param.const as cn


## ## We will use the param module to easily call common parameters.
## ## Simply call the system and pass the value function a string.

## Here we extract the lat20 lattice constant for lj from value()
lat20 = (lj.value('lat20') / lj.value('sigma')) * 1e-10

print lat20

## We can also use the latfit() for extrapolate lattice constants
lat4 = (lj.latfit(4) / lj.value('sigma')) * 1e-10

print lat4

## The param module also supports common constants

print cn.value('kb')

print 12.6 * cn.value('tocenti')
