## ## ## This script will provide examples of how to implement
## ## ## the ntplpy.strchg module

## First import the module
## Import the module
import sys
sys.path.append('/home/kevin/ntpl') # Needed to link the ntplpy module
import ntplpy.strchg as st

## ## We will use the sed command. First, create a dict of key value pairs.
## ## Each key is an original string and each value is the new string.

##	First, create a key:value dict.
mydict = dict({'key': 'lock', 'string': 'block'})
# The following methods are equivalent ways.
# mydict = dict(TEMP_PARAM='80', STRAIN_PARAM='0.02')
# mydict = dict(zip(('TEMP_PARAM', 'STRAIN_PARAM'), ('80', '0.02')))
# mydict = dict([['TEMP_PARAM', '80'], ['STRAIN_PARAM', '0.02']])

## Then simply call the function with the original and new file names.
st.sed(mydict, 'strchgtest.py', 'newfile.test')
## ## And you're done!!!

## If you want to edit the file in-place, i.e. not create a new file,
## then simply make the new file name a blank string.
new = ''
st.sed(mydict, orig, new)

## You can also use relative or absolute pathways for the file names.
mydict = dict({'lock': 'key', 'block': 'string'})
orig = '/home/kevin/ntpl/strchgmod/newfile.test'
new = './strchgtest.py'
st.sed(mydict, orig, new)

