## ## ## This script will provide examples of how to implement
## ## ## the ntpy.strchg module

ntplPATH = '/home/kevin/projects/ntpl/Codes'

## Import the module
import os
import sys
sys.path.append(ntplPATH) # Needed to link the ntplpy module
import ntpy.strchg as st

## ## We will use the sed command. First, create a dict of key value pairs.
## ## Each key is an original string and each value is the new string.

##	First, create a key:value dict.
mydict = dict({'sed': 'MARIO', 'dict': 'LUIGI'})
# The following methods are equivalent ways.
# mydict = dict(TEMP_PARAM='80', STRAIN_PARAM='0.02')
# mydict = dict(zip(('TEMP_PARAM', 'STRAIN_PARAM'), ('80', '0.02')))
# mydict = dict([['TEMP_PARAM', '80'], ['STRAIN_PARAM', '0.02']])

## Then simply call the function with the original and new file names.
st.sed(mydict, 'strchgtest.py', 'tmp.newfile')
### And you're done!!!

## If you want to edit the file in-place, i.e. not create a new file,
## then simply make the new file name a blank string.
st.sed(mydict, 'tmp.newfile', '')

## You can also use relative or absolute pathways for the file names.
mydict = dict({'sed': 'MARIO', 'dict': 'LUIGI'})
orig = ntplPATH + '/strchgmod/strchgtest.py'
new = './tmp.newfile'
st.sed(mydict, orig, new)

## ## Now we will try the grep command. This grep implementation uses
## ## the grep -A command.

## Simply pass it the search string, number of lines to pull, the original
## file name and the new file name.
st.grep('Import the module', 4, 'strchgtest.py', 'tmp.newfile')

#cleanup
os.system('rm tmp.newfile')
