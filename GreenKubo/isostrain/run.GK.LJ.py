## ## ## Kevin Parrish - 5/16/12 - runs in.GK.LJ.1-5.1-3 ## ## ##
import os

for istrain in range(17):
	for itemp in range(4):
		for iseed in range(10):
			## ## ## Change data in structures
			orig = 'PROC_NAME'
			new  = 'KDP.lammps.' + str(iseed) + '.' + str(itemp) + '.' + str(istrain)
			cmd1 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
		
			## ## ## Change output file
			orig = 'SHELL_NAME'
			new  = 'in.lammps.' + str(iseed) + '.' + str(itemp) + '.' + str(istrain)
			cmd2 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
		
#	print 'sed ' + cmd1 + cmd2 + 'lmp_GK_KDP.sh > ./lmp_GK_KDP_' + str(a + 1) + '.sh' #for debugging
		
			## ## ## Execute sed command
			os.system('sed ' + cmd1 + cmd2 + 'lammps_.sh > ./lammps_' + str(iseed) + '_' + str(itemp) + '_' + str(istrain) + '.sh')
			#sed -e 's/\<POS_IN_STRUCTURE\>/in_structure_1.pos/g' -e 's/\<OUT_CORR\>/out.GK.LJ.corr.1/g' lmp_GK_KDP.sh > ./lmp_GK_KDP_1.sh
			
			## ## ## Execute qsub shell script
			os.system('qsub -p -1000 lammps_' + str(iseed) + '_' + str(itemp) + '_' + str(istrain) + '.sh')
## ## ## END FILE
