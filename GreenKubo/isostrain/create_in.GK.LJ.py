## ## ## Kevin Parrish - 5/16/12 - in.GK.LJ editor ## ## ##
import os

temp = [20, 50, 80, 100]
strain = [-0.12, -0.10, -0.08, -0.06, -0.04, -0.02, -0.01, 0, 0.05, 0.01, 0.015, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12]

numSeeds = 10
numTemps = len(temp)
numStrains = len(strain)


for istrain in range(numStrains):
	for itemp in range(numTemps):
		for iseed in range(numSeeds):
			## ## ## Change data in structures
			orig = 'SEED'
			new  = str((iseed + 2) * 11111)
			cmd1 = '-e \'s/' + orig + '/' + new + '/g\' '
	
			## ## ## Change temp
			orig = 'TEMP_PARAM'
			new  = str(temp[itemp])
			cmd2 = '-e \'s/' + orig + '/' + new + '/g\' '
	
			## ## ## Change strain 
			orig = 'STRAIN_PARAM'
			new  = str(strain[istrain])
			cmd3 = '-e \'s/' + orig + '/' + new + '/g\' '
	
			## ## ## Change neighbor list for liquid
			if temp[itemp] == 100:
				orig = '^#LIQUID_PARAM'
				new  = str('neigh_modify every 1 delay 0')
				cmd4 = '-e \'s/' + orig + '/' + new + '/g\' '
			else:
				cmd4 = ''
	
			## ## ## Change input file
			orig = 'IN_DATA'
			new  = 'in.data.' + str(itemp)
			cmd5 = '-e \'s/' + orig + '/' + new + '/g\' '
		
			## ## ## Change output file
			orig = 'OUT_CORR'
			new  = 'out.lammps.corr.' + str(iseed) + '.' + str(itemp) + '.' + str(istrain)
			cmd6 = '-e \'s/' + orig + '/' + new + '/g\' '

			## ## ## Change output file
			orig = 'OUT_STRAIN'
			new  = 'out.lammps.strain.' + str(iseed) + '.' + str(itemp) + '.' + str(istrain)
			cmd7 = '-e \'s/' + orig + '/' + new + '/g\' '

			## ## ## Change volume file
			orig = 'OUT_VOLUME'
			new  = 'out.lammps.vol.' + str(iseed) + '.' + str(itemp) + '.' + str(istrain)
			cmd8 = '-e \'s/' + orig + '/' + new + '/g\' '
		
			## ## ## Change position file
			orig = 'OUT_POS'
			new  = 'out.lammps.pos.' + str(iseed) + '.' + str(itemp) + '.' + str(istrain) + '.xyz'
			cmd9 = '-e \'s/' + orig + '/' + new + '/g\' '
		
			## ## ## Change main log file
			orig = 'MAIN_LOG_FILE'
			new  = 'out.lammps.main.' + str(iseed) + '.' + str(itemp) + '.' + str(istrain)
			cmd10 = '-e \'s/' + orig + '/' + new + '/g\' '
		
			## ## ## Change kappa log file
			orig = 'KAPPA_LOG_FILE'
			new  = 'out.lammps.kappa.' + str(iseed) + '.' + str(itemp) + '.' + str(istrain)
			cmd11 = '-e \'s/' + orig + '/' + new + '/g\' '
		
			## ## ## Execute Command
			os.system('sed '+ cmd1+ cmd2+ cmd3+ cmd4+ cmd5+ cmd6+ cmd7+ cmd8+ cmd9+ cmd10+ cmd11+ 'in.lammps.temp > ./in.lammps.'+ str(iseed)+ '.'+ str(itemp)+ '.'+ str(istrain))

			#sed -e 's/\<POS_IN_STRUCTURE\>/in_structure_1.pos/g' -e 's/\<OUT_CORR\>/out.GK.LJ.corr.1/g' in.GK.LJ > ./in.GK.LJ.1
## ## ## END FILE
