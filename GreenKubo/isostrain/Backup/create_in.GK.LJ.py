## ## ## Kevin Parrish - 5/16/12 - in.GK.LJ editor ## ## ##
import os

temp = [20, 50, 80, 100]
strain = [-0.12, -0.10, -0.08, -0.06, -0.04, -0.02, -0.01, 0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12]

numSeeds = 10
numTemps = len(temp)
numStrains = len(strain)


for istrain in range(numStrains):
	for itemp in range(numTemps):
		for iseed in range(numSeeds):
			## ## ## Change data in structures
			orig = 'SEED'
			new  = str((iseed + 2) * 11111)
			cmd1 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
	
			## ## ## Change temp
			orig = 'TEMP_PARAM'
			new  = str(temp[itemp])
			cmd2 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
	
			## ## ## Change strain 
			orig = 'STRAIN_PARAM'
			new  = str(strain[istrain])
			cmd3 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
	
			## ## ## Change neighbor list for liquid
			if itemp == 3:
				orig = '#LIQUID_PARAM'
				new  = str(neigh_modify every 1 delay 0)
				cmd4 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
	
			## ## ## Change input file
			orig = 'IN_DATA'
			new  = 'in.data.' + str(itemp + 1)
			cmd5 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
		
			## ## ## Change output file
			orig = 'OUT_CORR'
			new  = 'out.GK.LJ.corr.' + str(iseed + 1) + '.' + str(itemp + 1) + '.' + str(istrain + 1)
			cmd6 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '

			## ## ## Change output file
			orig = 'OUT_STRAIN'
			new  = 'out.GK.LJ.strain.' + str(iseed + 1) + '.' + str(itemp + 1) + '.' + str(istrain + 1)
			cmd7 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '

			## ## ## Change volume file
			orig = 'OUT_VOLUME'
			new  = 'out.GK.LJ.vol.' + str(iseed + 1) + '.' + str(itemp + 1) + '.' + str(istrain + 1)
			cmd8 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
		
			## ## ## Change position file
			orig = 'OUT_POS'
			new  = 'out.GK.LJ.pos.' + str(iseed + 1) + '.' + str(itemp + 1) + '.' + str(istrain + 1)
			cmd9 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
		
			## ## ## Change main log file
			orig = 'MAIN_LOG_FILE'
			new  = 'out.GK.LJ.main.' + str(iseed + 1) + '.' + str(itemp + 1) + '.' + str(istrain + 1)
			cmd10 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
		
			## ## ## Change kappa log file
			orig = 'KAPPA_LOG_FILE'
			new  = 'out.GK.LJ.kappa.' + str(iseed + 1) + '.' + str(itemp + 1) + '.' + str(istrain + 1)
			cmd11 = '-e \'s/\<' + orig + '\>/' + new + '/g\' '
		
#	print 'sed ' + cmd1 + cmd2 + 'in.GK.LJ > ./in.GK.LJ.' + str(a + 1) 
	
			## ## ## Execute Command
			os.system('sed '+ cmd1+ cmd2+ cmd3+ cmd4+ cmd5+ cmd6+ cmd7+ cmd8+ cmd9+ cmd10+ cmd11+ 'in.GK.LJ > ./in.GK.LJ.'+ str(iseed + 1)+ '.'+ str(itemp + 1)+ '.'+ str(istrain + 1))
			#sed -e 's/\<POS_IN_STRUCTURE\>/in_structure_1.pos/g' -e 's/\<OUT_CORR\>/out.GK.LJ.corr.1/g' in.GK.LJ > ./in.GK.LJ.1
## ## ## END FILE
