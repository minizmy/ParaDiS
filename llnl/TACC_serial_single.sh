#!/bin/bash
#SBATCH -J TOR_D300_1			  	# job name
#SBATCH -o TOR_D300_1.o%j 			# output and error file name (%j expands to jobID)
#SBATCH -n 1               			# total number of mpi tasks requested
#SBATCH -N 1 	 	 			# request 1 nodes
#SBATCH -p normal          			# queue (partition) -- normal, development, etc.
#SBATCH -t 02:00:00     		   	# run time (hh:mm:ss) - 1mins
#SBATCH --mail-user=ill.ryu@utdallas.edu
#SBATCH --mail-type=begin  			# email me when the job starts
#SBATCH --mail-type=end   		 	# email me when the job finishes
#SBATCH  -A Fundamental-investig

## Run a command
#ibrun ./bin/FCC_TORSION_kMC_delocal_Coeff tests/FCC_TORSION_kMC_delocal_Coeff/D1000_1.ctrl > tests/FCC_TORSION_kMC_delocal_Coeff/D1000_1.out 
ibrun ./bin/FCC_TORSION_kMC_delocal_Coeff tests/FCC_TORSION_kMC_delocal_Coeff/D300_1.ctrl > tests/FCC_TORSION_kMC_delocal_Coeff/D300_1.out 
