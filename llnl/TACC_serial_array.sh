#!/bin/bash

#SBATCH -J TO_D${SAMPLE_SIZE}_${PARAM1}		# job name
#SBATCH -o TO_D${SAMPLE_SIZE}_${PARAM1}.o%j	# output and error file name (%j expands to jobID)
#SBATCH -n 1               			# total number of mpi tasks requested
#SBATCH -N 1   					# request 1 nodes
#SBATCH -p normal          			# queue (partition) -- normal, development, etc.
#SBATCH --mail-user=ill.ryu@utdallas.edu
#SBATCH --mail-type=begin  			# email me when the job starts
#SBATCH --mail-type=end   		 	# email me when the job finishes
#SBATCH  -A Fundamental-investig

## Run a command
# ibrun ./bin/FCC_TENSION_kMC_delocal_Coeff tests/FCC_TENSION_kMC_delocal_Coeff2/D${SAMPLE_SIZE}_${PARAM1}.ctrl > tests/FCC_TENSION_kMC_delocal_Coeff2/D${SAMPLE_SIZE}_${PARAM1}.out 
ibrun ./bin/FCC_TORSION_kMC_delocal_Coeff tests/FCC_TORSION_after_TENSION/D${SAMPLE_SIZE}_${PARAM1}.ctrl > tests/FCC_TORSION_after_TENSION/D${SAMPLE_SIZE}_${PARAM1}.out  

