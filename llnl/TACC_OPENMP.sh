#!/bin/bash
#SBATCH -J myMPI           # job name
#SBATCH -o myMPI.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 16              # total number of threads task requested
#SBATCH -N 1               # total number of nodes task requested
#SBATCH -p normal          # queue (partition) -- normal, development, etc.
#SBATCH -t 48:00:00        # run time (hh:mm:ss) - 2days
#SBATCH --mail-user=ill.ryu@utdallas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

## Run a command
ibrun runcyl tests/test_MOOSE/D150_1e5.ctrl 
