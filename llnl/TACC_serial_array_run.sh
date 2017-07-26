#!/bin/bash

for SAMPLE_SIZE in 1000; do
   for PARAM1 in 1 2 3 4 5 ;
      do
#
      echo "D=${SAMPLE_SIZE}nm, case=${PARAM1}"
      export PARAM1 SAMPLE_SIZE
#
      sbatch -t 48:00:00 \
	TACC_serial_array.sh
#
      sleep 1 # pause to be kind to the scheduler
   done
done
