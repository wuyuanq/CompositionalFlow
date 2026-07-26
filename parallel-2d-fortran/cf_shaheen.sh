#!/bin/bash

#SBATCH -N 1                           #Number of nodes
#SBATCH --ntasks=16             #Number of cores
#SBATCH --exclusive                      #We need to run on nodes where nobody will interfere with us
#SBATCH --time=24:00:00
#SBATCH -J cf          #Job name
#SBATCH -o cf.%j               #File to which standard out will be written
#SBATCH -e cf.%j.err                #File to which standard err will be written
date
time srun -n 16 --hint=nomultithread --ntasks-per-node=32 --ntasks-per-socket=16 --ntasks-per-core=1 ./CompositionalFlow_fortran_hpc
date