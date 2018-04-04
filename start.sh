#!/bin/bash -l
#
# allocate 16 nodes (64 CPUs) for 6 hours
#PBS -l nodes=4:ppn=40,walltime=01:00:00
#
# job name 
#PBS -N TestAval
#
# stdout and stderr files
#PBS -o aval.out -e aval.err
#
#PBS -q devel
#
# first non-empty non-comment line ends PBS options

# jobs always start in $HOME -
# change to a temporary job directory on $FASTTMP

# copy input file from location where job was submitted
EXEC="$HOME/Projects/hmn-percolation/percolation"

INPUTFILE='FileReaderTestInput.txt'
SIMPATH="$HOME/Projects/hmn-percolation/Results/alpha10/20/p045"

# run
HOSTN=$(hostname)
if [ ! "${HOSTN#ww8}" == "$(hostname)" ]; then

  module load intel
  cd $SIMPATH
  echo  "Hi:   $SIMPATH/$INPUTFILE"
  $EXEC "$SIMPATH/$INPUTFILE"
else
  echo  "Hi:   $SIMPATH/$INPUTFILE"
  mkdir ${HOME}/$PBS_JOBID
  cd ${HOME}/$PBS_JOBID
  module load intel64
  module load intelmpi
  module load mkl
  mpirun_rrze -npernode 20 $EXEC -in $SIMPATH/$INPUTFILE
  cp processed_data.dat $SIMPATH

  cd
  # get rid of the temporary job dir
  rm -rf ${HOME}/$PBS_JOBID
fi
# save output on parallel file system
#mkdir -p ${HOME}/codes/SIMSRC/$PBS_JOBID


