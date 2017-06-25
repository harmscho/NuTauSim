#!/bin/bash
#PBS -q verylongq
#PBS -l select=1:ncpus=1
#PBS -l walltime=192:00:00
#PBS -J 0-1754

# NOTE
# '#PBS' directives must immediately follow your shell initialization line '#!/bin/<shell>'
# '#PBS' directives must be consecutively listed without any empty lines in-between directives
# Reference the PBS Pro User Guide for #PBS directive options.
# To determine the appropriate queue (#PBS -q) and walltime (#PBS -l walltime) to use,
#  run (qmgr -c 'print server') on the login node.

# This is an example "Job Array" job script
# PBS Job Arrays are typically used for embarrassingly parallel codes
# The "#PBS -J" directive specifies an index range, which correlates to the amount of jobs spawned
# Input files can be named according to the index (input.1, input.2, ..etc)
# Reference the PBS Pro User Guide to learn about the power of Job Arrays.

# Set the output directory
# By default, PBS copies stdout and stderr files back to $PBS_O_WORKDIR
# When 'qsub' is run, PBS sets $PBS_O_WORKDIR to the directory where qsub is run.
# Change this environment variable if desired
#
export PBS_O_WORKDIR=/aurora_nobackup/eva/romerowo/nutau_outputs/20170327

mkdir $PBS_O_WORKDIR

# Set your input directory (optional)
#
export INPUT_DIR=/home/romerowo/nutau_sim/inputs/eventSim

# Change directories to the local compute node's scratch directory 
# Using each node's local scratch storage can greatly improve IO
# Any file created in $TMPDIR will be automatically cleaned up by PBS
# Aurora's and Halo's $TMPDIR by default is set to '/lscratch'
# Zodiac does not have local compute node storage
# 
export TMPDIR=/lscratch

cd $TMPDIR

# Load software modules
# Available modules can be found with the command 'module avail'
#
#module load <module_name>

# Run executable
# Job Arrays step through the index range specified by '#PBS -J'
# The value of $PBS_ARRAY_INDEX changes as each consecutive job is launched
# Using $TMPDIR will write to the local compute node's scratch storage, but remember to
#  move your files back to your output directory
#
/home/romerowo/nutau_sim/run_python_command.py $INPUT_DIR/input.$PBS_ARRAY_INDEX > output.$PBS_ARRAY_INDEX
#
## output currently in local node scratch
#
mv output.$PBS_ARRAY_INDEX $PBS_O_WORKDIR
#
## output moved to desired output directory
#
