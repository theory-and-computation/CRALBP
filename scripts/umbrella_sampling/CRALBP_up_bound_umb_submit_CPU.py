#!/usr/bin/python

"""submit.py: A simple script for submitting NAMD runs

  Usage: python submit.py <start> <end> """

import os, sys, re

#----------------------#
# Function Definitions #
#----------------------#

def gen_sub(jobstep):
 """ Generates an LSF submission script and returns the filename """
  
 script_name = "umb-%d.sh" % jobstep
 output = open(script_name, 'w') 
 output.write(r"""#!/bin/bash

#SBATCH -p free             # Queue (partition) name
#SBATCH -A alviny6_lab      # Account to charge
#SBATCH --nodes=8               # Total # of nodes 
#SBATCH --ntasks=48              # Total # of mpi tasks (4 * nnodes)
#SBATCH --cpus-per-task=8       # cores per task
#SBATCH -t 00:20:00             # Run time (hh:mm:ss)
#SBATCH --constraint="intel,avx512,mlx5_ib,fastscratch,nvme"
#SBATCH --exclusive

#SBATCH -J umb-p%d  # Job name
#SBATCH -o log/bound_upperpore-umb-%d.out          # Name of stdout output file

export window='%d'
export base_jobstep='00003'

module load gcc/11.2.0
module load openmpi/4.1.2/gcc.11.2.0
module load cuda/11.7.1

# Set environmental variables
export basename='umb-POPS_bound_smallChargedMembrane-ions-wats'
export nsteps='10000000'
export temp='309.2162'

export basedir='output'
export rstname=$basedir'/lp_'${basename}-${base_jobstep}'/lp_'${basename}-${base_jobstep}'.restart'
export rundir=$basedir'/up_'${basename}-${window}
export output=$rundir'/up_'${basename}-${window}

export firstjob='no'

export langPistPeriod='200'
export langPistDecay='100'
export rstFreq='10000'

# Set up directories
mkdir -p $rundir

# Store the output trajectory file names
echo $PWD'/'$output'.dcd' >> $basedir'/run.stk'

# Run NAMD
../software/ibverbs_NAMD_2.14_Source/Linux-x86_64-g++/charmrun      \
        +p336 ++mpiexec ++remote-shell          \
        mpirun ../software/ibverbs_NAMD_2.14_Source/Linux-x86_64-g++/namd2-ib       \
        ++ppn 7 \
        +pemap 2-14:2,3-15:2,18-30:2,19-31:2,34-46:2,35-47:2 \
        +commap 0,1,16,17,32,33 \
        umb_b_up-%d.conf > log/'up_'${basename}=${window}.log
""" % (jobstep, jobstep, jobstep, jobstep))
 output.close()
 os.system('chmod +x %s' % script_name) 
 return script_name 

def submit(fname, params=''):
 """ Submits the LSF job with fname, and returns the jobID """
  
 out = os.popen(f"sbatch {params} < {fname}").read()
 out = out.strip()
 print(out) 
 m = re.match('Submitted batch job (\d+)', out)
 return m.group(1) 

#--------------#
# Main Routine #
#--------------#

start = -8
end  = -8

for cycle in range(start, end + 1):

 script = gen_sub(cycle)
 jobid = submit(script)


