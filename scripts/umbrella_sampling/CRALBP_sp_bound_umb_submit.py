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

#SBATCH -p gpu             # Queue (partition) name
#SBATCH -A ALVINY6_LAB_GPU      # Account to charge
#SBATCH --nodes=1               # Total # of nodes 
#SBATCH --ntasks=1              # Total # of mpi tasks (4 * nnodes)
#SBATCH --cpus-per-task=40       # cores per task
#SBATCH --gres=gpu:V100:4       #Specify 4 V100 GPUs
#SBATCH -t 03:00:00             # Run time (hh:mm:ss)
#SBATCH --constraint="intel,avx512,mlx5_ib"

#SBATCH -J umb-p%d  # Job name
#SBATCH -o log/bound_sidepore-umb-%d.out          # Name of stdout output file

export window='%d'
export base_jobstep='00003'

module load gcc/11.2.0
module load openmpi/4.1.2/gcc.11.2.0
module load cuda/11.7.1

export OMP_NUM_THREADS=1
export OMPI_MCA_osc=^ucx
export UCX_TLS=rc,mm
export UCX_NET_DEVICES=mlx5_0:1

# Set environmental variables
export basename='umb-POPS_bound_smallChargedMembrane-ions-wats'
export nsteps='5000000'
export temp='309.2162'

export basedir='output'
export rstname=$basedir'/lp_'${basename}-${base_jobstep}'/lp_'${basename}-${base_jobstep}'.restart'
#export rstname=$basedir'/sp_'${basename}-${window}'/sp_'${basename}-${window}'.restart'

export rundir=$basedir'/sp_'${basename}-${window}
export output=$rundir'/sp_'${basename}-${window}

export firstjob='no'

export langPistPeriod='200'
export langPistDecay='100'
export rstFreq='10000'

# Set up directories
mkdir -p $rundir

# Store the output trajectory file names
echo $PWD'/'$output'.dcd' >> $basedir'/run.stk'

# Run NAMD
namd3 +p29 +pmepes 5 +setcpuaffinity +devices 0,1,2,3 umb_b_sp-%d.conf > log/'sp_'${basename}-${window}.log   # with 8 PEs per non-PME


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

start = -10
end  = 10

for cycle in range(start, end + 1):

 script = gen_sub(cycle)
 jobid = submit(script)

