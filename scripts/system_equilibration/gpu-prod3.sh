#!/bin/bash

#SBATCH -p gpu             # Queue (partition) name
#SBATCH -A ALVINY6_LAB_GPU      # Account to charge
#SBATCH --nodes=1               # Total # of nodes 
#SBATCH --ntasks=1              # Total # of mpi tasks (6 * nnodes)
#SBATCH --cpus-per-task=40      # cores per task
#SBATCH --gres=gpu:V100:4       # Specify 4 V100 GPUs
#SBATCH -t 00:10:00             # Run time (hh:mm:ss)
#SBATCH --constraint="intel,avx512,mlx5_ib"

#SBATCH -J rcm2_3  # Job name
#SBATCH -o log/rotatedSmallChargedMembrane_2-ions-wats-3.out          # Name of stdout output file

export jobstep='00003'
export prev_jobstep='00002'


module load gcc/11.2.0
module load openmpi/4.1.2/gcc.11.2.0
module load cuda/11.7.1

export OMP_NUM_THREADS=1
export OMPI_MCA_osc=^ucx
export UCX_TLS=rc,mm
export UCX_NET_DEVICES=mlx5_0:1

# Set environmental variables
export basename='rotatedSmallChargedMembrane_2-ions-wats'
export nsteps='50000'
export temp='310'

export basedir='output'
export rstname=$basedir'/'${basename}-${prev_jobstep}'/'${basename}-${prev_jobstep}'.restart'
export rundir=$basedir'/'${basename}-${jobstep}
export output=$rundir'/'${basename}-${jobstep}

export firstjob='no'

export langPistPeriod='200'
export langPistDecay='100'
export rstFreq='1000'

# Set up directories
mkdir -p $rundir

# Store the output trajectory file names
echo $PWD'/'$output'.dcd' >> $basedir'/run.stk'

# Run NAMD
#ibrun namd2 equil_npt.conf +ppn 4 +pemap 0-51 +commap 52-67 > log/$basename'-'$jobstep'.log'
#namd2-src +p 10 equil_npt.conf > log/$basename'-'$jobstep'.log'
namd3 +p29 +pmepes 5 +setcpuaffinity +devices 0,1,2,3 equil_rcmnpt.conf > log'/'${basename}-${jobstep}.log   # with 8 PEs per non-PME device, +p is set to 29 = 3*8 + 5
