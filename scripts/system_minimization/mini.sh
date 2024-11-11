#!/bin/sh

#SBATCH -p gpu          # Queue (partition) name
#SBATCH -A ALVINY6_LAB_GPU      # Account to charge
#SBATCH --nodes=1            # Total # of nodes 
#SBATCH --ntasks=1           # Total # of mpi tasks
#SBATCH --cpus-per-task=40   # Cores per task
#SBATCH --gres=gpu:V100:4    # Specify 4 V100 GPUs
#SBATCH -t 00:30:00          # Run time (hh:mm:ss)
#SBATCH --constraint="intel,avx512,mlx5_ib"

#SBATCH -J rsmc-mini	   # Job name
#SBATCH -o log/rotatedSmallChargedMembrane_2-ions-wats.out       # Name of stdout output file

export jobstep='00000'

module load gcc/11.2.0
module load openmpi/4.1.2/gcc.11.2.0
module load cuda/11.7.1

export OMP_NUM_THREADS=1
export OMPI_MCA_osc=^ucx
export UCX_TLS=rc,mm
export UCX_NET_DEVICES=mlx5_0:1

# Set environmental variables
export basename='rotatedSmallChargedMembrane_2-ions-wats'
export nsteps='100000'
export temp='310'

export basedir='output'
export rundir=$basedir'/'$basename
export output=$rundir'/'$basename

export firstjob='yes'

# Set up directories
mkdir -p $rundir

# Store the output trajectory file names
echo $PWD'/'$output'.dcd' >> $basedir'/run.stk'

# Run NAMD
#ibrun namd2 mini.conf +ppn 4 +pemap 0-51 +commap 52-67 > $basename'-'$jobstep'.log'
namd3 +p29 +pmepes 5 +setcpuaffinity +devices 0,1,2,3 mini.conf > log/$basename'-'$jobstep'.log'
