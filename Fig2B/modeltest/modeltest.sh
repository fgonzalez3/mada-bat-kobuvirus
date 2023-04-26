#!/bin/bash
#SBATCH --job-name=ModelTest-NG
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --mem=30gb
#SBATCH --output=ModelTest.out
#SBATCH --error=ModelTest.err

module load compiler/gcc/9 openmpi/4.1 modeltest-ng/0.1
module load vim/8.1 
module load java/12
module load emacs/27.2
module load cmake/3.5
module load python/3.9

export CC=`which gcc`
export CXX=`which c++`


modeltest-ng -d nt -i kobuvirus_alignment.fasta -t ml -p 8
