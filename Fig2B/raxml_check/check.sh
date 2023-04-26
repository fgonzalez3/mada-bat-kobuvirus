#!/bin/bash
#SBATCH --job-name=koburaxml
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=168:00:00
#SBATCH --output=koburaxml.out
#SBATCH --error=koburaxml.err

module purge
module load vim/8.1
module load emacs/27.2
module load python/3.9
module load java/12
module load cmake/3.20
module load compiler/gcc/9 
module load openmpi/4.1 
module load raxml-ng/1.1

raxml-ng-mpi --check --msa kobu_raxml.fasta --model GTR+I+G4 --prefix T1
