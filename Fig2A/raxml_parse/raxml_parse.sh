#!/bin/bash
#SBATCH --job-name=picornaraxml
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=168:00:00
#SBATCH --output=picornaraxml.out
#SBATCH --error=picornaraxml.err

module purge
module load vim/8.1
module load emacs/27.2
module load python/3.9
module load java/12
module load cmake/3.20
module load compiler/gcc/9 
module load openmpi/4.1 
module load raxml-ng/1.1

raxml-ng-mpi --parse --msa picorna_raxml.fasta --model TVM+I+G4 --prefix T2
