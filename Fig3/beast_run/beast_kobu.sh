#!/bin/bash
#SBATCH --job-name=beastkobu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --mem=30gb
#SBATCH --output=beastkobu.out
#SBATCH --error=beastkobu.err

module load vim/8.1
module load java/1.8
module load emacs/27.2
module load cmake/3.17
module load python/3.9
module load compiler/gcc/10
module load openmpi/4.1
module load beagle/5.2
module load beast2/2.6

export CC=`which gcc`
export CXX=`which c++`

export LD_LIBRARY_PATH=/work/cresslerlab/fgonzalez3/beagle-lib/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/work/cresslerlab/fgonzalez3/beagle-lib/lib/pkgconfig:$PKG_CONFIG_PATH
export BEAST_EXTRA_LIBS=/work/cresslerlab/fgonzalez3/beagle-lib/lib/:$BEAST_EXTRA_LIBS

beagle -Xms512m -Xmx10g nthreads=8
beast -threads 8 -beagle_GPU -seed 777 beast_kobu.xml
