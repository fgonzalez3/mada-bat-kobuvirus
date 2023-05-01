# Open Reading Frame (ORF) Kobuvirus Phylogeny (aa) 

This tutorial outlines methods for building a maximum-liklihood phylogenetic tree using RAxML. Here, we'll build an ORF phylongey of Kobviruses on [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) while fitting our new Madagascar whole genome sequence in. 

# Kobuvirus ORF Amino Acid Tree (Last Access: 6/23/2022)

1. To build this tree, I compiled the same list used for Fig3B. This list contained full- or near-full genome length Kobuvirus sequences, from both NCBI Virus and IDseq. A good handfull of these sequences, references included, were not previously annotated. I went ahead and used the 'predict ORF' function in GeneiousPrime to rule out sequences that did not contain full ORFs (and thus be removed for this analysis). 

Sequences containing ORFs were further annotated in GeneiousPrime using reference-based annotations. This isn't necessary to make this tree, but I did this in case there were any specific polyproteins that were common enough among sequences to make a supp. tree figure (which I didn't end up doing anyways). 

Once this list of sequences was compiled, ORFs were extracted and translated. Following translation, sequences were aligned using MAFFT under default parameters within GeneiousPrime. 

De-duplication was not necessary, since these were the same sequences used in our prior analysis for Fig2B. However, [ModelTest-NG](https://github.com/ddarriba/modeltest) was still ran in order to compare amino acid substitution models needed to construct our tree. The alignment used for ModelTest can be found under [Fig2C-ModelTest](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig2C/modeltest). The following shellscript was ran on the HCC: 

```
 #!/bin/bash
#SBATCH --job-name=ModelTest-NG
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
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


modeltest-ng -d aa -i ORF_alignment.fasta -t ml -p 8
```

The final list of sequences used in this analysis can be found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/blob/main/Fig2C/kobu_orf.csv). ModelTest-NG recovered a LG+G4+F amino acid subtitution model, which informed our RAxML tree construction. ModelTest results for Fig2C can be found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig2C/modeltest). RAxML was ran using the following script on the HCC: 

```

#!/bin/bash
#SBATCH --job-name=orfraxml
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1024  
#SBATCH --time=168:00:00
#SBATCH --output=orfraxml.out
#SBATCH --error=orfraxml.err

module purge
module load vim/8.1
module load emacs/27.2
module load python/3.9
module load java/12
module load cmake/3.20
module load compiler/gcc/9 
module load openmpi/4.1 
module load raxml-ng/1.1

raxml-ng-mpi --all --msa orf_raxml.fasta --model LG+G4+F --prefix T3  --seed 12 --threads 6 --bs-metric fbp,tbe

```

Once RAxML finished, the resulting tree was imported into R and visualized using [ggtree](https://yulab-smu.top/treedata-book/chapter9.html). Note that for quick viewing, you can import the RAxML output file into the open source program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). Rscript to create Fig2C can be found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig2C). 
