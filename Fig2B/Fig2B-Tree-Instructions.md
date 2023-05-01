# Full Genome Kobuvirus Phylogeny (nt)

This tutorial outlines methods for building a maximum-liklihood phylogentic tree using RAxML. Here, we'll build a full-genome length phylogeny of Kobuviruses on [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) while fitting our new Madagascar whole genome sequence in. 

# Full Genome Kobuvirus NCBI Virus Downloads (Last Accesss: 6/23/2022)

1. I first went into NCBI-virus and selected for the following:

- Minimum sequence length 5212bp and a maximum sequence length of 9000bp
- All sequences available for Kobuvirus (taxid: 194960; 109 hits) with complete nucleotide completness
- Additional Kobuvirus TaxID's did not yield more hits 

Because there are significantly less Kobuvirus sequences on NCBI, we did not limit our search to reference sequences as we did in Fig2A. 5.2kb was our smallest postulated Kobuvirus, so this was selected as our min. When this sequence was later found to be bacterial, methods were kept the same. Only the sequence was removed. 

2. Rat Rabovirus (Accession: NC_026314), was chosen as an outgroup. This sequence, though appearing as Kobuvirus using the above search parameters, was not a Kobuvirus according to further analysis. Additionally, bat Picornaviridae that clustered closely to known Kobuviruses in Fig2A were added to Fig2B. This was done to deduce whether these bat Picornaviruses may have been classified incorrectly. Those that met this criteria are listed below: 

- Miniopterus Picornavirus (NCBI Accession: MF352427)
- Bat Picornavirus (NCBI Accession: MN602325)

3. All .fasta files were imported into Geneious Prime and aligned using MAFFT under default parameters. Our alignment file was further trimmed to encompass a consensus region within all sequences, excluding gaps at either end which might interfere with tree analysis. This ultimately led to a full genome alignment, excluding the 3' and 5' UTRs. Our alignment file can be found under [Fig2B-RAxML](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig2B/raxml) as 'kobu_raxml.fasta'. I used R to edit the names of each sequence in MSA. RAxML won't accept spaces, periods, dashes, slashes, colon, semicolons, or parentheses in sequence names. Rscripts for this can be found in [TreePrep subfolder](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/TreePrep). 

4. I submitted the aligned .fasta file to the HCC to kick off [ModelTest-NG](https://github.com/ddarriba/modeltest), which simultaneously flags duplicated sequences and compares nucleotide substituion models. The following shell script was used, though yours might differ: 

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


modeltest-ng -d nt -i kobuvirus_alignment.fasta -t ml -p 8

```

After some pruning and figuring out that we really only had one full genome Kobuvirus, the sequences found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/blob/main/Fig2B/kobuvirus_manual2.csv) made the final cut. ModelTest-NG results for Fig2B can be found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig2B/modeltest). 

5. Once ModelTest-NG finished, a maximum likelihood tree using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) was built. ModelTest-NG support recovered a GTR+I+G4 nucleotide subtitution model. This was ran using the following script on the HCC: 

```

#!/bin/bash
#SBATCH --job-name=koburaxml
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1024  
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

raxml-ng-mpi --all --msa kobu_raxml.fasta --model GTR+I+G4 --prefix T3  --seed 12 --threads 6 --bs-metric fbp,tbe 

```

Once RAxML finished, the resulting tree was imported into R and visualized using [ggtree](https://yulab-smu.top/treedata-book/chapter9.html). Note that for quick viewing, you can import the RAxML output file into the open source program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). Rscript to create Fig2B can be found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig2B). 



