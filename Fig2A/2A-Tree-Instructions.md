# Full Genome Picornavirus Phylogeny (nt)

This tutorial outlines methods for building a maximum-likliehood phylogenetic tree using RAxML. Here, we'll build a full-genome length phylogeny of Picornaviruses on [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) while fitting our new Madagascar whole genome sequence in. 

# Full Genome Picornavirus NCBI Virus Downloads (Last Access 6/23/2022)

1. I first went into NCBI Virus and selected for the following: 

- Minimum sequence length of 7202bp and max of 9000bp 
- All reference sequences with complete nucleotide completeness for Picornaviridae (taxid 12058; 130 hits)
- TaxID's (unclassified Picornaviridae, tax id 478825; Picornaviridae sp. rodent/Ee/PicoV/NX2015, tax id 1917412) did not yield hits. 

At the time, our smallest postulated Kobuvirus sequence was 7202bp (later found to be a bacterial sequence). Choosing this minimum would eliminate any sequences smaller than ours, and make alignment steps easier later on. 

According to [ICTV](https://ictv.global/report/chapter/picornaviridae/picornaviridae), Picornaviruses cap out around 10kb while Kobuviruses cap out around 8.5kb. 9kb was a good cutoff. 

2. I then went back and did a non-refernece sequence search, selecting for those greater than 7.2kb. 17 unique hits were found in this manner, accounted for in the below TaxID's: 

- Bat Picornavirus, taxid 1281456
- Bat Picornavirus 1, taxid 1074863
- Bat Picornavirus 2, taxid 1074864
- Bat Picornavirus 3, taxid 1074865
- Washington bat picornavirus, taxid 1888309
- Aichivirus F, taxid 1986959
- Aichivirus F1, taxid 2758106
- Aichivirus F2, 2755032
- Miniopterus picornavirus 1, taxid 2184394
- Miniopterus picornavirus 2, taxid 2184395

I downloaded the above hits, in addition to our full genome hit (OP287812), as a .fasta file which would undergo de-duplication. 

3. Finally, I included a full-genome bat Coronavirus (Accession: NC_048212) as an outgroup. Coronaviruses and Picornaviruses are taxonomically in the same Class, distant enough to act as outgroups to one another. 

4. All .fasta files were  imported into Geneious Prime and aligned using MAFFT under default parameters. Our alignment file was further trimmed to encompass a consensus region within all sequences, excluding gaps at either end which might interfere with tree analysis. This ultimately led to a full genome alignment, excluding the 3' and 5' UTRs. This alignment file can be found under [Fig2A-RAxML](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig2A/raxml) as 'picorna_raxml.fasta'. I used R to edit the names of each sequence in MSA. RAxML won't accept spaces, periods, dashes, slashes, colon, semicolons, or parentheses in sequence names. Rscripts for this can be found [TreePrep subfolder](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/TreePrep). 

5. I submitted the aligned .fasta file to HCC to kick off [ModelTest-NG](https://github.com/ddarriba/modeltest), simultaneously flagging duplicated sequences and comparing nucleotide substition models. The following shellscript was used, though yours might differ: 

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


modeltest-ng -d nt -i picornavirus_alignment.fasta -t ml -p 8

```

After some pruning and figuring out that we really only had one full genome Kobuvirus, the sequences found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/blob/main/Fig2A/picorna_manual.csv) made the final cut. ModelTest-NG results for Fig2A can be found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig2A/modeltest). 

6. Once ModelTest-NG finished, a maximum liklihood tree using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) was built. ModelTest-NG support recovered a TVM+G4 nucleotide substitution model. This was ran using the following script on the HCC: 

```

#!/bin/bash
#SBATCH --job-name=picornaraxml
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1024  
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

raxml-ng-mpi --all --msa picorna_raxml.fasta --model TVM+I+G4 --prefix T3  --seed 12 --threads 8 --bs-metric fbp,tbe 

```

Once finished, the resulting tree was imported into R and visualized using [ggtree](https://yulab-smu.top/treedata-book/chapter9.html). Note that for quick viewing, you can visualize the RAxML output file in the open source program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). Rscript to create Fig2A can be found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig2A). 

