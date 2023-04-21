### Plotting seroprevalence using using full genome Kobuvirus as reference

1.First, I downloaded all non-host contigs derived from fecal samples, mapping to ANY taxon from IDseq in the 'RR034B1_feces' project. (Note that we did not take contigs from HeLa controls or water). This is done in bulk, manually, on IDseq.net in the top right-hand corner. 

Note that when you download all the non-host contigs, it will produce a folder with a separate fasta file for each sample, which lists the contigs by node number but does not include the sample ID. Before joining all the contigs (nodes) together, you need to distinguish them by sample ID. Cara wrote an Rscript that parses this for each filetype (rename-fastas-feces, rename-fastas-urine, rename-fastas-throat). To rename your files and the headers within them, copy the appropriate Rscript for the tissue type into your appropriate downloads folder, cd into that folder on the command line, and simply type (example here for feces):

*Note: The R script kept gettig hung up on contig file 193, so it was omitted. 

```
Rscript rename-fastas-feces.R
```

Once the file finishes running, concatenate all the abbreviated file names into a compiled file for downstream analyses. Note that these files are too big to include in the Github repo. Within working directory where all files are abbreviated files are located, run:

```
cat *.fasta > meta_fec.fasta
```

### Running CD-HIT

2.Next, upload to home cluster. Upload using OnDemand if it's offered run a script using the commmand _scp_: 

```
scp ./meta_fec.fasta fgonzalez3@crane.unl.edu:/work/cresslerlab/fgonzalez3
```

3. Now, because many of the contigs will overlap in sequences and slow our mapping down, we can deduplicate the compiled .fasta for each of the contigs within the fecal samples. 

You can run CD-HIT on your local computer or on the computing cluster of your choice. The latter choice is quickest, especially if your home cluster already has CD-HIT downloaded. Otherwise, you can use bioconda to load programs like CD-HIT onto an environment that you can utilize. Either of these reduce run time. Yale Ruddle luckily had CD-HIT downloaded, so I used the below shellscript but you can find instructions on how to create a bioconda environment in step 5.  

```
#!/bin/bash
#SBATCH --job-name=CD-HIT-fec
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --mem=30gb
#SBATCH --output=CD-HIT-fec.out
#SBATCH --error=CD-HIT-fec.err

module load CD-HIT/4.8.1-GCCcore-10.2.0

cd-hit-est -i meta_fec.fasta -c 0.95 -o meta_fec_DEDUP.fasta -M 0
```

4. You'll notice that a lot of the beginning steps included in this pipeline are the same as those found in our initial BLAST search. The biggest change occurs here, where I use the full Kobuvirus genome we found during our initial search as the reference for this search. 

### Running Blastn/Blastx using full genome as reference

5. After contigs are deduplicated, I downloaded my full genome Kobuvirus (genbank: OP287812) that I'm going to use as a reference sequence. I downloaded the nt and aa .fasta files. Below are the names I gave these files:

```
ref_nt.fasta 
ref_aa.fasta
```

6. I previously made a conda environment on Yale Ruddle, containing blast, for another project so I worked within this environment instead of installing NCBI-Blast directly on my computer. See Cara's workflow if you wish to install on your local computer ([here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/contig-blast-directions.md)). To create this environment, I followed the proceeding steps: 

5a. Ask for space off the login node 

```
salloc
```

If you have larger data to work with, ask for additional memory

```
salloc --mem=8G
```

6b. Load the Miniconda module (use 'ml' or 'module load')

```
ml miniconda
```

6c. Use Mamba to create an environment (this is where you will include all packages that you'll need for any future analyses as well). This is for a primer design project, hence the primer_design environment name and the packages contained within. Use another name and download any packages you might need for BLAST or other future analyses here. 

```
mamba create -yn primer_design primalscheme snippy
```

6d. Once built, end session and free resouces 

```
exit
```

6e. Activate your environment

```
conda activate primer_design
```

Once activated, your environment should look something like this if you're on the cluster terminal:

```
(primer_design)[flg9@c15n08 blastn]$
```

6f. From there you should be able to use 'conda install PACKAGENAME' inside the activated environment to install other packages in the environment if need be. If an error occurs at this step, delete the original environment and create a new one reflecting the packages you wish to use. To delete environment:

```
conda env remove -yn ENVNAME
```

7. After having activated our conda working environment and navigating to the working directory on Ruddle that contained my input refseq, I made a reference database from my downloaded sequences. To do this, I used the following commands on each of my full genome .fasta files:

```
makeblastdb –in ref_nt.fasta –dbtype nucl –parse_seqids -out KoV_nt
makeblastdb –in ref_aa.fasta –dbtype prot –parse_seqids -out KoV_aa
```

The commands above should generate a suite of files in the same folder that you are working within. 

8. Now that our reference sequence is in hand and the non-host contigs are deduplicated, kick off a command line blast for each contig subset on the two above databases (2 runs total). I ran these on Yale Ruddle. Basically, we are re-doing a more focused version of IDseq to see if any other "hits" to KoVs specifically shake out. 

You will run two BLASTs total: 1 nucleotide and 1 protein BLAST on each of the deduplicated contig samples above. First, run a BLASTn alignment of deduplicated set of contigs with the KoV-nt database, then run a BLASTx aligment of the deduplicated set of contigs with the KoV_aa database. Scripts for both are listed below:

```
blastn -word_size 10 -evalue 0.001 -query meta_throat_DEDUP.fasta -db KoV_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out Mada_Bat_KoV_blast_nt.txt

blastx -word_size 3 -evalue 0.001 -query meta_throat_DEDUP.fasta -db KoV_aa -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out Mada_Bat_KoV_blast_aa.txt
```
Again, these can be run locally but I used Ruddle to speed things up. 

9. From here, I followed the same protocol outline in my initial BLAST search pipeline ___. 

