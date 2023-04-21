### Finding positive Kobuvirus samples after mNGS + IDseq

1.First, we downloaded all non-host contigs derived from fecal samples, mapping to ANY taxon from IDseq in the 'RR034B1_feces' project. (Note that we did not take contigs from HeLa controls or water). This is done in bulk, manually, on IDseq.net in the top right-hand corner. 

Note that when you download all the non-host contigs, it will produce a folder with a separate fasta file for each sample, which lists the contigs by node number but does not include the sample ID. Before joining all the contigs (nodes) together, you need to distinguish them by sample ID. Cara wrote an Rscript that parses this for each filetype (rename-fastas-feces, rename-fastas-urine, rename-fastas-throat). To rename your files and the headers within them, copy the appropriate Rscript for the tissue type into your appropriate downloads folder, cd into that folder on the command line, and simply type (example here for feces):

```
Rscript rename-fastas-feces.R
```

2. Next, upload to home cluster. You can use the example script below: 

```
scp ./bat_feces_renamed.fasta fgonzalez3@crane.unl.edu:/work/cresslerlab/fgonzalez3
```
OR 

Upload the files to your cluster using OnDemand if it's offered. 

Once the file finishes running, concatenate all the abbreviated filenames into a compiled file for downstream analyses. Note that these files are too big to include in the GitHub repo.

3. Now, because many of the contigs will overlap in sequences and slow our mapping down, we can deduplicate the compiled .fasta for each of the contigs for all tissues using the program CD-HIT. You'll need to install this on your home computer first. If using MacOS, Cara reccomends using bioconda to do it. Check __ where you can find steps needed to create a conda environment. It can be run locally, but is extremely slow. 

The cluster I was using at the time (Nebraska HCC) had CD-HIT installed, so I went about this a little differently. I deduplicated my contig files via the following command line script in the same folder as the downloaded contigs (here shown just for the fecal contigs file -- you will need to run it for all three tissue types). 

Here's the bash script I used on the HCC to run the fecal example:

```
#!/bin/bash
#SBATCH --job-name=CD-HIT-feces
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=168:00:00
#SBATCH --mem=30gb
#SBATCH --output=CD-HIT.%J.out
#SBATCH --error=CD-HIT.%J.err

module load cd-hit/4.8

cd-hit-est -i bat_feces_renamed.fasta -o bat_feces_renamed_DEDUP.fasta -M 0 -T 0
```

4. Next, after the contigs are deduplicated (or simulataneously as you are doing this), you can download all the (a) nucleotide and (b) protein full genome reference sequences under Kobuvirus from NCBI virus. There were roughly ~10 taxid's when I did this (Fall 2021). Next, I concatenated the Kobuvirus nucleotide and protein files such that I had only two Kobuvirus reference fasta files: one for nucleotides and one for amino acids. 


5. I ran all my BLAST searches on the HCC, initially in shell, since BLAST was already installed. I later migrated to another cluster (Ruddle) where BLAST was not installed. There, I ran BLAST again using Bioconda. Either approach works, but do check __ if BLAST is not already installed for you. Provided below is my initial approach. 


6. The first step toward running an offline BLAST search is making a reference database from your downloaded NCBI Kobuvirus sequences.

6a. Using shell --  

For nucleotides: 

```
#!/bin/bash
#SBATCH --job-name=kobu_nt_DB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --mem=30gb
#SBATCH --output=kobu_blast_nt_db.%J.out
#SBATCH --error=kobu_blast_nt_db.%J.err

module load blast/2.10

makeblastdb -in kobu_nt_seq.fasta -dbtype nucl -parse_seqids -out Kobu_nt_db
```

6b. 

For proteins: 

```
#!/bin/bash
#SBATCH --job-name=kobu_pt_db
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --mem=30gb
#SBATCH --output=kobu_blast_pt_db.%J.out
#SBATCH --error=kobu_blast_pt_db.%J.err

module load blast/2.10

makeblastdb -in kobu_protein_seq.fasta -dbtype prot -parse_seqids -out Kobu_prot_db
```

The commands above should generate a suite of files in the same folder that have the prefix specified after the "out" command. NOTE: generate two additional directories, containing your (1) non-host (a) nucleotide and (b) amino acid contigs and (2) the files generated above. Navigate to these when running the below BLAST searches. So it should look like this: 

BLASTx Directory --

- aa_non_host_contigs
- .Kobu_prot_db files
- BLASTx shellscript below

BLASTn Directory -- 

- nt_non_host_contigs 
- .Kobu_nt_db files
- BLASTn shellscript below

7. Now that the reference sequence is in hand and the non-host contigs deduplicated, you can kick off a command line blast for each contig subset on the two above databases (so two runs in total). Again, I initially ran these on HCC. Basically, you are re-doing a more focused version of IDseq to see if any other "hits" to Kobuvirus, that may have not been detected by IDseq, shake out. 

You will run two BLASTs in total: 1 nucleotide and 1 protein BLAST, one for each of the deduplicated contigs above. First, run a "blastn" analysis of deduplicated set of contigs with the Kobuvirus nucleotide database generated above, then run a "blastx" alignment of the deduplicated set of contigs with the Kobuvirus amino acid database generated above. Scripts for both are listed below:

Here is the bash script for the deduplicated fecal contigs for the nucleotide blast:

```
#!/bin/bash
#SBATCH --job-name=BlastN_kobu_feces
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=168:00:00
#SBATCH --mem=30gb
#SBATCH --output=BlastN_kobu_feces.%J.out
#SBATCH --error=BlastN_kobu_feces.%J.err

module load cmake/3.20
module load compiler/gcc/10	
module load blast/2.10
module load biodata/1.0

cp -r Kobu_nt/ /scratch/
cp bat_feces_renamed_DEDUP.fasta /scratch/

blastn -word_size 10 -query /scratch/bat_feces_renamed_DEDUP.fasta -db /scratch/Kobu_nt/Kobu_nt_db -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out /scratch/Mada_Bat_Kobu_blastn_feces_nt.txt -num_threads $SLURM_NTASKS_PER_NODE -evalue 1e-10

cp /scratch/Mada_Bat_Kobu_blastn_feces_nt.txt .
```

Here is the bash script for the deduped fecal contigs for the protein blast:

```
#!/bin/bash
#SBATCH --job-name=BlastX_kobu_feces
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=168:00:00
#SBATCH --mem=30gb
#SBATCH --output=BlastX_kobu_feces.%J.out
#SBATCH --error=BlastX_kobu_feces.%J.err

module load cmake/3.20
module load compiler/gcc/10	
module load blast/2.10
module load biodata/1.0

cp -r Kobu_prot/ /scratch/
cp bat_feces_renamed_DEDUP.fasta /scratch/

blastx -word_size 10 -query /scratch/bat_feces_renamed_DEDUP.fasta -db /scratch/Kobu_prot/Kobu_prot_db -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out /scratch/Mada_Bat_Kobu_blastx_feces.txt -num_threads $SLURM_NTASKS_PER_NODE -evalue 1e-10

cp /scratch/Mada_Bat_Kobu_blastx_feces.txt .
```

For an alternative, and I think easier way to run BLAST searches using Bioconda, see ___. 

8. After the blast finishes, you'll want to curate a bit to the high quality hits. After Cara's lead, I went ahead and parsed for alignments that show alignment length > 100 aa and bit score > 100.

Here's the script for the nt parse (feces as example):

```
cat Mada_Bat_Kobu_blast_fecal_nt.txt | awk -F\t '($4>99 && $5<0.00001) {print $1,$3, $4, $5, $8,$9}' > Mada_Bat_Kobu_blast_fecal_nt_results.txt
```

And for the protein parse (feces again):

```
cat Mada_Bat_Kobu_blast_feces_prot.txt | awk -F\t '($4>99 && $6>99) {print $1,$3, $4,$5,$8,$9}' > Mada_Bat_Kobu_blast_feces_prot_results_100len100bit.txt
```

9.Once the blast results have beeen sub-selected a bit, you can summarize them to link back the hits to the samples of interest. Within the same folder as your output, try the following script to save the unique contig IDs which align to KobuV (example here for blastn alignment of fecal samples):

```
cat Mada_Bat_Kobu_blast_fecal_nt_results.txt | awk '{print $1}' | sort | uniq > Mada_Bat_Kobu_unique_contigs_feces_nt_hiqual.txt
```

And here to save the unique sample IDs for the same example (nucl):

```
cat Mada_Bat_Kobu_unique_contigs_feces_nt_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > Mada_Bat_CoV_unique_sampleID_feces_nt_hiqual.txt
```

And for the protein parse (feces again):

```
cat Mada_Bat_Kobu_blast_feces_prot_results_100len100bit.txt | awk '{print $1}' | sort | uniq > Mada_Bat_Kobu_unique_contigs_feces_prot_hiqual.txt
```

And here to save the unique sample IDs for the same example (prot):

```
cat Mada_Bat_Kobu_unique_contigs_feces_prot_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > Mada_Bat_Kobu_unique_sampleID_feces_prot_hiqual.txt
```

10.And do the same for the other sample types and for the blastx outputs. Load the outputs into R and determine which samples with meta-data were infected at various times/places. I put a bash script ('curate-blast-output.txt') in the 'blast-output' folder that runs through steps 6 through 8 for all of the blast output from this project. It saves summary files from both the curated sample set (those ending in "hiqual") and those from any kind of hit (lacking the suffix). You can run it with:

```
sh -e curate-blast-output.txt
```

You can see from investigating the "unique_contigs" folders that there is LOTS of cool stuff going on. We have a whole bunch of big contigs that are full genome of the virus(es):

-RR034B_010_NODE_1_length_29122 (Pteropus rufus, sample AMB130, 2/26/2018)
-RR034B_232_NODE_1_length_28926 (Rousettus madagascariensis, sample MIZ178, 4/14/2018)
-RR034B_288_NODE_2_length_28980 (Rousettus madagascariensis, sample MIZ240, 9/11/2018)

There are Eidolon hits in the urine, but none at full genome. The three above are HUGE contributions.I think the paper will have two trees: one small sequence RdRp tree of our samples and closely related CoVs and one full genome phylogeny of all CoVs that includes the three megasamples above. #Cara's comment 

11.Now, import the curated contigs into R and compare them against the IDseq hits for what is KobuV positive and how it maps into the meta-data. It looks like no throat samples were CoV hits, as is consistent with what is recovered from IDseq. See R script ('CoV-hits-vs-metadata.R') for further comparison of manual pipeline hits for fecal and urine samples. There are differences, so I am going through them individually IDseq to check. #Cara's comment, update later

12.For calling positives in cases where there was a discrepancy between this (stringent) pipeline and IDseq (see spreadsheet here), we will accept them as positive hits if (and only if!) the reads from that sample assembled into one or more contigs. In this case, contigs should only be acceptable if the average read depth at that contig is 2 or more reads (per Amy's rule). So, in manually curating any positive samples from IDseq, check the broad (not hiqual) contig summary file for that sample (for feces, "20210721_Mada_Bat_CoV_unique_contigs_feces_nt.txt") and only call it as positive if it has at least one contig with >2 reads for average coverage.
