### Finding initial positive Kobuvirus samples after mNGS + IDseq

1.First, I downloaded all non-host contigs derived from fecal samples, mapping to ANY taxon from IDseq in the 'RR034B1_feces' project. (Note that we did not take contigs from HeLa controls or water). This is done in bulk, manually, on IDseq.net in the top right-hand corner. 

Note that when you download all the non-host contigs, it will produce a folder with a separate fasta file for each sample, which lists the contigs by node number but does not include the sample ID. Before joining all the contigs (nodes) together, you need to distinguish them by sample ID. Cara wrote an Rscript that parses this for each filetype (rename-fastas-feces, rename-fastas-urine, rename-fastas-throat). To rename your files and the headers within them, copy the appropriate Rscript for the tissue type into your appropriate downloads folder, cd into that folder on the command line, and simply type (example here for feces):

*NOTE: Rscript got hung up on contig file 193, so it was omitted. 

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

8. Once BLAST finishes, create a broad catch-all summary file from the BLAST outputs. There are too many hits, but you'll want to keep them as references. We'll report the number of initial hits on our publication. 

Contigs:

```
cat initial_KoV_blast_nt.txt | awk '{print $1}' | sort | uniq > KoV_unique_contigs_feces_nt.txt
```

Sample IDs:

```
cat initial_KoV_blast_nt.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > KoV_unique_sampleID_feces_nt.txt
```

9. I parsed for high quality hits. After Cara's lead, I went ahead and parsed for alignments that show alignment length >100 + e-value <0.00001 and alignment length >100 + bit score >100. Scripts below:

Script for nt parse using e-value:

```
cat initial_KoV_blast_nt.txt | awk -F\t '($4>99 && $5<0.00001)' > KoV_blast_nt_100len5eval.txt
```


Script for nt parse using bit score: 

```
cat initial_KoV_blast_nt.txt | awk -F\t '($4>99 && $6>99)' > KoV_blast_nt_100len100bit.txt
```


Script for aa parse using e-value: 

```
cat initial_KoV_blast_aa.txt | awk -F\t '($4>99 && $5<0.00001)' > KoV_blast_nt_100len5eval.txt
```

Script for aa parse using bit score: 

```
cat initial_KoV_blast_aa.txt | awk -F\t '($4>99 && $6>99)' > KoV_blast_nt_100len100bit.txt
```


10.Once the blast results have beeen sub-selected a bit, you can summarize them to link back the hits to the samples of interest. Within the same folder as your output, try the following script to save the unique contig IDs which align to Kobuviruses.

Example here shows the nt parse for e-value as query:

```
cat KoV_blast_nt_100len5eval.txt | awk '{print $1}' | sort | uniq > KoV_unique_contigs_feces_nt_hiqual.txt
```

And here to save the unique sample IDs for the same example:

```
cat KoV_unique_contigs_feces_nt_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > KoV_unique_sampleID_feces_nt_hiqual.txt
```

I did the same for the other three parses, giving me high-quality sample IDs for:

- nt parse (100 aln & 100 bit)
- nt parse (100 aln & 5 eval)
- aa parse (100 aln & 100 bit)
- aa parse (100 aln & 5 eval)

11. Do the same for other sample types if you are studying others. Load the outputs into R and determine which samples with meta-data were infected at various times/places. Cara wrote a bash script ('curate-blast-output.txt'), attached in the 'blast-output' folder, that runs through through some of the above steps. It saves summary files from both the curated sample set (those ending in "hiqual") and those from any kind of hit (lacking the suffix). You can run it with:

```
sh -e curate-blast-output.txt
```

After investigating initial BLASTn and BLASTx results, there's some good data. Two contigs listed here are common hits across all four of the above parses. This is where we find our full genome, NODE_4_length_8263 (NCBI Accession: OP287812). I used it as a reference to run an additional offline BLAST upon its acceptance into NCBI, from which I made Figure 1 and determined true positives. 

- RR034B_052_NODE_4_length_2077_cov_943.883284
- RR034B_145_NODE_4_length_8263_cov_106.105872

Since these contigs are the only two that appeared in all four of my BLAST analyses, listed at the end of step 10, I did not choose to look at any other contigs further. I considered the two above contigs as my "initial positives". You can find a summary of this analysis in Supplementary Figure __. 

12. My consensus-based method differs from Cara's positive calling. Considering there are few Kobuviruses on NCBI, I did not want to risk calling false-positives. 
