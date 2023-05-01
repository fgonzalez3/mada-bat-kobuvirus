# Building A Phylogenetic Using BEAST2

# Gathering Background Sequences From GenBank

The goal of this analysis was to estimate the timing of the most recent common ancestor (MRCA) of all Kobuviruses, as well as to estimate the time of divergence of the Madagascar bat Kobuvirus from other Aichivirus A lineages defined in our paper. To do this, we compiled Kobuviruses from each Aichivirus sub-clade, according to their known hosts. Reference sequences were prioritized, though we opted for non-reference sequences if none existed. We also included two additional Picornaviridae sequences found in Fig2B-2C. A summary table can be found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/blob/main/Fig3/beast_kobuvirus_metadata_manual%20copy.csv). To estimate the MRCA of all Kobuviruses, we opted to construct a Bayesian TimeTree.

# Alignment and Model Selection

After compiling sequences for this tree, we renamed the sequences in BEAST format. See [script](https://github.com/fgonzalez3/mada-bat-kobuvirus/blob/main/TreePrep/pre_beast_name.R). This included the accession number and the collection date. For cases where onlny a collection year was reported, we set the middle of the year (July 15th) as the corresponding date. We aligned these sequences using [MAFFT](https://mafft.cbrc.jp/alignment/server/) in GeneiousPrime under default parameters and evaluated the optimal nucleotide substitution model using [ModelTest-NG](https://github.com/ddarriba/modeltest). This corresponded to TVM+I+G4 for the dataset in hand. 

The MSA and ModelTest-NG outputs can all be found [here](https://github.com/fgonzalez3/mada-bat-kobuvirus/tree/main/Fig3). 

# Building a Phylogenetic Tree Using BEAST2

The BEAST and BEAST2 comunity maintains a number of extremely helpful tutorials that you should practice before getting BEAST2 running. Found [here](https://taming-the-beast.org/tutorials/Introduction-to-BEAST2/). 

After model selection, we built a Bayesian phylogenetic tree for the Kobuvirus phylogeny in BEAST2, using nucleotide substitution model suggested by ModelTest-NG above. The first step in this process requires generation of a .xml file for input into BEAST, using the program BEAUTi. Some substitution models are easily specified in BEAUTi. For specifying less common models in BEAUTi, see this [blog post](https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/), the reccomendations [here](https://groups.google.com/g/ggplot2/c/H50aGubqt2U), or [here](http://www.iqtree.org/doc/Substitution-Models). It is also possible to generate .xml files outside of BEAUTi, though this approach was not needed for this scenario. 

To prepare the .xml file, we used the following parameters in the tab inputs at the top of the screen in BEAUTi:

- **Tip Dates**: We used the date of the sample collection as the "Tip Date". For any sample from GenBank which only listed year of collection, we set the tip date to July 15 of the year of collection. Each alignment was uplodaed to BEAUTi with sequence names arranged so as to easily decipher the date. 

