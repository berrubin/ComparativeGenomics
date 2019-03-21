# Selection Pipeline

These scripts make up a pipeline used for running tests for signatures of selection in comparative genomics studies. Filtering and alignment of coding sequences are the most basic features. These alignments can be fed into PAML and HyPhy to more easily run tests across all genes present in a set of genomes.

### Alignment

#### Required inputs

It is necessary that species designations are four characters long and that they preface gene names.

1. A file with orthogroups. Theoretically, any method can be used to identify orthogroups. Currently, the default is to read files of the output format of OrthoFinder, e.g.:
OG0000000: TAX1_00001 TAX2_00002

However, OrthoMCL and ProteinOrtho can also be read if `-w orthomcl` or `-w proteinortho` is specified when calling `selection_pipeline.py`. If converting from another format, it is required to maintain at least five zeros in the OG designation.

2. Files of CDS for each gene included in the orthogroups file. These files are listed in a two-column params file where the first column is the four character taxon abbreviation and the second column is the path to the fasta file of CDS for that taxon.

#### Running alignment

1. Write fastas for the orthogroups that meet the filtering requirements:

```selection_pipeline.py -a write_orthos -b [output directory] -o [output prefix] -r [orthogroup file] -t [min. taxa to include locus] -d [params file]```

`-a` is the "action" to perform. In this case, we want to `write_orthos`.

`-b` is the base output directory where all output will be written.

`-o` is the prefix to be used for labeling output for this project.

`-r` is the path to the OrthoFinder formatted orthogroup file.

`-t` is the minimum number of taxa represented in an orthogroup in order to include that locus.

`-d` is the path to the params file listing the locations of the CDS files for each species. The format is given above.

Several filtering steps are hardcoded in at this stage:
1. Any sequence less than half the median length of all sequences in an orthogroup is removed.
2. Any species with more than three sequences in an orthogroup is completely removed.
3. Orthogroups with more than 1.5x the number of sequences as the number of species in an orthogroup is completely discarded.

These can be modified by changing the `remove_shorter_seqs()` function in `utils.py` if you want to.

Within the base output directory, the result of this command is a directory ```[output_prefix]_orthos``` with fasta files of the CDS for every orthogroup and a directory ```[output_prefix]_orthos_prots``` with fasta files of the amino acid sequences for every orthogroup. An index file is also created ```[output_prefix]_ortho.index``` that lists the numbers of taxa and numbers of sequences present in each orthogroup. This file acts as a reference for future processing and choosing of orthogroups.

In order to align the written orthogroups, the pipeline uses `FSA` with the `--nucprot` option to align coding sequence. In order to run the alignment, a command such as:

```selection_pipeline.py -a align_coding -p 16 -b [output directory] -o [output prefix] -r [orthogroup file] -t [min. taxa to include locus] -d [params file]```

Note that we use the additional `-p` parameter here to specify the number of threads to use. The alignment can take awhile so using multiple threads is a good idea. It is embarrassingly parallelized so just splits all of the orthogroups among the available threads and aligns each one on a single thread.

This command will create a new directory, `[output_prefix]_fsa_coding` with unaligned (`*.fa`) and aligned (`*.afa`) fasta files for every locus. If the pipeline crashes before finishing all of the alignments, not to worry. Simply rerun the above command. It checks whether alignments for each OG have already been completed and will not redo anything that is already done.

An aligned concatenated matrix of proteins from all orthogroups with a single sequence from each species will also be produced in `[output_prefix].afa`. These sequences are from the `trimAl` filtered alignments -- more conservative than the Jarvis filtering but less conservative than `Gblocks`. This matrix can be used to infer a phylogeny if you like. It is formatted to work with RAxML (i.e., `raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAWAG -s [output_prefix].afa -n [output_prefix].tree -T 20`).

### Alignment filtering

Finally, we want to filter the alignments, removing both low-confidence positions and low confidence or low information sequences. By default, both Gblocks and trimAl are run on each OG and output is included. However, these two methods can be quite stringent, removing potentially valuable informative sequence. Therefore, we use several other post-alignment filters. These can be applied using the `alignment_filter` action. If you run this action then several steps will be executed by default:
1. Columns in the alignment with fewer than a minimum number of non-gap characters can be removed. The default is 8 and this can be adjusted with `--nogap_min_count`.
2. Columns in the alignment with less than a proportion of sequences with known nucleotides (anything other than "-" and "N") can be removed. Default behavior is to require at least 30% of sequences be non-gap characters but this can be adjusted with ```--nogap_min_prop```.
3. Columns with sequences from fewer than a minimum number of species can be removed. Default behavior requires at least 4 species have non-gap characters but this can be changed with `--nogap_min_species`.
4. Then we used the Jarvis et al. (Science 346: 1320-1331) Avian Phylogenomics Project scripts for filtering amino acid alignments, which mask over poorly aligning regions of individual sequences (rather than omitting entire alignment columns). The scripts ```spotProblematicSeqsModules.py``` and ```spotProblematicSeqsModules-W12S4.py``` were downloaded from ftp://parrot.genomics.cn/gigadb/pub/10.5524/101001_102000/101041/Scripts.tar.gz on 1/31/2019. These scripts were incorporated into the ```selection_pipeline.py``` pipeline through the ```jarvis_filtering()``` function.
5. These masked sequences are then run through the first three filters again.
6. Then, sequences are removed completely if the number of known nucleotides in a sequence after filtering was less than half the length of that sequence originally prefiltering. This ratio can be changed with `--min_seq_prop_kept`.
7. Sequences are also removed if the sequence was more than 50% unknown (gap or masked) sequence. This can be changed with `--max_seq_prop_gap`.
8. Finally, sequences with fewer than 300 known nucleotides in a sequence after filtering are removed. This minimum length can be changed with `--min_cds_len`.

Several new directories are created during the filtering process. The aligned fasta files from the first three filtering steps are including in the `[output_prefix]_fsa_coding_columnfilt` directory. The output from the Jarvis filter in step 4 is in `[output_prefix]_coding_jarvis` and the column-filtered Jarvis-masked alignments (step 5) are in `[output_prefix]_coding_jarvis_columnfilt`. The results from steps 6-8 are in `[output_prefix]_coding_jarvis_columnfilt_seqfilt`. The Jarvis filter can take awhile so you will want to use many threads for this step as well.

### Compiling data for RERconverge

`RERconverge` is a relatively new piece of software that is built to find signatures of convergent rate shifts in genes on a genome-wide scale. It has been and continues to be developed here: https://github.com/nclark-lab/RERconverge. The required input for running RERconverge is a set of phylogenies, one from each gene of interest, with branch lengths inferred from protein sequence divergence. `selection_pipeline.py` can produce such a file using AAML to infer branch lengths using the `-a rer_converge` action command. 

This will require as input a species tree, which can be specified using `-e`. Note that branch lengths don't matter for this species tree as it only provides the topology for branch lengths to be estimated on for each gene.

If you need to filter loci to make sure that particular taxa are present or need to exclude particular taxa from an analysis, you can use the `--taxa_inclusion` parameter. Here, you should specify the path to a tab-delimited file with your taxonomic requirements. There are a few options for formatting.

Below would require that at least 1 of HLIG, HQUA, and HRUB as well as at least 2 of AAUR, APUR, and AVIR be present in a locus in order for it to be included in the analyis:
`1    HLIG,HQUA,HRUB
2     AAUR,APUR,AVIR`

Sometimes we need particular sets of taxa to be present together. For example, below would require that at least one of the pairs specified (i.e., LMAR and LFIG, LZEP and LVIE, or LPAU and LOEN) be present for a locus to be included. A locus with LMAR, LFIG, and LZEP would be included but a locus with LMAR, LZEP, and LPAU, would not.
`1	(LMAR,LFIG),(LZEP,LVIE),(LPAU,LOEN)`

You can also exclude taxa from the analysis. Below would remove HLIG from all loci before running the analysis:
`-1 	HLIG`

The `--taxa_inclusion` file can have as many lines and specifications as you like. 

Not that RERconverge works on the species levels so requires only a single sequence per species per locus. Duplicated genes or orthologous groups that otherwise have more than one sequence per species can't be included. Therefore, all such sequences are removed for this analysis. So if multiple sequences are present for a particular species, all of the sequences for that species are removed before proceeding with the analysis.

Finally, you will also need to specify an output file name using `--outputfile`. The file will be placed in a predetermined directory so do not specify the whole path, just the name of the file itself.

So a command for producing input for `RERconverge` could look like this:

`selection_pipeline.py -a rer_converge -p 16 -b [output directory] -o [output prefix] -t [min. taxa to include locus] --outputfile [output file name] --taxa_inclusion [taxon requirement file] -e [newick species tree file]`

This will create a directory named `aaml_compiled/` in your base directory with the output file in it. This file will have the locus name as the first column and the species tree with branch lengths for that locus in the second column and can be passed directly to RERconverge.

### GO enrichment

I run a lot of GO enrichment analyses so incorporated the GOATOOLS (https://github.com/tanghaibao/goatools) `find_enrichment.py` script into `selection_pipeline.py`. To run this type of enrichment, you will need three files. First, you need GO annotations for all of your gene orthogroups. These need to be in a tab-delimited 2-column format where the first column is the orthogroup ID and the second column is a semicolon separated list of GO terms, i.e.,:
`OG_10001   GO:0003674;GO:0065007;GO:2001141
OG_10002    GO:0000226`

Then you will need two files, one listing the set of orthogroup IDs that are of interest and the other representing the background set of orthogroup IDs to which the first set should be compared. GOATOOLS performs many tests by default but I trim these down to make the multiple test correction a bit more manageable. `selection_pipeline.py` only includes Biological Process terms and only tests those terms present at least 3 times in the target set of genes. 

`selection_pipeline.py -a goatools -b [output directory] -o [output prefix] --go_database [GO database file] --goa_forefile [file with target set of orthogroups] --goa_backfile [file with background set of orthogroups]`

This will produce a `[base directory]/[output prefix]_goatools/` directory. In that directory, another directory will be created with the name of `--go_forefile` without the path or file extension. GOATOOLS will do its work in that directory and produce an output file with the same name as the created directory plus `_ps.txt`. The file with that base name but ending in `_filtered_ps.txt` is the processed result excluding rare GO terms and non-Biological Process terms. The columns in the output file should be relatively self-explanatory.

### GO enrichment in RERconverge output

I look for GO enrichment in RERconverge output files a lot so I included a function to take in this type of output directly. RERconverge finds genes evolving both significantly faster and significantly slower in the taxa of interest and `selection_pipeline.py` examines these sets of genes separately for the purposes of GO term enrichment analysis. The command to use is `-a rer_goatools`, i.e.,:

`selection_pipeline.py -a rer_goatools -b [output directory] -o [output prefix] --go_database [GO database file] -j [RERconverge output file]`

`-j` specifies the file produced by RERconverge. Running this command will produce a directory `[base directory]/RER_goatools`. Inside this directory, two directories will be produced with the format `rer_0.05_[slow/fast]er_go_[RERconverge output rootname]`. The RERconverge output rootname is just the name of the RERconverge output file without the file extension. "slow" and "fast" indicate whether the foreground set of genes was found to be evolving significantly faster or significantly slower with a p-value cutoff of 0.05. `_filtered_ps.txt` files will be produced for each foreground set as explained above in the GO enrichment section.

### What about gene trees?

Sometimes, gene trees are different than species trees. This could cause problems when estimating branch lengths so one may want to be aware of possible discordance between gene trees and species trees so I built this into `selection_pipeline.py`. First, we need to build gene trees using the `-a nopara_gene_trees`. Since we are looking for gene tree/species tree discordance we again remove sequences from taxa represented more than once in each locus.

`selection_pipeline.py -a nopara_gene_trees -p 16 -b [output directory] -o [output prefix] -t [min. taxa to include locus]`

This command will make gene trees from the nucleotide sequences of each locus using RAxML with the GTRGAMMA model.

Then, we use the `-a check_discordance` action command to check for signatures of discordance, i.e.,:

`selection_pipeline.py -a check_discordance -p 16 -b [output directory] -o [output prefix] -t [min. taxa to include locus] -e [newick species tree file]`

Here, we need to also include the species tree (`-e`) that we are comparing the gene trees to. This will use `FastTree` to compare the likelihoods of the gene tree and species tree topologies given the nucleotide sequence alignment for that locus and uses `CONSEL` to determine if there is significant evidence of discordance. Individual output files for each locus will be put in `[base directory]/[output prefix]_discordance/`. A summary file will also be produced at `[base directory]/consel_consistency.txt` with p-values indicating the significance of the differences between the species tree and gene tree topologies. FDR-corrected p-values are also included in this table, as well as the number of taxa in each locus and the number of bases in the alignment for each locus. To be clear, low p-values mean a higher chance of discordance.

### HyPhy RELAX

`selection_pipeline.py` also includes an implementation of the HyPhy RELAX test for relaxation or intensification of selection in target sets of taxa. The command for running HyPhy RELAX is of a similar form as for RERconverge. `--taxa_inclusion` should be used to indicate minimum taxon representation and those taxa to exclude and a species tree needs to be included with `-e`. In addition, a comma-separated list (i.e., HLIG,HQUA,AAUR) of taxa to be treated as the target taxa needs to be included with `-c`. Relaxation of selection will be identified in the specified taxa. I have also implemented this to work on non-terminal lineages but I'm not going to bother explaining that right now because it will probably never get used.

`selection_pipeline.py -a hyphy_relax -p 16 -b [output directory] -o [output prefix] -t [min. taxa to include locus] --taxa_inclusion [taxon requirement file] -e [newick species tree file] -c [comma-separated list of focal taxa]`

This will create a directory at `[base directory]/[output prefix]_[focal taxa list]_RELAX`. All files required to run HyPhy RELAX and all results will be created in this directory. Output files will be in `og_[current OG]_relax_unlabeledback.txt`. A summary file will also be produced at `[base directory]/compiled_relax_[focal taxa list].txt`. This file will list the OG, the K parameter estimated by HyPhy RELAX, the p-value inferred by HyPhy RELAX, and an FDR-corrected p-value. The "interpretation" column indicates whether a relaxation of selection or intensification of selection is indicated (K < 1 indicated relaxation and K > 1 indicates intensification).

HyPhy RELAX takes substantial computation time. If this process crashes or is cancelled for some reason, it can be restarted with the same command and work will not be repeated.

## Population genomics

MK-test functionality is built into selection_pipeline.py. Several filters are built into the reading of the genotypes provided. First, if a site has genotypes for fewer than 4 individuals, it is not included in the analysis. Second, if all individuals have heterozygous genotypes or are homozygous for an alternative allele, that site is not included. Incomplete genotype calls with only a single allele from an individual called are also excluded before making the determination that all individuals are either homozygous for an alternative allele or heterozygous. So an individual with "0/." is not counted as either heterozygous or homozygous for the reference allele. Anything with a minor allele frequency of less than or equal to 0.1 is not counted as variable as has been done previously (http://www.genetics.org/content/158/3/1227.long, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3666677/). Also, any alternative allele that only appears once is not considered. However, if just a single individual is homozygous for an alternative allele that is considered.

### Non-coding sequence

NCARs (from https://github.com/berrubin/BeeGenomeAligner) are Non-Coding Alignable Regions. While not statistically conserved, they are conserved enough to be alignable across your taxa of interest. I have built functionality into `selection_pipeline.py` so that it is relatively straightforward to perform tests for selection on NCARs. However, any non-coding sequence could be examined in a similar way if you have GFFs with their coordinates, their sequences in fasta format, and the genome sequences of the taxa being examined. An index file, very similar to the file indicating orthologous groups of genes is also required.

The approach for examining selective pressures on non-coding sequences is simply an MK-test. All bases in the NCAR are considered to be selected sites. Neutral sites are the fourfold degenerate sites in the 20 genes neighboring the NCAR. Both alpha and the neutrality index are calculated. The selected sites are also tested for a difference to the neutral sites using Fisher exact tests.

