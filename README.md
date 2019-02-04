# Selection Pipeline

These scripts make up a pipeline used for running tests for signatures of selection in comparative genomics studies. Filtering and alignment of coding sequences are the most basic features. These alignments can be fed into PAML and HyPhy to more easily run tests across all genes present in a set of genomes.

### Alignment

#### Required inputs

It is necessary that species designations are four characters long and that they preface gene names.

1. A file with orthogroups. Theoretically, any method can be used to identify orthogroups. Currently, the pipeline is only fully functional with the output format of OrthoFinder, e.g.:
OG0000000: TAX1_00001 TAX2_00002

If converting from another format, it is required to maintain at least five zeros in the OG designation.

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

This command will create a new directory, `[output_prefix]_fsa_coding` with unaligned (`*.fa`) and aligned (`*.afa`) fasta files for every locus.

Finally, we want to filter the alignments, removing both low-confidence positions and low confidence or low information sequences.

Post-alignment, a number of additional filters can be applied using the `alignment_filter` action. If you run this action then several steps will be executed by default:
1. Columns in the alignment with fewer than a minimum number of non-gap characters can be removed. The default is 8 and this can be adjusted with `--nogap_min_count`.
2. Columns in the alignment with less than a proportion of sequences with known nucleotides (anything other than "-" and "N") can be removed. Default behavior is to require at least 30% of sequences be non-gap characters but this can be adjusted with ```--nogap_min_prop```.
3. Columns with sequences from fewer than a minimum number of species can be removed. Default behavior requires at least 4 species have non-gap characters but this can be changed with `--nogap_min_species`.
4. Then we used the Jarvis et al. (Science 346: 1320-1331) Avian Phylogenomics Project scripts for filtering amino acid alignments, which mask over poorly aligning regions of individual sequences (rather than omitting entire alignment columns). The scripts ```spotProblematicSeqsModules.py``` and ```spotProblematicSeqsModules-W12S4.py``` were downloaded from ftp://parrot.genomics.cn/gigadb/pub/10.5524/101001_102000/101041/Scripts.tar.gz on 1/31/2019. These scripts were incorporated into the ```selection_pipeline.py``` pipeline through the ```jarvis_filtering()``` function.
5. These masked sequences are then run through the first three filters again.
6. If the number of known nucleotides in a sequence after filtering was less than half the length of that sequence originally prefiltering.
7. If the sequence was more than 50% unknown (gap or masked) sequence.
8. Sequences with fewer than 300 known nucleotides in a sequence after filtering.

Several new directories are created during the filtering process. The aligned fasta files from the first three filtering steps are including in the `[output_prefix]_fsa_coding_columnfilt` directory. The output from the Jarvis filter in step 4 is in '[output_prefix]_coding_jarvis` and the column-filtered Jarvis-masked alignments (step 5) are in '[output_prefix]_coding_jarvis_columnfilt`. The results from steps 6-8 are in '[output_prefix]_coding_jarvis_columnfilt_seqfilt`. The Jarvis filter can take awhile so you will want to use many threads for this step as well.