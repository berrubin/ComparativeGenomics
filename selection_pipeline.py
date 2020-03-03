#!/usr/local/bin/python
import os
import paml_tests
import utils
from optparse import OptionParser
import sys

parser = OptionParser()

parser.add_option("-p", "--num_threads", dest = "num_threads", type = int, default = 1, help = "Number of cores to use.")
parser.add_option("-m", "--min_og_group", dest = "min_og_group", type = int, default = 0, help = "For limiting the number of OG's examined. Don't analyze those with OG numbers less than this number.")
parser.add_option("-x", "--max_og_group", dest = "max_og_group", type = int, default = 10000000, help = "For limiting the number of OG's examined. Don't analyze those with OG numbers more than this number.")
parser.add_option("-o", "--prefix", dest = "prefix", type = str, help = "String used at the beginning of output directories and files.")
parser.add_option("-b", "--base_dir", dest = "base_dir", type = str, help = "Output directory.")
parser.add_option("-t", "--min_taxa", dest = "min_taxa", type = int, default = 4)
parser.add_option("-r", "--ortho_file", dest = "ortho_file", type = str, default = "Orthogroups.txt", help = "File of orthologous groups.")
parser.add_option("--ncar_ortho_file", dest = "ncar_ortho_file", type = str, help = "File of orthologous NCARs.")
parser.add_option("-e", "--tree_file", dest = "tree_file", type = str, default = "species.tree", help = "Phylogeny of species examined.")
parser.add_option("-a", "--action", dest = "action", type = str, help = "Analysis to do", default = "paml")
parser.add_option("-i", "--inspecies", dest = "inspecies", type = str, help = "Species of interest in pairwise analyses")
parser.add_option("-u", "--outspecies", dest = "outspecies", type = str, help = "Outgroup species in pairwise analyses")
parser.add_option("-y", "--timetree", dest = "timetree", type = str, help = "Time calibrated tree")
parser.add_option("-z", "--traittree", dest = "traittree", type = str, help = "Trait tree")
parser.add_option("-f", "--noncoding_seq_type", dest = "noncoding_seq_type", type = str, help = "Sequence type for HKA tests (flank, intron, or first_intron)", default = "flank")
parser.add_option("-g", "--no_gblocks", dest = "use_gblocks", action = "store_false", default = False, help = "Should Gblocks be used on alignments?")
parser.add_option("--no_paralogs", dest = "no_paralogs", action = "store_true", default = False, help = "Should all paralogs be discarded immediately?")
parser.add_option("-c", "--foreground", dest = "foreground", type = str, default = "", help = "Foreground taxa for selection test")
parser.add_option("-d", "--param_file", dest = "param_file", type = str, default = "halictids.params")
parser.add_option("--gff_params", dest = "gff_params", type = str, default = "halictids_gff.params", help = "File with paths to coding GFF files")
parser.add_option("--ncar_gff_params", dest = "ncar_gff_params", type = str, default = "halictids_ncar_gff.params", help = "File with paths to ncar GFF files")
parser.add_option("--genome_params", dest = "genome_params", type = str, default = "halictid_genomes.params", help = "File with paths to genome files")
parser.add_option("-j", "--rerconverge_output", dest = "rerconverge_output", type = "str", default = "", help = "Output file from RERconverge")
parser.add_option("-w", "--orthofile_format", dest = "orthofile_format", type = "str", default = "orthofinder", help = "Format of orthofile.")
parser.add_option("--nogap_min_count", dest = "nogap_min_count", type = int, default = 8, help = "The minimum number of sequences in an alignment that are not gaps in an individual column.")
parser.add_option("--nogap_min_prop", dest = "nogap_min_prop", type = float, default = 0.3, help = "The minimum proportion of sequences in an alignment that are not gaps in an individual column.")
parser.add_option("--nogap_min_species", dest = "nogap_min_species", type = int, default = 4, help = "The minimum number of species in an alignment that are not gaps in an individual column.")
parser.add_option("--min_seq_prop_kept", dest = "min_seq_prop_kept", type = float, default = 0.5, help = "The minimum fraction of a sequence that has to remain after every filtering step in order to keep that sequence in the alignment.")
parser.add_option("--max_seq_prop_gap", dest = "max_seq_prop_gap", type = float, default = 0.5, help = "The maximum fraction of a sequence that is gap characters. Otherwise that sequence is discarded from the alignment.")
parser.add_option("--min_cds_len", dest = "min_cds_len", type = int, default = 300, help = "The minimum length of a coding sequence to keep.")
parser.add_option("--paths_file", dest = "paths_file", type = str, default = "pathsfile.params", help = "File containing paths to executables.")
parser.add_option("--outputfile", dest = "outputfile", type = str, default = "output.txt", help = "Name of output file.")
parser.add_option("--taxa_inclusion", dest = "taxa_inclusion", type = str, help = "File with taxa requirements.")
parser.add_option("--go_database", dest = "go_database", type = str, help = "File with GO terms mapped to orthogroup names")
parser.add_option("--goa_forefile", dest = "goa_forefile", type = str, help = "File listing focal orthogroups to be tested against the background for GO enrichment")
parser.add_option("--goa_backfile", dest = "goa_backfile", type = str, help = "File listing background orthogroups to be used for GO enrichment")
parser.add_option("--pickle_dir", dest = "pickle_dir", type = str, help = "Directory for storing gff/vcf data pickles")
parser.add_option("--hyper_pop", dest = "hyper_pop", type = str, help = "All OGs under consideration (e.g. all OGs with human orthologs).")
parser.add_option("--hyper_pop_cond", dest = "hyper_pop_cond", type = str, help = "All OGs with condition (e.g. all OGs with human orthologs with autism association")
parser.add_option("--hyper_targets", dest = "hyper_targets", type = str, help = "Focal OGs (e.g. OGs evolving faster in social taxa).")
parser.add_option("--hyper_targets_back", dest = "hyper_targets_back", type = str, help = "Focal OGs background (e.g. OGs included in test of rate changes).")
parser.add_option("--og_list_file", dest = "og_list_file", type = str, help = "List of OGs to use for analyses")

(options, args) = parser.parse_args()

#PATHS_DIC = utils.read_exec_paths(options.paths_file)

HALICTUS = ["HRUB", "HQUA", "HLIG"]
AUGOCHLORINES = ["AAUR", "APUR", "MGEN"]
LASIOGLOSSUM = ["LLEU", "LMAR", "LFIG", "LZEP", "LVIE", "LALB", "LCAL", "LMAL", "LPAU", "LOEN"]

def main():
    if not os.path.isdir(options.base_dir):
        os.mkdir(options.base_dir) #create working directory
    print options.action

    if options.action == "write_orthos":
        cds_dic = utils.read_params(options.param_file)
        ortho_dic = utils.read_orthofile(options.orthofile_format, options.ortho_file)
        index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
        seq_dic = utils.get_cds_files(cds_dic)
        if options.no_paralogs:
            utils.write_orthos(options.ortho_file, seq_dic, "%s/%s_orthos" % (options.base_dir, options.prefix), index_file)
        else:
            utils.write_orthoparagroups(ortho_dic, seq_dic, "%s/%s_orthos" % (options.base_dir, options.prefix), index_file, options.min_taxa)
        print "Orthogroups written to %s/%s_orthos" % (options.base_dir, options.prefix)
        print "Exiting"
        sys.exit()

    if options.action == "write_ncars":
        ncar_dic = utils.read_params(options.param_file)
        ortho_dic = utils.ncar_ortho_dic(options.ncar_ortho_file, options.min_taxa) #this needs to be "filtered_loci.index" from the NCAR pipeline
        seq_dic = utils.get_cds_files(ncar_dic)
        index_file = "%s/%s_ncar_ortho.index" % (options.base_dir, options.prefix)
        utils.write_ncars(ortho_dic, seq_dic, "%s/%s_ncars" % (options.base_dir, options.prefix), options.min_taxa, index_file)
        sys.exit()

    if options.action == "write_cnees":
        ncar_dic = utils.read_params(options.param_file)
        ortho_dic = utils.ncar_ortho_dic(options.ncar_ortho_file, options.min_taxa) #this needs to be "filtered_loci.index" from the NCAR pipeline
        seq_dic = utils.get_cds_files(ncar_dic)
        index_file = "%s/%s_ncar_ortho.index" % (options.base_dir, options.prefix)
        utils.write_ncar_cnees(ortho_dic, seq_dic, "%s/%s_ncars" % (options.base_dir, options.prefix), options.min_taxa, index_file)
        sys.exit()


###Align coding sequences and concatenate all protein sequences into 
###an aligned matrix that can be input into RAxML to make a phylogeny.
    if options.action == "align_coding":
        cds_dic = utils.read_params(options.param_file)
        index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
        paras_allowed = True 
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed) #Gets list of OGs that meet minimum taxa requirement. If paras_allowed is False then will not return any OGs with any paralogs in them.
        iscoding = True
        utils.fsa_coding_align(og_list, "%s/%s_orthos/" % (options.base_dir, options.prefix), "%s/%s_fsa_coding" % (options.base_dir, options.prefix), options.num_threads, iscoding)
        print "Orthogroups aligned using FSA and output written to %s/%s_fsa_coding" % (options.base_dir, options.prefix)
        paras_allowed = False
        og_list = utils.read_ortho_index(index_file, len(cds_dic.keys()), paras_allowed) #Gets only those OGs that have a single sequence for every species in the study. This is for making a sequence matrix that can be used for phylogenetics.
        utils.concatenate_for_raxml("%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s.afa" % (options.base_dir, options.prefix), og_list, cds_dic.keys())
        print "If you would like to run a phylogenetic analysis, a concatenated amino acid sequence matrix of all orthogroups including all %s of the species in your study has been written to %s/%s.afa" % (len(cds_dic.keys()), options.base_dir, options.prefix)
        print "Exiting"
        sys.exit()

    if options.action == "fourfold_matrix":
        cds_dic = utils.read_params(options.param_file)
        index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
        paras_allowed = False
        og_list = utils.read_ortho_index(index_file, len(cds_dic.keys()), paras_allowed) #Gets only those OGs that have a single sequence for every species in the study. This is for making a sequence matrix that can be used for phylogenetics.
        utils.concatenate_fourf_for_raxml("%s/%s_gene_ancestral" % (options.base_dir, options.prefix), "%s/%s_fourfold.afa" % (options.base_dir, options.prefix), og_list, cds_dic.keys())
        sys.exit()

    if options.action == "align_ncars":
        iscoding = False
        ortho_dic = utils.ncar_ortho_dic(options.ncar_ortho_file, options.min_taxa) #this needs to be "filtered_loci.index" from the NCAR pipeline
        ncar_list = ortho_dic.keys()
        utils.fsa_ncar_align(ncar_list, "%s/%s_ncars" % (options.base_dir, options.prefix), "%s/%s_fsa_ncar" % (options.base_dir, options.prefix), options.num_threads, iscoding)
        sys.exit()
        
    if options.action == "ncar_ancestor":
        test_type = "ancestral"
        foreground = "ncar"
        ortho_dic = utils.ncar_ortho_dic(options.ncar_ortho_file, options.min_taxa)
        ncar_list = ortho_dic.keys()
#        ncar_list = [10510, 10513, 17170, 17388, 17389, 23400]
        utils.paml_test(ncar_list, foreground, test_type,"%s/%s_fsa_ncar" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads, options.use_gblocks, options.min_taxa, [])
        sys.exit()
    
    if options.action == "coding_ancestor":
        test_type = "ancestral"
        foreground = "gene"
        index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
        paras_allowed = True 
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        utils.remove_aligned_paras(og_list, "%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_fsa_coding_noparas" % (options.base_dir, options.prefix), "%s/%s_ortho_noparas.index" % (options.base_dir, options.prefix))
        index_file = "%s/%s_ortho_noparas.index" % (options.base_dir, options.prefix)
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding_noparas" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads, options.use_gblocks, options.min_taxa, [])
        sys.exit()

    if options.action == "ncar_mk":
        ncar_gff_dic = utils.read_params(options.ncar_gff_params)
        ncar_gff_file = ncar_gff_dic[options.inspecies]
        coding_gff_dic = utils.read_params(options.gff_params)
        coding_gff_file = coding_gff_dic[options.inspecies]
        genome_dic = utils.read_params(options.genome_params)
        genome_file = genome_dic[options.inspecies]
        ortho_dic = utils.ncar_ortho_dic(options.ncar_ortho_file, 4)
        coding_ortho_dic = utils.read_orthofile("orthofinder", options.ortho_file) 
        exclude_paras = True
        ncar_list = utils.ncar_min_taxa_membership({(options.inspecies, options.outspecies) : 2}, {}, [], "%s/%s_ncar_ortho.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        good_ncar_ortho_dic = {}
        for ncar in ncar_list:
            good_ncar_ortho_dic[ncar] = ortho_dic[ncar]
#        print good_ncar_ortho_dic["19393"]
        if options.outspecies == "DNOV":
            coding_outspecies = "Dnov"
        else:
            coding_outspecies = options.outspecies
        og_list = utils.min_taxa_membership({(options.inspecies, coding_outspecies) : 2}, {}, [], "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        good_coding_ortho_dic = {}
        for og in og_list:
            good_coding_ortho_dic[og] = coding_ortho_dic[og]
        pickle_dir = options.pickle_dir
        utils.ncar_mk_test(options.inspecies, options.outspecies, good_ncar_ortho_dic, "%s/%s_fsa_ncar" % (options.base_dir, options.prefix), "%s/%s_ncar_ancestral" % (options.base_dir, options.prefix), options.base_dir, options.num_threads, options.min_taxa, ncar_gff_file, genome_file, coding_gff_file, pickle_dir, good_coding_ortho_dic, "%s/%s_gene_ancestral" % (options.base_dir, options.prefix), "%s/%s_fsa_coding" % (options.base_dir, options.prefix))
        sys.exit()

    if options.action == "coding_mk":
        coding_gff_dic = utils.read_params(options.gff_params)
        coding_gff_file = coding_gff_dic[options.inspecies]
        genome_dic = utils.read_params(options.genome_params)
        genome_file = genome_dic[options.inspecies]
        coding_ortho_dic = utils.read_orthofile("orthofinder", options.ortho_file) 
        exclude_paras = True
        og_list = utils.min_taxa_membership({(options.inspecies, options.outspecies) : 2}, {}, [], "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        good_coding_ortho_dic = {}
        for og in og_list:
            good_coding_ortho_dic[og] = coding_ortho_dic[og]
        pickle_dir = options.pickle_dir
        utils.mk_test(options.inspecies, options.outspecies, good_coding_ortho_dic, "%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_gene_ancestral" % (options.base_dir, options.prefix), options.base_dir, options.num_threads, options.min_taxa, genome_file, coding_gff_file, pickle_dir)
        sys.exit()

    if options.action == "fixed_v_shared":
        coding_ortho_dic = utils.read_orthofile("orthofinder", options.ortho_file) 
        exclude_paras = True
        og_list = utils.min_taxa_membership({(options.inspecies, options.outspecies) : 2}, {}, [], "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        good_coding_ortho_dic = {}
#        og_list = [10644] #11419, 12394, 11141, 11231, 11334, 11341]
        for og in og_list:
            good_coding_ortho_dic[og] = coding_ortho_dic[og]

        pickle_dir = options.pickle_dir
        utils.gather_fixed_v_shared(options.inspecies, options.outspecies, good_coding_ortho_dic, "%s/%s_fsa_coding" % (options.base_dir, options.prefix), options.base_dir, options.num_threads, options.min_taxa, pickle_dir)
        sys.exit()

    if options.action == "fixed_v_shared":
        coding_ortho_dic = utils.read_orthofile("orthofinder", options.ortho_file) 
        exclude_paras = True
        og_list = utils.min_taxa_membership({(options.inspecies, options.outspecies) : 2}, {}, [], "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        good_coding_ortho_dic = {}
#        og_list = [10644] #11419, 12394, 11141, 11231, 11334, 11341]
        for og in og_list:
            good_coding_ortho_dic[og] = coding_ortho_dic[og]

        pickle_dir = options.pickle_dir
        utils.gather_fixed_v_shared(options.inspecies, options.outspecies, good_coding_ortho_dic, "%s/%s_fsa_coding" % (options.base_dir, options.prefix), options.base_dir, options.num_threads, options.min_taxa, pickle_dir)
        sys.exit()

    if options.action == "pairs_coding_div":
        coding_ortho_dic = utils.read_orthofile("orthofinder", options.ortho_file) 
        exclude_paras = True
        og_list = utils.min_taxa_membership({(options.inspecies, options.outspecies) : 2}, {}, [], "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        good_coding_ortho_dic = {}
#        og_list = [10644, 11419, 12394, 11141, 11231, 11334, 11341]
        for og in og_list:
            good_coding_ortho_dic[og] = coding_ortho_dic[og]

        pickle_dir = options.pickle_dir
        utils.pairs_coding_div(options.inspecies, options.outspecies, good_coding_ortho_dic, "%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), options.base_dir, options.num_threads, options.min_taxa)
        sys.exit()


    if options.action == "alignment_filter":
        index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
        paras_allowed = True
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        utils.alignment_column_filtering("%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_fsa_coding_columnfilt" % (options.base_dir, options.prefix), og_list, options.nogap_min_count, options.nogap_min_prop, options.nogap_min_species, {}, options.num_threads)
        print "First iteration of column filtering done. Results written to %s/%s_fsa_coding_columnfilt" % (options.base_dir, options.prefix)
        print "Starting Jarvis filter."
        utils.jarvis_filtering(og_list, "%s/%s_fsa_coding_columnfilt" % (options.base_dir, options.prefix), "%s/%s_fsa_coding_jarvis" % (options.base_dir, options.prefix), options.min_cds_len, options.num_threads)
        print "Jarvis filtering done. Results written to %s/%s_fsa_coding_jarvis" % (options.base_dir, options.prefix)
        utils.alignment_column_filtering("%s/%s_fsa_coding_jarvis" % (options.base_dir, options.prefix), "%s/%s_fsa_coding_jarvis_columnfilt" % (options.base_dir, options.prefix), og_list, options.nogap_min_count, options.nogap_min_prop, options.nogap_min_species, {}, options.num_threads)
        print "Second iteration of column filtering done. Results written to %s/%s_fsa_coding_jarvis_columnfilt" % (options.base_dir, options.prefix)
        utils.sequence_gap_filtering("%s/%s_fsa_coding_jarvis_columnfilt" % (options.base_dir, options.prefix), "%s/%s_fsa_coding_jarvis_columnfilt_seqfilt" % (options.base_dir, options.prefix), "%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), "%s/%s_orthos" % (options.base_dir, options.prefix), og_list, options.min_seq_prop_kept, options.max_seq_prop_gap, options.min_cds_len, "%s/%s_filtered.index" % (options.base_dir, options.prefix))
        print "Filtering of whole sequences based on gap content done. Results written to %s/%s_fsa_coding_jarvis_columnfilt_seqfilt" % (options.base_dir, options.prefix)
        print "Exiting"
        sys.exit()

    if options.action == "rer_converge":
        test_type = "aaml_blengths"
        foreground = "aaml_blengths"
        exclude_paras = True
        manda_taxa, multi_taxa, remove_list = utils.make_taxa_dic(options.taxa_inclusion)
        og_list = utils.min_taxa_membership(manda_taxa, multi_taxa, remove_list, "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        print len(og_list)
        utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, options.outputfile.split(".")[0]), options.tree_file, options.num_threads, options.use_gblocks, options.min_taxa, remove_list)
        utils.read_aaml_phylos(og_list, "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, options.outputfile.split(".")[0]), "%s/aaml_compiled" % (options.base_dir), options.outputfile, options.min_taxa)
        sys.exit()

    if options.action == "ncar_rer_converge": #untested (5/18/19)
        exclude_paras = True
        foreground = "baseml_blengths"
        manda_taxa, multi_taxa, remove_list = utils.make_taxa_dic(options.taxa_inclusion)
        ncar_list = utils.ncar_min_taxa_membership(manda_taxa, multi_taxa, remove_list, "%s/%s_ncar_ortho.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)

        print len(ncar_list)
#        ncar_list = ["ce99895"]
        utils.baseml_blengths(ncar_list, "%s/%s_fsa_ncar" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, options.outputfile.split(".")[0]), options.tree_file, options.num_threads, remove_list)
        utils.read_baseml_phylos(ncar_list, "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, options.outputfile.split(".")[0]), "%s/baseml_compiled" % (options.base_dir), options.outputfile)
        sys.exit()



    if options.action == "acacia_rescue":
        index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
        test_type = "aaml_blengths"
        foreground = "aaml_blengths"
        
        og_list = utils.read_ortho_index(index_file, options.min_taxa, True)
#        utils.read_aaml_phylos(og_list, "/Genomics/kocherlab/berubin/acacias/ant_genomes/comparative/acacia_selection/acacias_aaml_blengths_aaml_blengths", "%s/aaml_compiled" % (options.base_dir), options.outputfile, options.min_taxa)
        utils.read_aaml_phylos(og_list, "/Genomics/kocherlab/berubin/ants/bees/bees12_selection/bees_aaml_blengths_aaml_blengths", "%s/aaml_compiled" % (options.base_dir), options.outputfile, options.min_taxa)
        sys.exit()

    if options.action == "nopara_gene_trees":
        constrained = False
        paras_allowed = True
        include_paras = False
        og_list = utils.read_ortho_index("%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, paras_allowed)
        cur_og_list = og_list
        utils.gene_trees(cur_og_list, "%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), "%s/%s_nopara_nucl_gene_trees" % (options.base_dir, options.prefix), constrained, options.tree_file, options.num_threads, "nucs")
        sys.exit()

    if options.action == "check_discordance":
        paras_allowed = True
        og_list = utils.read_ortho_index("%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, paras_allowed)
        utils.discordance(og_list, "%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), "%s/%s_nopara_nucl_gene_trees" % (options.base_dir, options.prefix), "%s/%s_discordance" % (options.base_dir, options.prefix), options.tree_file, options.num_threads)
        utils.read_discordance("%s/%s_discordance" % (options.base_dir, options.prefix), og_list, options.base_dir)
        sys.exit()

    if options.action == "rer_goatools":
        if not os.path.exists("%s/RER_goatools" % options.base_dir):
            os.mkdir("%s/RER_goatools" % options.base_dir)
        rerconverge_output = options.rerconverge_output
        short_outputname = rerconverge_output.split("/")[-1][0:-4]
        utils.rer_goatools(rerconverge_output, rerconverge_output, options.go_database, "%s/RER_goatools/rer_0.05_slower_go_%s" % (options.base_dir, short_outputname), 3, 0.05, "slow")
        utils.rer_goatools(rerconverge_output, rerconverge_output, options.go_database, "%s/RER_goatools/rer_0.05_faster_go_%s" % (options.base_dir, short_outputname), 3, 0.05, "fast")
        sys.exit()

    if options.action == "goatools":
        outbase = options.goa_forefile.split("/")[-1].rsplit(".", 1)[0]
        if not os.path.exists("%s/%s_goatools/" % (options.base_dir, options.prefix)):
            os.mkdir("%s/%s_goatools/" % (options.base_dir, options.prefix))
        outdir = "%s/%s_goatools/%s" % (options.base_dir, options.prefix, outbase)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        utils.og_list_goatools(options.goa_forefile, options.goa_backfile, options.go_database, outdir)
        sys.exit()
        
    if options.action == "hyphy_relax":                
        test_type = "RELAX"
        exclude_paras = True
        og_list = []
        manda_taxa, multi_taxa, remove_list = utils.make_taxa_dic(options.taxa_inclusion)

        if options.og_list_file:

            reader = open(options.og_list_file, 'rU')
            for line in reader:
                cur_og = int(line.strip())
                og_list.append(cur_og)
        else:
        
            og_list = utils.min_taxa_membership(manda_taxa, multi_taxa, remove_list, "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)

        print len(og_list)
        if options.foreground == "INTREE":
            fore_list = "INTREE"
#            cur_og_list = og_list
        elif options.foreground.startswith("DAUGHTERS"):
            fore_list = options.foreground.split(",")
        else:
            fore_list = options.foreground.split(",")
#        og_list = [10724, 11488, 12704, 13036, 13879, 15282]
        print options.foreground
        utils.paml_test(og_list, fore_list, test_type, "%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, options.foreground, test_type), options.tree_file, options.num_threads, options.use_gblocks, options.min_taxa, remove_list)
        utils.read_hyphy_relax(og_list, "%s/%s_%s_%s" % (options.base_dir, options.prefix, options.foreground, test_type), options.base_dir, options.foreground)
        sys.exit()

    if options.action == "hyphy_absrel":
        test_type = "aBSREL"
        exclude_paras = True
        manda_taxa, multi_taxa, remove_list = utils.make_taxa_dic(options.taxa_inclusion)
        og_list = []
        if options.og_list_file:

            reader = open(options.og_list_file, 'rU')
            for line in reader:
                cur_og = int(line.strip())
                og_list.append(cur_og)
        else:
        
            og_list = utils.min_taxa_membership(manda_taxa, multi_taxa, remove_list, "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)

        
        print len(og_list)
#        print og_list[0:10]
        og_list = utils.limit_list(og_list, options.min_og_group, options.max_og_group)
#        og_list = og_list[0:10]
        utils.paml_test(og_list, [], test_type, "%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, "all", test_type), options.tree_file, options.num_threads, options.use_gblocks, options.min_taxa, remove_list)
        utils.read_hyphy_absrel(og_list, "%s/%s_%s_%s" % (options.base_dir, options.prefix, "all", test_type), options.base_dir)
        sys.exit()

    if options.action == "hyphy_ncar":
        ncar_gff_dic = utils.read_params(options.ncar_gff_params)
        coding_gff_dic = utils.read_params(options.gff_params)
        ncar_gff_file = ncar_gff_dic[options.inspecies]
        coding_gff_file = coding_gff_dic[options.inspecies]
        ncar_gff = utils.basic_gff(ncar_gff_file)
        coding_gff = utils.basic_gff(coding_gff_file)
        coding_ortho_dic = utils.read_orthofile("orthofinder", options.ortho_file) 
        ncar_ortho_dic = utils.ncar_ortho_dic(options.ncar_ortho_file, options.min_taxa)
        inspecies_ncar_dic = utils.inspecies_ncar_translate(ncar_ortho_dic, options.inspecies, ncar_gff)
        
        test_type = "hyphy_ncar"
        exclude_paras = True
        manda_taxa, multi_taxa, remove_list = utils.make_taxa_dic(options.taxa_inclusion)
        
        ncar_list = utils.ncar_min_taxa_membership(manda_taxa, multi_taxa, remove_list, "%s/%s_ncar_ortho.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        ncar_list = utils.limit_list_ncars(ncar_list, options.min_og_group, options.max_og_group)
        print len(ncar_list)
#        ncar_list = ["ce128829"]
#        ncar_list = ncar_list[0:10]
#        print ncar_list
#        ncar_list = ['ce10308', 'ce128603', 'ce141120', 'ce142916', 'ce14937', 'ce15561', 'ce156349', 'ce173808', 'ce20265', 'ce203710', 'ce225300', 'ce225303', 'ce225306', 'ce225308', 'ce225309', 'ce226206', 'ce226239', 'ce232728', 'ce245727', 'ce252796', 'ce254313', 'ce26150', 'ce26152', 'ce26153', 'ce26154', 'ce26159', 'ce271998', 'ce300387', 'ce324440', 'ce327441', 'ce328975', 'ce340768', 'ce350279', 'ce365105', 'ce376507', 'ce384958', 'ce387164', 'ce401670', 'ce401671', 'ce401672', 'ce426751', 'ce431990', 'ce457543', 'ce489792', 'ce493831', 'ce519625', 'ce525591', 'ce58859', 'ce87743', 'ce91377']

        if len(ncar_list) == 0:
            sys.exit()

#        utils.hyphy_noncoding(ncar_list, "%s/%s_fsa_ncar" % (options.base_dir, options.prefix), "%s/%s_all_%s" % (options.base_dir, options.prefix, test_type), options.tree_file, options.num_threads, remove_list, "%s/%s_gene_ancestral" % (options.base_dir, options.prefix), inspecies_ncar_dic, coding_gff, coding_ortho_dic, options.outspecies)
        utils.read_hyphy_noncoding(ncar_list, "%s/%s_all_%s" % (options.base_dir, options.prefix, test_type), options.base_dir)

        sys.exit()
        

    if options.action == "hypergeom":
        utils.hypergeom_test(options.hyper_pop, options.hyper_pop_cond, options.hyper_targets, options.hyper_targets_back)
        sys.exit()

    if options.action == "rer_hypergeom":
        utils.rer_hypergeom_test(options.hyper_pop, options.hyper_pop_cond, options.rerconverge_output, 0, 0.05, "fast")
        utils.rer_hypergeom_test(options.hyper_pop, options.hyper_pop_cond, options.rerconverge_output, 0, 0.05, "slow")
        utils.rer_hypergeom_test(options.hyper_pop, options.hyper_pop_cond, options.rerconverge_output, 0, 0.01, "fast")
        utils.rer_hypergeom_test(options.hyper_pop, options.hyper_pop_cond, options.rerconverge_output, 0, 0.01, "slow")
        sys.exit()

    if options.action == "hka":
        utils.hka_test(options.inspecies, options.outspecies, options.noncoding_seq_type, ortho_dic, options.base_dir, options.num_threads,"%s/%s_fsa_coding" % (options.base_dir, options.prefix))    
        sys.exit()

    if options.action == "godatabase":
        gaf_file = "/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/AMEL.gaf"
        gaf_file = "/Genomics/kocherlab/berubin/annotation/trinotate/PGRA/PGRA.gaf"
        gaf_file = "/Genomics/kocherlab/berubin/annotation/trinotate/ACEP/ACEP.gaf"
        gaf_file = "/Genomics/kocherlab/berubin/annotation/hic/trinotate/LALB/LALB.gaf"
        utils.make_go_database(ortho_dic, ipr_taxa_list, "%s/%s" % (options.base_dir, options.prefix), gaf_file)
        sys.exit()

    if options.action == "yn_dnds":
#        paras_allowed = True
        exclude_paras = True
        manda_taxa, multi_taxa, remove_list = utils.make_taxa_dic(options.taxa_inclusion)
        og_list = utils.min_taxa_membership(manda_taxa, multi_taxa, remove_list, "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        print len(og_list)
#        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        utils.yn_estimates(og_list, "%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), "%s/%s_yn" % (options.base_dir, options.prefix), options.tree_file, options.min_taxa, options.use_gblocks, remove_list)
        sys.exit()

    if options.action == "gc_content":
        if not os.path.isdir("%s/%s_gc_content" % (options.base_dir, options.prefix)):
            os.mkdir("%s/%s_gc_content" % (options.base_dir, options.prefix))
        paras_allowed = True
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        utils.gc_content(og_list, "%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_gc_content" % (options.base_dir, options.prefix))
        sys.exit()

    if options.action == "free_ratios":
        test_type = "free"
        foreground = "free"
        get_dn_ds = True
        exclude_paras = True
        manda_taxa, multi_taxa, remove_list = utils.make_taxa_dic(options.taxa_inclusion)
        if options.og_list_file:
            reader = open(options.og_list_file, 'rU')
            for line in reader:
                cur_og = int(line.strip())
                og_list.append(cur_og)
        else:
            og_list = utils.min_taxa_membership(manda_taxa, multi_taxa, remove_list, "%s/%s_filtered.index" % (options.base_dir, options.prefix), options.min_taxa, exclude_paras)
        print len(og_list)
#        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)

        utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads, options.use_gblocks, options.min_taxa, remove_list)
        cur_og_list = og_list
        utils.read_frees("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s_results" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), get_dn_ds, options.tree_file, cur_og_list)
        sys.exit()

    if options.action == "dnds_matrix":
        test_type = "free"
        foreground = "free"
        cds_dic = utils.read_params(options.param_file)
        index_file = "%s/%s_ortho_noparas.index" % (options.base_dir, options.prefix)
        paras_allowed = False
        og_list = utils.read_ortho_index(index_file, len(cds_dic.keys()), paras_allowed) #Gets only those OGs that have a single sequence for every species in the study. This is for making a sequence matrix that can be used for phylogenetics.
        og_list = utils.aligned_og_completeness(og_list, "%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), len(cds_dic.keys()))
        og_list = og_list[0:10]
        utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding_jarvis_columnfilt_seqfilt_noparas" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, options.outputfile.split(".")[0]), options.tree_file, options.num_threads, options.use_gblocks, len(cds_dic.keys()), [])
#        utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads, options.use_gblocks)
#        utils.concatenate_fourf_for_raxml("%s/%s_gene_ancestral" % (options.base_dir, options.prefix), "%s/%s_fourfold.afa" % (options.base_dir, options.prefix), og_list, cds_dic.keys())
        sys.exit()


    if options.action == "termfinder":
#        for species in ["AFLO", "AMEL", "MQUA", "MQUA_apis"]:
        for species in ["BIMP", "BTER"]:
#            utils.og_list_termfinder("/Genomics/kocherlab/berubin/alignment/10bees/revisions/leaveoneout/%s_out_advanced_slow.txt" % (species), "/Genomics/kocherlab/berubin/alignment/10bees/revisions/leaveoneout/%s_out_advanced_back.txt" % (species), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_leaveoneout/%s_out_advanced_slow_go" % (options.base_dir, species), -9, -9)
#            utils.og_list_termfinder("/Genomics/kocherlab/berubin/alignment/10bees/revisions/leaveoneout/%s_out_advanced_fast.txt" % (species), "/Genomics/kocherlab/berubin/alignment/10bees/revisions/leaveoneout/%s_out_advanced_back.txt" % (species), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_leaveoneout/%s_out_advanced_fast_go" % (options.base_dir, species), -9, -9)
            utils.og_list_termfinder("/Genomics/kocherlab/berubin/alignment/10bees/revisions/leaveoneout/apis_%s_fast.txt" % (species), "/Genomics/kocherlab/berubin/alignment/10bees/revisions/leaveoneout/apis_%s_back.txt" % (species), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_leaveoneout/apis_%s_fast_go" % (options.base_dir, species), -9, -9)
            utils.og_list_termfinder("/Genomics/kocherlab/berubin/alignment/10bees/revisions/leaveoneout/apis_%s_slow.txt" % (species), "/Genomics/kocherlab/berubin/alignment/10bees/revisions/leaveoneout/apis_%s_back.txt" % (species), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_leaveoneout/apis_%s_slow_go" % (options.base_dir, species), -9, -9)
        sys.exit()
        for index in range(0,100):
#            utils.og_list_termfinder("/Genomics/kocherlab/berubin/alignment/10bees/revisions/random_advanced_rers_genes/random_fore_rers_%s_fast.txt" % (index), "/Genomics/kocherlab/berubin/alignment/10bees/revisions/random_advanced_rers_genes/random_fore_rers_%s_back.txt" % (index), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_advanced_random_termfinder/random_fore_rers_%s_fast_go" % (options.base_dir, index), -9, -9)
#            utils.og_list_termfinder("/Genomics/kocherlab/berubin/alignment/10bees/revisions/random_advanced_rers_genes/random_fore_rers_%s_slow.txt" % (index), "/Genomics/kocherlab/berubin/alignment/10bees/revisions/random_advanced_rers_genes/random_fore_rers_%s_back.txt" % (index), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_advanced_random_termfinder/random_fore_rers_%s_slow_go" % (options.base_dir, index), -9, -9)
            utils.og_list_termfinder("/Genomics/kocherlab/berubin/alignment/10bees/revisions/random_ncar_sets/random_ncar_set_%s_fore.txt" % (index), "/Genomics/kocherlab/berubin/alignment/10bees/revisions/random_ncar_sets/random_ncar_set_%s_back.txt" % (index), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_random_ncar_set_termfinder/random_ncar_set_%s_go" % (options.base_dir, index), -9, -9)
#        utils.og_list_termfinder("/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/sig_midgut_genes.txt", "/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/Background_genes.txt","/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/AMEL.gaf", ortho_dic, "/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/midgut_go", -9, 0.05)
#        utils.og_list_termfinder("/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/sig_fat_genes.txt", "/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/Background_genes.txt","/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/AMEL.gaf", ortho_dic, "/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/fat_go", -9, 0.05)
#        utils.og_list_termfinder("/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/BrainDEG.txt", "/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/Background_genes.txt","/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/AMEL.gaf", ortho_dic, "/Genomics/kocherlab/berubin/annotation/trinotate/AMEL/rittschof/brain_go", -9, 0.05)
        sys.exit()
#        utils.og_list_termfinder("%s/lasihali_augochlorine_overlap_branch_sigs.txt" % (options.base_dir), "%s/lasihali_augochlorine_overlap_branch_tests.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/lasihali_augochlorine_overlap_branch_sigs_go" % (options.base_dir), -9, 0.05)
#        utils.og_list_termfinder("%s/lasihali_augochlorine_overlap_branch_slower_sigs.txt" % (options.base_dir), "%s/lasihali_augochlorine_overlap_branch_tests.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/lasihali_augochlorine_overlap_branch_sigs_slower_go" % (options.base_dir), -9, 0.05)
#        utils.og_list_termfinder("%s/lasihali_augochlorine_overlap_branch_faster_sigs.txt" % (options.base_dir), "%s/lasihali_augochlorine_overlap_branch_tests.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/lasihali_augochlorine_overlap_branch_sigs_faster_go" % (options.base_dir), -9, 0.05)
#        utils.og_list_termfinder("%s/gill_GNIG_branch.lrt_faster.lrt" % (options.base_dir), "%s/gill_GNIG_branch.lrt.lrt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/gill_GNIG_branch_go_faster" % (options.base_dir), -9, 0.05)
#        utils.og_list_termfinder("%s/gill_GNIG_branch.lrt_slower.lrt" % (options.base_dir), "%s/gill_GNIG_branch.lrt.lrt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/gill_GNIG_branch_go_slower" % (options.base_dir), -9, 0.05)
#        utils.og_list_termfinder("%s/gill_GNIG_bs.lrt.lrt" % (options.base_dir), "%s/gill_GNIG_bs.lrt.lrt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/gill_GNIG_bs_go_0.05" % (options.base_dir), -9, 0.05)
#        utils.og_list_termfinder("%s/gill_GNIG_bs.lrt.lrt" % (options.base_dir), "%s/gill_GNIG_bs.lrt.lrt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/gill_GNIG_bs_go_0.01" % (options.base_dir), -9, 0.01)
#        utils.og_list_termfinder("%s/bart_BNIG_branch.lrt_faster.lrt" % (options.base_dir), "%s/bart_BNIG_branch.lrt.lrt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/bart_BNIG_branch_go_faster" % (options.base_dir), -9, 0.05)
#        utils.og_list_termfinder("%s/bart_BNIG_branch.lrt_slower.lrt" % (options.base_dir), "%s/bart_BNIG_branch.lrt.lrt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/bart_BNIG_branch_go_slower" % (options.base_dir), -9, 0.05)
#        utils.og_list_termfinder("%s/bart_BNIG_bs.lrt.lrt" % (options.base_dir), "%s/bart_BNIG_bs.lrt.lrt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/bart_BNIG_bs_go_0.05" % (options.base_dir), -9, 0.05)
#        utils.og_list_termfinder("%s/bart_BNIG_bs.lrt.lrt" % (options.base_dir), "%s/bart_BNIG_bs.lrt.lrt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/bart_BNIG_bs_go_0.01" % (options.base_dir), -9, 0.01)
#        utils.og_list_termfinder("%s/../thoe_missing.txt" % (options.base_dir), "%s/../bnig_present.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/thoe_missing_go" % (options.base_dir), -9, -9)
#        utils.og_list_termfinder("%s/12bees_min10_gb.txt" % (options.base_dir), "%s/12bees_min10_gb.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/rer_0.01_slower_go" % (options.base_dir), 3, 0.01)
#        utils.og_list_termfinder("%s/12bees_min10_gb_cut0.001_bumble_advanced.txt" % (options.base_dir), "%s/12bees_min10_gb_cut0.001_bumble_advanced.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/rer_0.05_bumble_advanced_faster_go" % (options.base_dir), 3, 0.05)
#        utils.og_list_termfinder("%s/acacias_gb_cut0.001.txt_noNAs" % (options.base_dir), "%s/acacias_gb_cut0.001.txt_noNAs" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/rer_0.05_mutualism_slower_go_noNAs" % (options.base_dir), 3, 0.05)
#        utils.og_list_termfinder("%s/leafs_medium_gb_cut0.001.txt_noNAs" % (options.base_dir), "%s/leafs_medium_gb_cut0.001.txt_noNAs" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/rer_0.05_medium_slower_go_noNAs" % (options.base_dir), 3, 0.05)
#        utils.og_list_termfinder("%s/leafs_medium_gb_cut0.001.txt_noNAs" % (options.base_dir), "%s/leafs_medium_gb_cut0.001.txt_noNAs" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/rer_0.05_medium_faster_go_noNAs" % (options.base_dir), 3, 0.05)
        test_type = "mk"
#        for foreground in ["AAUR", "APUR", "HLIG", "LCAL", "LFIG", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP"]:
#            utils.og_list_termfinder("%s/mk_tests_noanc/%s_pos_conservative.lst" % (options.base_dir, foreground), "%s/mk_tests_noanc/%s_tested_conservative.lst" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/mk_tests_noanc/%s_%s_%s_conservative_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
#        for foreground in ["soc", "sol"]:
#            utils.og_list_termfinder("%s/mk_tests_noanc/og_list_%s.txt" % (options.base_dir, foreground), "%s/mk_tests_noanc/soc_sol_background.lst" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/mk_tests_noanc/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
        test_type = "noncoding"
#        utils.og_list_termfinder("%s/noncoding_outlier_tests/nonc_time_calibrated_nooverlaps_besthit_fastest100ogs.txt" % (options.base_dir), "%s/noncoding_outlier_tests/nonc_time_calibrated_nooverlaps_besthit_fastest100ogs_unsig.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_outlier_tests/nonc_time_calibrated_nooverlaps_besthit_fastest100ogs_go" % (options.base_dir), -9, -9)
#        utils.og_list_termfinder("%s/noncoding_outlier_tests/nonc_time_calibrated_nooverlaps_besthit_slowest100ogs.txt" % (options.base_dir), "%s/noncoding_outlier_tests/nonc_time_calibrated_nooverlaps_besthit_slowest100ogs_unsig.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_outlier_tests/nonc_time_calibrated_nooverlaps_besthit_slowest100ogs_go" % (options.base_dir), -9, -9)
#        utils.og_list_termfinder("%s/coding_outlier_tests/cod_nucs_calibrated_nooverlaps_fastest100.txt" % (options.base_dir), "%s/coding_outlier_tests/cod_nucs_calibrated_nooverlaps_fastest100_unsig.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/coding_outlier_tests/cod_nucs_calibrated_nooverlaps_fastest100_go" % (options.base_dir), -9, -9)
#        utils.og_list_termfinder("%s/coding_outlier_tests/cod_nucs_calibrated_nooverlaps_slowest100.txt" % (options.base_dir), "%s/coding_outlier_tests/cod_nucs_calibrated_nooverlaps_slowest100_unsig.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/coding_outlier_tests/cod_nucs_calibrated_nooverlaps_slowest100_go" % (options.base_dir), -9, -9)
#        utils.og_list_termfinder("%s/coding_outlier_tests/lowcod_highnonc.txt" % (options.base_dir), "%s/coding_outlier_tests/coding_and_noncoding.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/coding_outlier_tests/lowcod_highnonc_go" % (options.base_dir), -9, -9)
#        utils.og_list_termfinder("%s/coding_outlier_tests/highcod_lownonc.txt" % (options.base_dir), "%s/coding_outlier_tests/coding_and_noncoding.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/coding_outlier_tests/highcod_lownonc_go" % (options.base_dir), -9, -9)
        # utils.og_list_termfinder("%s/noncoding_cladespecific_tests/obligate_specific_ogs.txt" % (options.base_dir), "%s/noncoding_cladespecific_tests/nonc_time_calibrated_nooverlaps_besthit_all_ogs.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_cladespecific_tests/obligate_v_all_go" % (options.base_dir), -9, -9)
        # utils.og_list_termfinder("%s/noncoding_cladespecific_tests/corbiculate_specific_ogs.txt" % (options.base_dir), "%s/noncoding_cladespecific_tests/nonc_time_calibrated_nooverlaps_besthit_all_ogs.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_cladespecific_tests/corbiculate_v_all_go" % (options.base_dir), -9, -9)
        # utils.og_list_termfinder("%s/noncoding_cladespecific_tests/habplus_specific_ogs.txt" % (options.base_dir), "%s/noncoding_cladespecific_tests/nonc_time_calibrated_nooverlaps_besthit_all_ogs.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_cladespecific_tests/habplus_v_all_go" % (options.base_dir), -9, -9)
        # utils.og_list_termfinder("%s/noncoding_cladespecific_tests/megplus_specific_ogs.txt" % (options.base_dir), "%s/noncoding_cladespecific_tests/nonc_time_calibrated_nooverlaps_besthit_all_ogs.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_cladespecific_tests/megplus_v_all_go" % (options.base_dir), -9, -9)
        # utils.og_list_termfinder("%s/noncoding_cladespecific_tests/complex_specific_ogs.txt" % (options.base_dir), "%s/noncoding_cladespecific_tests/nonc_time_calibrated_nooverlaps_besthit_all_ogs.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_cladespecific_tests/complex_v_all_go" % (options.base_dir), -9, -9)
        # utils.og_list_termfinder("%s/noncoding_cladespecific_tests/carpenter_specific_ogs.txt" % (options.base_dir), "%s/noncoding_cladespecific_tests/nonc_time_calibrated_nooverlaps_besthit_all_ogs.txt" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/noncoding_cladespecific_tests/carpenter_v_all_go" % (options.base_dir), -9, -9)

        sys.exit()
        for foreground in ["all", "downstream", "upstream", "intron", "promoter"]:
#            utils.og_list_termfinder("%s/RER_noncoding/solloss_faster_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding/solloss_faster_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding/solloss_faster_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
#            utils.og_list_termfinder("%s/RER_noncoding/solloss_slower_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding/solloss_slower_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding/solloss_slower_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
#            utils.og_list_termfinder("%s/RER_noncoding/socgain_faster_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding/socgain_faster_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding/socgain_faster_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
#            utils.og_list_termfinder("%s/RER_noncoding/socgain_slower_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding/socgain_slower_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding/socgain_slower_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
#            utils.og_list_termfinder("%s/RER_noncoding/socgain_faster_solloss_slower_%s_overlap.txt" % (options.base_dir, foreground), "%s/RER_noncoding/socgain_faster_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding/socgain_faster_solloss_slower_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
#            utils.og_list_termfinder("%s/RER_noncoding/socgain_slower_solloss_faster_%s_overlap.txt" % (options.base_dir, foreground), "%s/RER_noncoding/socgain_slower_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding/socgain_slower_solloss_faster_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)

            utils.og_list_termfinder("%s/RER_noncoding_besthit/advanced_p0.05_faster_bests_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding_besthit/advanced_p0.05_faster_bests_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding_besthit/advanced_p0.05_faster_bests_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
            utils.og_list_termfinder("%s/RER_noncoding_besthit/advanced_p0.05_slower_bests_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding_besthit/advanced_p0.05_slower_bests_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding_besthit/advanced_p0.05_slower_bests_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
            utils.og_list_termfinder("%s/RER_noncoding_besthit/bumble_advanced_p0.05_faster_bests_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding_besthit/bumble_advanced_p0.05_faster_bests_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding_besthit/bumble_advanced_p0.05_faster_bests_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
            utils.og_list_termfinder("%s/RER_noncoding_besthit/bumble_advanced_p0.05_slower_bests_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding_besthit/bumble_advanced_p0.05_slower_bests_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding_besthit/bumble_advanced_p0.05_slower_bests_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
            utils.og_list_termfinder("%s/RER_noncoding_besthit/primitive_bumble_advanced_p0.05_faster_bests_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding_besthit/primitive_bumble_advanced_p0.05_faster_bests_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding_besthit/primitive_bumble_advanced_p0.05_faster_bests_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)
            utils.og_list_termfinder("%s/RER_noncoding_besthit/primitive_bumble_advanced_p0.05_slower_bests_sig_%s.txt" % (options.base_dir, foreground), "%s/RER_noncoding_besthit/primitive_bumble_advanced_p0.05_slower_bests_unsig_%s.txt" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_noncoding_besthit/primitive_bumble_advanced_p0.05_slower_bests_%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), -9, -9)

        sys.exit()


    if options.action == "rer_termfinder":
        if not os.path.exists("%s/RER_termfinder" % options.base_dir):
            os.mkdir("%s/RER_termfinder" % options.base_dir)
        rerconverge_output = options.rerconverge_output
        short_outputname = rerconverge_output.split("/")[-1][0:-4]
        utils.rer_goatools(rerconverge_output, rerconverge_output, options.go_database, "%s/RER_termfinder/rer_0.05_slower_go_%s" % (options.base_dir, short_outputname), 3, 0.05, "slow")
        utils.rer_goatools(rerconverge_output, rerconverge_output, options.go_database, "%s/RER_termfinder/rer_0.05_faster_go_%s" % (options.base_dir, short_outputname), 3, 0.05, "fast")

        # utils.rer_termfinder(rerconverge_output, rerconverge_output, "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_termfinder/rer_0.05_slower_go_%s" % (options.base_dir, short_outputname), 3, 0.05, "slow")
        # utils.rer_termfinder(rerconverge_output, rerconverge_output, "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_termfinder/rer_0.05_faster_go_%s" % (options.base_dir, short_outputname), 3, 0.05, "fast")
        # utils.rer_termfinder(rerconverge_output, rerconverge_output, "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_termfinder/rer_0.01_slower_go_%s" % (options.base_dir, short_outputname), 3, 0.01, "slow")
        # utils.rer_termfinder(rerconverge_output, rerconverge_output, "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_termfinder/rer_0.01_faster_go_%s" % (options.base_dir, short_outputname), 3, 0.01, "fast")
        # utils.og_list_termfinder(rerconverge_output, rerconverge_output, "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_termfinder/rer_0.01_all_go_%s" % (options.base_dir, short_outputname), 3, 0.01)
        # utils.og_list_termfinder(rerconverge_output, rerconverge_output, "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_termfinder/rer_0.05_all_go_%s" % (options.base_dir, short_outputname), 3, 0.05)
        sys.exit()
    if options.action == "rer_random_termfinder":
        if not os.path.exists("%s/RER_random_termfinder" % options.base_dir):
            os.mkdir("%s/RER_random_termfinder" % options.base_dir)
        rerconverge_output = options.rerconverge_output
        short_outputname = rerconverge_output.split("/")[-1][0:-4]
        utils.rer_termfinder(rerconverge_output, rerconverge_output, "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_random_termfinder/rer_0.05_slower_go_%s" % (options.base_dir, short_outputname), 3, 0.05, "slow")
        utils.rer_termfinder(rerconverge_output, rerconverge_output, "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/RER_random_termfinder/rer_0.05_faster_go_%s" % (options.base_dir, short_outputname), 3, 0.05, "fast")
        sys.exit()

    if options.action == "time_aamls":
        paras_allowed = True
        foreground = "aaml_blengths"
        test_type = "aaml_blengths"
        fore_list = options.foreground.split(",")
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
#        og_list = [6627]
        utils.aaml_time_phylos(og_list, "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_aaml_time_calibrated" % (options.base_dir, options.prefix), options.timetree, fore_list)
        sys.exit()
        
    if options.action == "ds_correlations":
        paras_allowed = True
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
#        og_list = og_list[0:100]
        bootstrap_taxa = True
        categorical = False
        utils.bootstrapping_ds_time_correlations(og_list, "%s/%s_free_free" % (options.base_dir, options.prefix), "%s/%s_ds_corrs" % (options.base_dir, options.prefix), options.timetree, options.traittree, bootstrap_taxa, categorical)
        sys.exit()

    if options.action == "rer_autism":
        if not os.path.exists("%s/RER_autism" % options.base_dir):
            os.mkdir("%s/RER_autism" % options.base_dir)
        rerconverge_output = options.rerconverge_output
        short_outputname = rerconverge_output.split("/")[-1][0:-4]
        autism_dir = "/Genomics/kocherlab/berubin/annotation/hic/autism/LALB"
        reference_species = "LALB"
        utils.rer_hypergeom("%s/LALB_human_bestreciprocal.txt" % autism_dir, "%s/LALB_sfaris.txt" % autism_dir, rerconverge_output, ortho_dic, -9, 0.05, "fast", "%s/RER_autism/rer_0.05_faster_%s_allsfari.txt" % (options.base_dir, short_outputname))
        utils.rer_hypergeom("%s/LALB_human_bestreciprocal.txt" % autism_dir, "%s/LALB_sfaris.txt" % autism_dir, rerconverge_output, ortho_dic, -9, 0.05, "slow", "%s/RER_autism/rer_0.05_slower_%s_allsfari.txt" % (options.base_dir, short_outputname))
        utils.rer_hypergeom("%s/LALB_human_bestreciprocal.txt" % autism_dir, "%s/LALB_sfaris.txt" % autism_dir, rerconverge_output, ortho_dic, -9, 0.05, "all", "%s/RER_autism/rer_0.05_all_%s_allsfari.txt" % (options.base_dir, short_outputname))
        utils.rer_hypergeom("%s/LALB_human_bestreciprocal.txt" % autism_dir, "%s/LALB_neuros.txt" % autism_dir, rerconverge_output, ortho_dic, -9, 0.05, "fast", "%s/RER_autism/rer_0.05_faster_%s_neuro.txt" % (options.base_dir, short_outputname))
        utils.rer_hypergeom("%s/LALB_human_bestreciprocal.txt" % autism_dir, "%s/LALB_neuros.txt" % autism_dir, rerconverge_output, ortho_dic, -9, 0.05, "slow", "%s/RER_autism/rer_0.05_slower_%s_neuro.txt" % (options.base_dir, short_outputname))
        utils.rer_hypergeom("%s/LALB_human_bestreciprocal.txt" % autism_dir, "%s/LALB_neuros.txt" % autism_dir, rerconverge_output, ortho_dic, -9, 0.05, "all", "%s/RER_autism/rer_0.05_all_%s_neuro.txt" % (options.base_dir, short_outputname))
        for syndrome in ["1", "2", "3", "4"]:
            utils.rer_hypergeom("%s/LALB_human_schizo_bestreciprocal.txt" % autism_dir, "%s/LALB_schizos_%s.txt" % (autism_dir, syndrome), rerconverge_output, ortho_dic, -9, 0.05, "fast", "%s/RER_autism/rer_0.05_faster_%s_%sschizo.txt" % (options.base_dir, short_outputname, syndrome))
            utils.rer_hypergeom("%s/LALB_human_schizo_bestreciprocal.txt" % autism_dir, "%s/LALB_schizos_%s.txt" % (autism_dir, syndrome), rerconverge_output, ortho_dic, -9, 0.05, "slow", "%s/RER_autism/rer_0.05_slower_%s_%sschizo.txt" % (options.base_dir, short_outputname, syndrome))
            utils.rer_hypergeom("%s/LALB_human_schizo_bestreciprocal.txt" % autism_dir, "%s/LALB_schizos_%s.txt" % (autism_dir, syndrome), rerconverge_output, ortho_dic, -9, 0.05, "all", "%s/RER_autism/rer_0.05_all_%s_%sschizo.txt" % (options.base_dir, short_outputname, syndrome))

        for syndrome in ["0", "1", "2", "3", "4", "5", "6", "S", "S1", "S12", "S123", "S1234", "S12345", "S123456"]:
            utils.rer_hypergeom("%s/LALB_human_bestreciprocal.txt" % autism_dir, "%s/LALB_sfaris_%s.txt" % (autism_dir, syndrome), rerconverge_output, ortho_dic, -9, 0.05, "fast", "%s/RER_autism/rer_0.05_faster_%s_%ssfari.txt" % (options.base_dir, short_outputname, syndrome))
            utils.rer_hypergeom("%s/LALB_human_bestreciprocal.txt" % autism_dir, "%s/LALB_sfaris_%s.txt" % (autism_dir, syndrome), rerconverge_output, ortho_dic, -9, 0.05, "slow", "%s/RER_autism/rer_0.05_slower_%s_%ssfari.txt" % (options.base_dir, short_outputname, syndrome))
            utils.rer_hypergeom("%s/LALB_human_bestreciprocal.txt" % autism_dir, "%s/LALB_sfaris_%s.txt" % (autism_dir, syndrome), rerconverge_output, ortho_dic, -9, 0.05, "all", "%s/RER_autism/rer_0.05_all_%s_%ssfari.txt" % (options.base_dir, short_outputname, syndrome))
        sys.exit()
    if options.action == "rer_stratigraphy":
        if not os.path.exists("%s/RER_stratigraphy" % options.base_dir):
            os.mkdir("%s/RER_stratigraphy" % options.base_dir)
        rerconverge_output = options.rerconverge_output
        short_outputname = rerconverge_output.split("/")[-1][0:-4]
        reference_species = "LALB"
        for syndrome in ["3", "4", "14", "15", "16", "17", "18", "19", "21"]:
            utils.rer_hypergeom("/Genomics/kocherlab/berubin/stratigraphy/LALB_all_levels.txt", "/Genomics/kocherlab/berubin/stratigraphy/LALB_level_%s.txt" % syndrome, rerconverge_output, ortho_dic, -9, 0.05, "fast", "%s/RER_stratigraphy/rer_0.05_faster_%s_level_%s.txt" % (options.base_dir, short_outputname, syndrome))
            utils.rer_hypergeom("/Genomics/kocherlab/berubin/stratigraphy/LALB_all_levels.txt", "/Genomics/kocherlab/berubin/stratigraphy/LALB_level_%s.txt" % syndrome, rerconverge_output, ortho_dic, -9, 0.05, "slow", "%s/RER_stratigraphy/rer_0.05_slower_%s_level_%s.txt" % (options.base_dir, short_outputname, syndrome))
            utils.rer_hypergeom("/Genomics/kocherlab/berubin/stratigraphy/LALB_all_levels.txt", "/Genomics/kocherlab/berubin/stratigraphy/LALB_level_%s.txt" % syndrome, rerconverge_output, ortho_dic, -9, 0.05, "all", "%s/RER_stratigraphy/rer_0.05_all_%s_level_%s.txt" % (options.base_dir, short_outputname, syndrome))
        sys.exit()

    if options.action == "autism_enrich":
#        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_free_free_results/soc_larger.txt" % (options.base_dir, options.prefix), ortho_dic)
#        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_lasihali_bs.lrt" % (options.base_dir, options.prefix), ortho_dic)
        print "\nAbdomen sigs"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/expression/abd_sigs.txt" % (options.base_dir), ortho_dic, -9)
        print "\nAntennae sigs"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/expression/ant_sigs.txt" % (options.base_dir), ortho_dic, -9)
        print "\nHead sigs"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/expression/head_sigs.txt" % (options.base_dir), ortho_dic, -9)
        print "\nSolitary marine 250"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/marine_tests_sol.txt" % (options.base_dir, options.prefix), ortho_dic, 250)
        print "\nSolitary marine 100"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/marine_tests_sol.txt" % (options.base_dir, options.prefix), ortho_dic, 100)
        print "\nSocial marine 250"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/marine_tests_soc.txt" % (options.base_dir, options.prefix), ortho_dic, 250)
        print "\nSocial marine 100"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/marine_tests_soc.txt" % (options.base_dir, options.prefix), ortho_dic, 100)
        print "\nSlower social Clark significant"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/soc_slower_clark.txt" % (options.base_dir, options.prefix), ortho_dic, -1)
        print "\nFaster social Clark significant"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/soc_faster_clark.txt" % (options.base_dir, options.prefix), ortho_dic, -1)
        print "\nSlower solitary Clark significant"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/sol_slower_clark.txt" % (options.base_dir, options.prefix), ortho_dic, -1)
        print "\nFaster solitary Clark significant"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/sol_faster_clark.txt" % (options.base_dir, options.prefix), ortho_dic, -1)
        print "\nAll solitary Clark significant"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/sol_sig_clark.txt" % (options.base_dir, options.prefix), ortho_dic, -1)
        print "\nAll social Clark significant"
        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_aaml_blengths_aaml_blengths_results/soc_sig_clark.txt" % (options.base_dir, options.prefix), ortho_dic, -1)
        pairs_list = [("LMAR", "LFIG"), ("LZEP", "LVIE"), ("LPAU", "LOEN"), ("AAUR", "APUR")]
#        for pair in pairs_li
#        utils.og_list_termfinder("%s/%s_aaml_blengths_aaml_blengths_results/sol_sig_clark.txt" % (options.base_dir, options.prefix), "%s/%s_aaml_blengths_aaml_blengths_results/compiled_aaml_phylos.txt" % (options.base_dir, options.prefix),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_aaml_blengths_aaml_blengths_results/sol_sig_clark_go" % (options.base_dir, options.prefix))
#        utils.og_list_termfinder("%s/%s_aaml_blengths_aaml_blengths_results/soc_sig_clark.txt" % (options.base_dir, options.prefix), "%s/%s_aaml_blengths_aaml_blengths_results/compiled_aaml_phylos.txt" % (options.base_dir, options.prefix),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_aaml_blengths_aaml_blengths_results/soc_sig_clark_go" % (options.base_dir, options.prefix), -1, -1)
#        utils.og_list_termfinder("%s/expression/abd_sigs.txt" % (options.base_dir), "%s/expression/abd.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/expression/abd_go" % (options.base_dir), -9, 0.05)
        utils.og_list_termfinder("%s/expression/ant_sigs.txt" % (options.base_dir), "%s/expression/ant.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/expression/ant_go" % (options.base_dir), -9, 0.05)
        utils.og_list_termfinder("%s/expression/head_sigs.txt" % (options.base_dir), "%s/expression/head.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/expression/head_go" % (options.base_dir), -9, 0.05)
        utils.og_list_termfinder("%s/%s_aaml_blengths_aaml_blengths_results/marine_tests_sol.txt" % (options.base_dir, options.prefix), "%s/%s_aaml_blengths_aaml_blengths_results/compiled_aaml_phylos.txt" % (options.base_dir, options.prefix),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_aaml_blengths_aaml_blengths_results/marine_tests_sol_go" % (options.base_dir, options.prefix), 1, 0.05)
        utils.og_list_termfinder("%s/%s_aaml_blengths_aaml_blengths_results/marine_tests_soc.txt" % (options.base_dir, options.prefix), "%s/%s_aaml_blengths_aaml_blengths_results/compiled_aaml_phylos.txt" % (options.base_dir, options.prefix),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_aaml_blengths_aaml_blengths_results/marine_tests_soc_go" % (options.base_dir, options.prefix), 1, 0.05)
        pairs_list = [("LMAR", "LFIG"), ("LZEP", "LVIE"), ("LPAU", "LOEN"), ("AAUR", "APUR")]
        for pair in pairs_list:
            outlier_num = 100
            print pair
            utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/hka_tests/%s_%s_flank_hka_direct_pairwise_results_correct_slower.txt" % (options.base_dir, pair[0], pair[1]), ortho_dic, outlier_num)
            utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/hka_tests/%s_%s_flank_hka_direct_pairwise_results_correct_slower.txt" % (options.base_dir, pair[1], pair[0]), ortho_dic, outlier_num)
#        for terminal in ["LVIE", "LZEP", "LPAU", "LOEN", "LMAR", "LMAL", "LLEU", "LFIG", "HLIG", "AAUR", "APUR"]:
#            print terminal
#            utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/%s_terminals/%s_%s_branch_slower.lrt" % (options.base_dir, options.prefix, options.prefix, terminal), ortho_dic, -9)

#        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/hka_tests/LZEP_LVIE_flank_hka_direct_pairwise_results_correct_slower.txt" % (options.base_dir), ortho_dic, outlier_num)
#        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/hka_tests/LVIE_LZEP_flank_hka_direct_pairwise_results_correct_slower.txt" % (options.base_dir), ortho_dic, outlier_num)
#        utils.hypergeom("/Genomics/kocherlab/berubin/annotation/autism/LALB_human_bestreciprocal.txt", "/Genomics/kocherlab/berubin/annotation/autism/LALB_sfaris.txt", "%s/hka_tests/LPAU_LOEN_flank_hka_direct_pairwise_results_correct_slower.txt" % (options.base_dir), ortho_dic, outlier_num)
        sys.exit()
#        utils.og_list_termfinder("%s/mk_tests_noanc/og_list_soc.txt" % (options.base_dir), "%s/mk_tests_noanc/soc_sol_background.lst" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/mk_tests_noanc/%s_soc_%s_go" % (options.base_dir, options.prefix, test_type))

#    seq_dic = utils.get_cds() #get coding sequences from all species
#    seq_dic = utils.get_cds("/Genomics/kocherlab/berubin/ants/ant_data/renamed_cds_seqs")
    #store name of orthology index file
    index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
    #write fastas for ALL orthologous groups
    #utils.write_orthos(options.ortho_file, seq_dic, True, "%s/%s_orthos" % (options.base_dir, options.prefix), index_file)
    paras_allowed = False #do not include OG's with paralogs
    #get a list of the first 100 genes with min_taxa (should always be all
    #of the available species) for making the phylogeny
    og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)#[0:1000]
    iscoding = True
#    og_list = ["vha"]
#    og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
#    utils.fsa_coding_align(og_list, "%s/%s_orthos/" % (options.base_dir, options.prefix), "%s/%s_fsa_coding" % (options.base_dir, options.prefix), options.num_threads, iscoding)
    use_backbone = False #do not use a phylogeny to align these first 100 genes
    #align first 100 genes
    #utils.prank_align(og_list, "%s/%s_orthos/" % (options.base_dir, options.prefix), "%s/%s_prank_no_backbone" % (options.base_dir, options.prefix), use_backbone, "nophylogeny", options.num_threads)
    #concatenate the aligned 100 genes
#    utils.concatenate_for_raxml("%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s.afa" % (options.base_dir, options.prefix), og_list)

    #then run raxml to create a backbone phylogeny
    #raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s %s/%s.afa % (options.base_dir, options.prefix) -T 20
    #get list of all genes without paralogs
    
#    og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)

#    utils.yn_estimates(og_list, "%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_yn" % (options.base_dir, options.prefix))
#    og_list = utils.limit_list(og_list, options.min_og_group, options.max_og_group)

    use_backbone = True #use phylogeny to align these genes
    #align non-paralogous genes
#    utils.prank_align(og_list, "%s/%s_orthos/" % (options.base_dir, options.prefix), "%s/%s_prank" % (options.base_dir, options.prefix), use_backbone, options.tree_file, options.num_threads)
    test_type = "ancestral"
    foreground = "dummy"
#    og_list = ["vha"]
#    cur_og_list = ["vha"]
#    utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
    constrained = True
#    utils.gene_trees(og_list, "%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_gene_trees_constrained" % (options.base_dir, options.prefix), constrained, options.tree_file, options.num_threads)
#    utils.compile_gene_trees(og_list, "%s/%s_gene_trees_constrained" % (options.base_dir, options.prefix), options.tree_file, "%s/%s_gene_trees_bls" % (options.base_dir, options.prefix))
#    target_taxa = ["SPRA", "SPIE", "SGLO"] #["ECOL", "SFLE", "SPRA", "SPIE"]
#    cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
#    cur_og_list = og_list
#    utils.compile_ancestrals(cur_og_list, "%s/%s_%s" % (options.base_dir, options.prefix, test_type), "%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_%s" % (options.base_dir, options.prefix, "substitution_distances_free"))

#    og_list = utils.external_list("%s/%s_free_free_results/sol_larger.txt" % (options.base_dir, options.prefix))
#    utils.compile_ancestrals(og_list, "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_%s" % (options.base_dir, options.prefix, "substitution_distances_free_sol"))
#    og_list = utils.external_list("%s/%s_free_free_results/soc_larger.txt" % (options.base_dir, options.prefix))
#    utils.compile_ancestrals(og_list, "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_%s" % (options.base_dir, options.prefix, "substitution_distances_free_soc"))


    test_type = "mk"
    """
#    for species in ["AAUR", "APUR", "HLIG", "LCAL", "LFIG", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP"]:
#        foreground = species
#        utils.og_list_termfinder("%s/mk_tests_noanc/%s_pos_adjusted.lst" % (options.base_dir, foreground), "%s/mk_tests_noanc/%s_tested_conservative.lst" % (options.base_dir, foreground), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/mk_tests_noanc/%s_%s_%s_adjusted_go" % (options.base_dir, options.prefix, foreground, test_type))

    utils.og_list_termfinder("%s/mk_tests_noanc/og_list_soc.txt" % (options.base_dir), "%s/mk_tests_noanc/soc_sol_background.lst" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/mk_tests_noanc/%s_soc_%s_go" % (options.base_dir, options.prefix, test_type))
    utils.og_list_termfinder("%s/mk_tests_noanc/og_list_sol.txt" % (options.base_dir), "%s/mk_tests_noanc/soc_sol_background.lst" % (options.base_dir), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/mk_tests_noanc/%s_sol_%s_go" % (options.base_dir, options.prefix, test_type))
    """
    if options.action == "branch_test":
        test_type = "branch"
        foreground = options.foreground
#        foreground = "BNIG"
#        foreground = "GNIG"
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        cur_og_list = og_list
#        target_taxa = ["BNIG", "THOE"]
#        target_taxa = ["GNIG", "GAPI"]
#        cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
        utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
        utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s.lrt" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
        sys.exit()

    if options.action == "bs_test":
        test_type = "bs"
#        foreground = "BNIG"
#        foreground = "GNIG"
        foreground = "solitary"
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        cur_og_list = og_list
#        target_taxa = ["GNIG", "GAPI"]
#        target_taxa = ["BNIG", "THOE"]
#        cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
        foreground = options.foreground
        utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads, options.use_gblocks, options.min_taxa)
        utils.test_lrt("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s.lrt" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))

        sys.exit()

    """
    test_type = "branch"
    foreground = "solitary"
    solitary_taxa = ["LLEU", "LFIG", "LVIE", "LOEN", "APUR"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, solitary_taxa, og_list, 2)
    utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#    utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s.lrt" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
    
    test_type = "bs"
    foreground = "solitary"
    solitary_taxa = ["LLEU", "LFIG", "LVIE", "LOEN", "APUR"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, solitary_taxa, og_list, 2)
    utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads
)#    utils.test_lrt("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s.lrt" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
    test_type = "branch"
    foreground = "social"
    utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#    utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s.lrt" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
    test_type = "bs"
    foreground = "social"
    utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#    utils.test_lrt("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s.lrt" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
    """
    
    """
    test_type = "branch"
    foreground = "terminals"
    target_taxa = ["LLEU", "LFIG", "LVIE", "LOEN", "APUR", "AAUR", "LMAR", "LZEP", "LMAL", "LPAU", "HLIG"]
    solitary_taxa = ["LLEU", "LFIG", "LVIE", "LOEN", "APUR"]
    social_taxa = ["AAUR", "LMAR", "LZEP", "LMAL", "LPAU", "HLIG"]
    pairs_list = [("AAUR", "APUR"), ("LMAR", "LFIG"), ("LZEP", "LVIE"), ("LPAU", "LOEN")]
    
    """
    """
    cur_og_list = utils.min_taxa_membership(ortho_dic, solitary_taxa, og_list, 3)
    cur_og_list = utils.min_taxa_membership(ortho_dic, social_taxa, cur_og_list, 4)
    for foreground in target_taxa:
#        utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_terminals/%s_%s_%s" % (options.base_dir, options.prefix, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
        utils.test_lrt_branch("%s/%s_terminals/%s_%s_%s" % (options.base_dir, options.prefix, options.prefix, foreground, test_type), "%s/%s_terminals/%s_%s_%s" % (options.base_dir, options.prefix, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_terminals/%s_%s_%s_go" % (options.base_dir, options.prefix,options.prefix, foreground, test_type))
    
    utils.terminal_test_overlap("%s/%s_terminals/" % (options.base_dir, options.prefix), options.prefix, test_type, target_taxa, social_taxa, solitary_taxa, pairs_list)
    """
    """
    test_type = "bs"
    foreground = "terminals"
    target_taxa = ["LLEU", "LFIG", "LVIE", "LOEN", "APUR", "AAUR", "LMAR", "LZEP", "LMAL", "LPAU", "HLIG"]
    solitary_taxa = ["LLEU", "LFIG", "LVIE", "LOEN", "APUR"]
    social_taxa = ["AAUR", "LMAR", "LZEP", "LMAL", "LPAU", "HLIG"]
    pairs_list = [("AAUR", "APUR"), ("LMAR", "LFIG"), ("LZEP", "LVIE"), ("LPAU", "LOEN")]
    cur_og_list = utils.min_taxa_membership(ortho_dic, solitary_taxa, og_list, 3)
    cur_og_list = utils.min_taxa_membership(ortho_dic, social_taxa, cur_og_list, 4)
    for foreground in target_taxa:
        utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_terminals/%s_%s_%s" % (options.base_dir, options.prefix, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#        utils.test_lrt_branch("%s/%s_terminals/%s_%s_%s" % (options.base_dir, options.prefix, options.prefix, foreground, test_type), "%s/%s_terminals/%s_%s_%s" % (options.base_dir, options.prefix, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_terminals/%s_%s_%s_go" % (options.base_dir, options.prefix,options.prefix, foreground, test_type))
#    utils.terminal_test_overlap("%s/%s_terminals/" % (options.base_dir, options.prefix, options.prefix), test_type, target_taxa, social_taxa, solitary_taxa, pairs_list)
    """
    test_type = "bs"
    foreground = "lasiaugo_clade"
#    target_taxa = ["AVIR"]
#    og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)

#    outgroup_taxa = ["Dnov", "Nmel", "AAUR", "APUR", "AVIR"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, outgroup_taxa, og_list, 2)
#    outgroup_taxa = ["Dnov", "Nmel", "AVIR"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, outgroup_taxa, cur_og_list, 1)
#    print len(cur_og_list)
#    utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#    utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
    
    test_type = "bs"
    foreground = "lasihali_clade"
#    target_taxa = ["AVIR"]
#    og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)

#    outgroup_taxa = ["Dnov", "Nmel", "AAUR", "APUR", "AVIR"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, outgroup_taxa, og_list, 2)
#    outgroup_taxa = ["Dnov", "Nmel", "AVIR"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, outgroup_taxa, cur_og_list, 1)
#    print len(cur_og_list)
#    utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#    utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
    
    test_type = "bs"
    foreground = "augochlorine_clade"
#    target_taxa = ["APUR", "AAUR"]
#    cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
#    outgroup_taxa = ["Dnov", "Nmel"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, outgroup_taxa, cur_og_list, 1)
#    print "Number OGs: %s" % len(og_list)
#    utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#    utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
    
    test_type = "bs"
    foreground = "halictus_clade"
#    target_taxa = ["HLIG", "HRUB"]
#    cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
#    outgroup_taxa = ["Dnov", "Nmel", "AAUR", "APUR"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, outgroup_taxa, cur_og_list, 2)
#    utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#    utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))

    test_type = "bs"
    foreground = "lasioglossum_clade"
#    target_taxa = ["AVIR"]
#    og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
#    outgroup_taxa = ["Dnov", "Nmel", "AAUR", "APUR", "AVIR"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, outgroup_taxa, og_list, 2)
#    ingroup_taxa = ["LLEU", "LFIG", "LMAR", "LVIE", "LZEP", "LCAL", "LALB", "LMAL", "LOEN", "LPAU"]
#    cur_og_list = utils.min_taxa_membership(ortho_dic, ingroup_taxa, cur_og_list, 5)
#    utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#    utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
    
    
    test_type = "free"
    foreground = "free"
    get_dn_ds = True
#    utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)

#    utils.read_frees("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s_results" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), get_dn_ds, options.timetree)

    test_type = "aaml_blengths"
    foreground = "aaml_blengths"
#    utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
    phylo_dic = utils.read_aaml_phylos(og_list, "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s_results" % (options.base_dir, options.prefix, foreground, test_type))
    fore_list = ["APUR", "LFIG", "LOEN", "LVIE", "LLEU"]
#    fore_list = ["AAUR", "LMAR", "LPAU", "LZEP", "LMAL", "HLIG"]
    exclude_list = ["LCAL", "LALB", "HRUB", "NMEL", "DNOV", "AVIR"]
#    utils.marine_test(phylo_dic, options.tree_file, fore_list, exclude_list, "%s/%s_%s_%s_results" % (options.base_dir, options.prefix, foreground, test_type))
    target_taxa = ["APUR", "LFIG", "LOEN", "LVIE", "AAUR", "LMAR", "LPAU", "LZEP"]
    cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
    utils.codon_bias(options.base_dir, SPECIES_LIST, cur_og_list, ortho_dic)
    
    """
    paras_allowed = True #allow paralogs
    #get list of all OG's including those with paralogs
    og_list_with_paras = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
    og_list_with_paras = utils.limit_list(og_list_with_paras, options.min_og_group, options.max_og_group)
    #this should be just the OG's with paralogs in them
    para_list = list(set(og_list_with_paras) - set(og_list))
    use_backbone = False # turn backbone phylogeny off
    #align all OG's with paralogs -- can't use backbone phylogeny for this
    utils.prank_align(para_list, "%s/%s_orthos/" % (options.base_dir, options.prefix), "%s/%s_prank" % (options.base_dir, options.prefix), use_backbone, options.tree_file, options.num_threads)
    paml_test(og_list, foreground, test_type,"%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), phylogeny_file, options.num_threads)
    foreground = "social"
    test_type = "model_d"
    test_type = "bs"
#    test_type = "free"
#    paml_test(og_list, foreground, test_type,"%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), phylogeny_file, options.num_threads)
#    read_frees("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type))
#    test_lrt("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type))
"""

if __name__ == '__main__':
    main()
