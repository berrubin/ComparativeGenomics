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
        print og_list
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

        if len(ncar_list) == 0:
            sys.exit()

        utils.hyphy_noncoding(ncar_list, "%s/%s_fsa_ncar" % (options.base_dir, options.prefix), "%s/%s_all_%s" % (options.base_dir, options.prefix, test_type), options.tree_file, options.num_threads, remove_list, "%s/%s_gene_ancestral" % (options.base_dir, options.prefix), inspecies_ncar_dic, coding_gff, coding_ortho_dic, options.outspecies)
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
        index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
        fore_list = options.foreground.split(",")
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
#        og_list = [3576]
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

    if options.action == "branch_test":
        test_type = "branch"
        foreground = options.foreground
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        cur_og_list = og_list
#        cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
        utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
        utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s.lrt" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
        sys.exit()

    if options.action == "bs_test":
        test_type = "bs"
        foreground = "solitary"
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        cur_og_list = og_list
        foreground = options.foreground
        utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads, options.use_gblocks, options.min_taxa)
        utils.test_lrt("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s.lrt" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))

        sys.exit()


if __name__ == '__main__':
    main()
