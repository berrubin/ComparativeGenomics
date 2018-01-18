#!/usr/local/bin/python
import os
import paml_tests
import utils
from optparse import OptionParser
import sys

parser = OptionParser()

parser.add_option("-p", "--num_threads", dest = "num_threads", type = int, default = 1, help = "Number of cores to use.")
parser.add_option("-m", "--min_og_group", dest = "min_og_group", type = int, default = 0, help = "For limiting the number of OG's examined. Don't analyze those with OG numbers less than this number.")
parser.add_option("-x", "--max_og_group", dest = "max_og_group", type = int, default = 100000, help = "For limiting the number of OG's examined. Don't analyze those with OG numbers more than this number.")
parser.add_option("-o", "--prefix", dest = "prefix", type = str, help = "String used at the beginning of output directories and files.")
parser.add_option("-b", "--base_dir", dest = "base_dir", type = str, help = "Output directory.")
parser.add_option("-t", "--min_taxa", dest = "min_taxa", type = int)
parser.add_option("-r", "--ortho_file", dest = "ortho_file", type = str, default = "/Genomics/kocherlab/berubin/annotation/orthology/proteinortho3.proteinortho", help = "File of orthologous groups.")
parser.add_option("-e", "--tree_file", dest = "tree_file", type = str, default = "/Genomics/kocherlab/berubin/annotation/orthology/sc_15_taxa/RAxML_bestTree.sc_15_taxa_100_genes.tree", help = "Phylogeny of species examined.")
parser.add_option("-a", "--action", dest = "action", type = str, help = "Analysis to do", default = "paml")
parser.add_option("-i", "--inspecies", dest = "inspecies", type = str, help = "Species of interest in pairwise analyses")
parser.add_option("-u", "--outspecies", dest = "outspecies", type = str, help = "Outgroup species in pairwise analyses")
parser.add_option("-y", "--timetree", dest = "timetree", type = str, help = "Time calibrated tree")
parser.add_option("-f", "--noncoding_seq_type", dest = "noncoding_seq_type", type = str, help = "Sequence type for HKA tests (flank, intron, or first_intron)", default = "flank")
(options, args) = parser.parse_args()

STOP_CODONS = ["TAA", "TAG", "TGA"]
SPECIES_LIST = ["APUR", "HLIG", "LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LALB", "AAUR", "LFIG", "HRUB", "AVIR"]
SOCIAL = ["HLIG","LMAL", "LMAR", "LPAU", "LZEP", "AAUR", "LCAL", "LALB"]
REV_SOLITARY = ["APUR", "LLEU", "LOEN", "LVIE", "LFIG"]
ANC_SOLITARY = ["Nmel", "Dnov"]
POLYMORPHIC = ["LCAL", "LALB"]

def main():
    if not os.path.isdir(options.base_dir):
        os.mkdir(options.base_dir)     #create working directory
#    ortho_dic = utils.ortho_reader("/Genomics/kocherlab/berubin/annotation/orthology/proteinortho3.proteinortho")
    ortho_dic = utils.ortho_reader(options.ortho_file)
    index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
    if options.action == "mk":
        utils.mk_test(options.inspecies, options.outspecies, ortho_dic, "%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_dummy_ancestral" % (options.base_dir, options.prefix), options.base_dir, options.num_threads, options.min_taxa)
        sys.exit()
    if options.action == "hka":
        utils.hka_test(options.inspecies, options.outspecies, options.noncoding_seq_type, ortho_dic, options.base_dir, options.num_threads,"%s/%s_prank" % (options.base_dir, options.prefix))    
        sys.exit()
    if options.action == "godatabase":
        ipr_taxa_list = ["AAUR", "APUR", "AVIR", "HLIG", "HRUB", "LCAL", "LFIG", "LLEU", "LMAL"]
        ipr_taxa_list = ["GNIG"]
        ipr_taxa_list = ["BNIG"]
        utils.make_go_database(ortho_dic, ipr_taxa_list, "%s/%s_thoe_missing" % (options.base_dir, options.prefix))
        sys.exit()
    if options.action == "align_coding":
        paras_allowed = True #because of the way that write_orthos works
        #this just allows OG's to be aligned even if they initially include
        #paralogs because the paralogous sequences are excluded 
        #in write_orthos
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        iscoding = True
#        target_taxa = ["AMEL"]
#        cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
#        cur_og_list = [368]
        cur_og_list = og_list
        utils.fsa_coding_align(cur_og_list, "%s/%s_orthos/" % (options.base_dir, options.prefix), "%s/%s_fsa_coding" % (options.base_dir, options.prefix), options.num_threads, iscoding)
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        utils.concatenate_for_raxml("%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s.afa" % (options.base_dir, options.prefix), og_list)
        sys.exit()
    if options.action == "write_orthos":
#        seq_dic = utils.get_cds("/Genomics/kocherlab/berubin/acacias/selection/gilliamella/renamed_cds_seqs")
        seq_dic = utils.get_cds("/Genomics/kocherlab/berubin/acacias/selection/bartonella/renamed_cds_seqs")
#        seq_dic = utils.get_cds("/Genomics/kocherlab/berubin/ants/termites/termite_data/renamed_cds_seqs")
#        seq_dic = utils.get_cds("/Genomics/kocherlab/berubin/ants/bees/bee_data/renamed_cds_seqs")
#        seq_dic = utils.get_bee_cds()
#        seq_dic = utils.get_bumble_cds()
        #write fastas for ALL orthologous groups
        utils.write_orthos(options.ortho_file, seq_dic, True, "%s/%s_orthos" % (options.base_dir, options.prefix), index_file)
        sys.exit()
    if options.action == "yn_dnds":
        paras_allowed = True
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        utils.yn_estimates(og_list, "%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_yn" % (options.base_dir, options.prefix))
        sys.exit()
    if options.action == "free_ratios":
        test_type = "free"
        foreground = "free"
#        get_dn_ds = True
        get_dn_ds = False
        paras_allowed = True
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
#        og_list = [3050]
#        utils.paml_test(og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
#        target_taxa = ["AAUR", "APUR", "LFIG", "LMAR", "LZEP", "LVIE", "LOEN", "LPAU"]
#        cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
        cur_og_list = og_list
        utils.read_frees("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s_results" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type), get_dn_ds, options.timetree, cur_og_list)
        sys.exit()
    if options.action == "termfinder":
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
        utils.og_list_termfinder("%s/../thoe_missing.txt" % (options.base_dir), "%s/../bnig_present.txt" % (options.base_dir),"%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/thoe_missing_go" % (options.base_dir), -9, -9)


        sys.exit()
    if options.action == "gene_trees":
        constrained = True
        paras_allowed = True
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        cur_og_list = og_list

#        utils.gene_trees(cur_og_list, "%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_gene_trees_constrained" % (options.base_dir, options.prefix), constrained, options.tree_file, options.num_threads)
        target_taxa = ["LALB", "DNOV"]
        cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)

        utils.compile_gene_trees(cur_og_list, "%s/%s_gene_trees_constrained" % (options.base_dir, options.prefix), options.tree_file, "%s/%s_gene_trees_bls" % (options.base_dir, options.prefix))
        sys.exit()
    if options.action == "clark_test":
        test_type = "aaml_blengths"
        foreground = "aaml_blengths"
        paras_allowed = True
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        cur_og_list = utils.min_taxa_membership(ortho_dic, ["LOEN", "LVIE", "LFIG", "APUR"], og_list, 3)
#        cur_og_list = [368]
#        utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
        phylo_dic = utils.read_aaml_phylos(cur_og_list, "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s_results" % (options.base_dir, options.prefix, foreground, test_type))
        fore_list = ["APUR", "LFIG", "LOEN", "LVIE", "LLEU"]
    #    fore_list = ["AAUR", "LMAR", "LPAU", "LZEP", "LMAL", "HLIG"]
        exclude_list = ["LCAL", "LALB", "HRUB", "NMEL", "DNOV", "AVIR"]
    #    utils.marine_test(phylo_dic, options.tree_file, fore_list, exclude_list, "%s/%s_%s_%s_results" % (options.base_dir, options.prefix, foreground, test_type))

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
        foreground = "BNIG"
#        foreground = "GNIG"
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        target_taxa = ["BNIG", "THOE"]
#        target_taxa = ["GNIG", "GAPI"]
        cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
        utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
        utils.test_lrt_branch("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s_%s_%s.lrt" % (options.base_dir, options.prefix, foreground, test_type), "%s/%s.gaf" % (options.base_dir, options.prefix), ortho_dic, "%s/%s_%s_%s_go" % (options.base_dir, options.prefix, foreground, test_type))
        sys.exit()

    if options.action == "bs_test":
        test_type = "bs"
#        foreground = "BNIG"
        foreground = "GNIG"
        og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
        target_taxa = ["GNIG", "GAPI"]
#        target_taxa = ["BNIG", "THOE"]
        cur_og_list = utils.target_taxa_in_og(ortho_dic, target_taxa, og_list)
        utils.paml_test(cur_og_list, foreground, test_type,"%s/%s_fsa_coding" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
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
