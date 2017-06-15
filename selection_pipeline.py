#!/usr/local/bin/python
import os
import paml_tests
import utils
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-p", "--num_threads", dest = "num_threads", type = int, default = 1, help = "Number of cores to use.")
parser.add_option("-m", "--min_og_group", dest = "min_og_group", type = int, help = "For limiting the number of OG's examined. Don't analyze those with OG numbers less than this number.")
parser.add_option("-x", "--max_og_group", dest = "max_og_group", type = int, help = "For limiting the number of OG's examined. Don't analyze those with OG numbers more than this number.")
parser.add_option("-o", "--prefix", dest = "prefix", type = str, help = "String used at the beginning of output directories and files.")
parser.add_option("-b", "--base_dir", dest = "base_dir", type = str, help = "Output directory.")
parser.add_option("-t", "--min_taxa", dest = "min_taxa", type = int)
parser.add_option("-r", "--ortho_file", dest = "ortho_file", type = str, default = "/Genomics/kocherlab/berubin/annotation/orthology/proteinortho3.proteinortho", help = "File of orthologous groups.")
parser.add_option("-e", "--tree_file", dest = "tree_file", type = str, default = "/Genomics/kocherlab/berubin/annotation/orthology/sc_15_taxa/RAxML_bestTree.sc_15_taxa_100_genes.tree", help = "Phylogeny of species examined.")
(options, args) = parser.parse_args()

STOP_CODONS = ["TAA", "TAG", "TGA"]
SPECIES_LIST = ["APUR", "HLIG", "LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LALB", "Nmel", "Dnov", "AAUR", "LFIG"]
SOCIAL = ["HLIG","LMAL", "LMAR", "LPAU", "LZEP", "AAUR", "LCAL", "LALB"]
REV_SOLITARY = ["APUR", "LLEU", "LOEN", "LVIE", "LFIG"]
ANC_SOLITARY = ["Nmel", "Dnov"]
POLYMORPHIC = ["LCAL", "LALB"]

def main():
    if not os.path.isdir(options.base_dir):
        os.mkdir(options.base_dir)     #create working directory
#    ortho_dic = utils.ortho_reader("/Genomics/kocherlab/berubin/annotation/orthology/proteinortho3.proteinortho")
    ortho_dic = utils.ortho_reader(options.ortho_file)
#    utils.mk_test("LALB", "LCAL", ortho_dic, "%s/%s_prank_no_backbone" % (options.base_dir, options.prefix), options.base_dir, options.num_threads)
#    utils.hka_test("LZEP", "LVIE", "flank", ortho_dic, options.base_dir, options.num_threads,"%s/%s_prank" % (options.base_dir, options.prefix))    
#    utils.hka_test("LPAU", "LLEU", "flank", ortho_dic, options.base_dir, options.num_threads,"%s/%s_prank" % (options.base_dir, options.prefix))
#    utils.hka_test("LALB", "LCAL", "flank", ortho_dic, options.base_dir, options.num_threads,"%s/%s_prank_no_backbone" % (options.base_dir, options.prefix))
#    utils.hka_test("LLEU", "LMAL", "flank", ortho_dic, options.base_dir, options.num_threads)
    
    seq_dic = utils.get_cds() #get coding sequences from all species
    #store name of orthology index file
    index_file = "%s/%s_ortho.index" % (options.base_dir, options.prefix)
    #write fastas for ALL orthologous groups
#    utils.write_orthos(options.ortho_file, seq_dic, True, "%s/%s_orthos" % (options.base_dir, options.prefix), index_file)
    paras_allowed = False #do not include OG's with paralogs
    #get a list of the first 100 genes with min_taxa (should always be all
    #of the available species) for making the phylogeny
    og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)[0:100]
    use_backbone = False #do not use a phylogeny to align these first 100 genes
    #align first 100 genes
#    utils.prank_align(og_list, "%s/%s_orthos/" % (options.base_dir, options.prefix), "%s/%s_prank_no_backbone" % (options.base_dir, options.prefix), use_backbone, "nophylogeny", options.num_threads)
    #concatenate the aligned 100 genes
#    utils.concatenate_for_raxml("%s/%s_prank_no_backbone" % (options.base_dir, options.prefix), "%s/%s.afa" % (options.base_dir, options.prefix))

    #then run raxml to create a backbone phylogeny
    #raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s %s/%s.afa % (options.base_dir, options.prefix) -T 20
    #get list of all genes without paralogs
    
    og_list = utils.read_ortho_index(index_file, options.min_taxa, paras_allowed)
    og_list = utils.limit_list(og_list, options.min_og_group, options.max_og_group)

    use_backbone = True #use phylogeny to align these genes
    #align non-paralogous genes
#    utils.prank_align(og_list, "%s/%s_orthos/" % (options.base_dir, options.prefix), "%s/%s_prank" % (options.base_dir, options.prefix), use_backbone, options.tree_file, options.num_threads)
    test_type = "ancestral"
    foreground = "dummy"
    utils.paml_test(og_list, foreground, test_type,"%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type), options.tree_file, options.num_threads)
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
