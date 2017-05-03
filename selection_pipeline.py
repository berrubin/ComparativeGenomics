import os
import shutil
import multiprocessing
from multiprocessing import Pool
from potpour import Worker
from Bio import SeqIO
import subprocess
from optparse import OptionParser
from glob import glob
from Bio.Phylo.PAML import codeml
from Bio.Seq import Seq
import ete3
from ete3 import PhyloTree
from Bio.Phylo.PAML.chi2 import cdf_chi2
from scipy.stats import chisqprob
import statsmodels.stats.multitest as smm
import paml_tests
import utils

parser = OptionParser()

parser.add_option("-p", "--num_threads", dest = "num_threads", type = int, default = 1)
parser.add_option("-m", "--min_og_group", dest = "min_og_group", type = int)
parser.add_option("-x", "--max_og_group", dest = "max_og_group", type = int)
parser.add_option("-o", "--prefix", dest = "prefix", type = str)
parser.add_option("-b", "--base_dir", dest = "base_dir", type = str)
parser.add_option("-t", "--min_taxa", dest = "min_taxa", type = int)
parser.add_option("-r", "--ortho_file", dest = "ortho_file", type = str, default = "/Genomics/kocherlab/berubin/annotation/orthology/proteinortho3.proteinortho")
parser.add_option("-e", "--tree_file", dest = "tree_file", type = str, default = "/Genomics/kocherlab/berubin/annotation/orthology/sc_15_taxa/RAxML_bestTree.sc_15_taxa_100_genes.tree")
(options, args) = parser.parse_args()

STOP_CODONS = ["TAA", "TAG", "TGA"]
SPECIES_LIST = ["APUR", "HLIG", "LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LALB", "Nmel", "Dnov", "AAUR", "LFIG"]
SOCIAL = ["HLIG","LMAL", "LMAR", "LPAU", "LZEP", "AAUR", "LCAL", "LALB"]
REV_SOLITARY = ["APUR", "LLEU", "LOEN", "LVIE", "LFIG"]
ANC_SOLITARY = ["Nmel", "Dnov"]
POLYMORPHIC = ["LCAL", "LALB"]

def main():
    utils.mk_test("LMAL", "LLEU", "%s/%s_prank" % (options.base_dir, options.prefix))
    if not os.path.isdir(options.base_dir):
        os.mkdir(options.base_dir)
    reader = open(options.ortho_file, 'rU')
    seq_dic = utils.get_cds()
#    utils.write_orthos(reader, seq_dic, True, "%s/%s_orthos" % (options.base_dir, options.prefix))
    paras_allowed = False
#    og_list = utils.read_ortho_index(options.min_taxa, paras_allowed)[:100]
#    print og_list
    use_backbone = False
#    utils.prank_align(og_list, "%s/%s_taxa_%s/" % (options.base_dir, options.prefix, options.min_taxa), "%s/%s_prank_no_backbone" % (options.base_dir, options.prefix), use_backbone)
#    utils.concatenate_for_raxml("%s/%s_prank_no_backbone" % (options.base_dir, options.prefix), "%s/%s.afa" % (options.base_dir, options.prefix))
    #then run raxml to create a backbone phylogeny
    og_list = utils.read_ortho_index(options.min_taxa, paras_allowed)
    og_list = [2110]
#    og_list = utils.limit_list(og_list, options.min_og_group, options.max_og_group)
    use_backbone = True
    utils.prank_align(og_list, "%s/%s_orthos/" % (options.base_dir, options.prefix), "%s/%s_prank" % (options.base_dir, options.prefix), use_backbone)

    foreground = "social"
    test_type = "model_d"
    test_type = "bs"
#    test_type = "free"
#    paml_test(og_list, foreground, test_type,"%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type))
#    read_frees("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type))
#    test_lrt("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type))

if __name__ == '__main__':
    main()
