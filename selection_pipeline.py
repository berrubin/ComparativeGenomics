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

parser = OptionParser()

parser.add_option("-p", "--num_threads", dest = "num_threads", type = int)
parser.add_option("-m", "--min_og_group", dest = "min_og_group", type = int)
parser.add_option("-x", "--max_og_group", dest = "max_og_group", type = int)
parser.add_option("-o", "--prefix", dest = "prefix", type = str)
parser.add_option("-b", "--base_dir", dest = "base_dir", type = str)
parser.add_option("-t", "--min_taxa", dest = "min_taxa", type = int)
parser.add_option("-r", "--ortho_file", dest = "ortho_file", type = str, defauly = "/Genomics/kocherlab/berubin/annotation/orthology/proteinortho3.proteinortho")
parser.add_option("-e", "--tree_file", dest = "tree_file", type = str, default = "/Genomics/kocherlab/berubin/annotation/orthology/sc_15_taxa/RAxML_bestTree.sc_15_taxa_100_genes.tree")
(options, args) = parser.parse_args()

STOP_CODONS = ["TAA", "TAG", "TGA"]
SPECIES_LIST = ["APUR", "HLIG", "LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LALB", "Nmel", "Dnov", "AAUR", "LFIG"]
SOCIAL = ["HLIG","LMAL", "LMAR", "LPAU", "LZEP", "AAUR", "LCAL", "LALB"]
REV_SOLITARY = ["APUR", "LLEU", "LOEN", "LVIE", "LFIG"]
ANC_SOLITARY = ["Nmel", "Dnov"]
POLYMORPHIC = ["LCAL", "LALB"]

def check_for_stop(seq):
    #looks for stop codons in a coding sequence (in that frame)
    x = 0
    while x < len(seq):
        if seq[x:x+3] in STOP_CODONS:
            return x
        x += 3
    return False

def trim_phylo(taxa_list, fore_list, orthogroup, outdir):
    #trims taxa and adds foreground tags for PAML analysis tree
    tree = PhyloTree(options.tree_file)
    tree.prune(taxa_list)
    tree.unroot()
    tree_str = tree.write(format = 5)
    for tax in fore_list:
        tree_str = tree_str.replace(tax, "%s#1" % tax)
    outfile = open("%s/og_%s.tree" % (outdir, orthogroup), 'w')
    outfile.write(tree_str)
    outfile.close()

def rename_tree(seqfile, outname):
    #trims taxa and renames tips to match sequence names (mostly for prank)
    name_dic = {}
    reader = SeqIO.parse(seqfile, format = 'fasta')
    for rec in reader:
        name_dic[rec.id[0:4]] = rec.id
    tree = PhyloTree(options.tree_file)
    tree.prune(name_dic.keys())
    tree_str = tree.write(format = 5)
    for k, v in name_dic.items():
        tree_str = tree_str.replace(k, v)
    outfile = open(outname, 'w')
    outfile.write(tree_str)
    outfile.close()

def prep_paml_files(orthogroup, indir, outdir, foreground):
    #formats fasta and tree files for PAML analysis
    tree_prep = True
    fore_list = []
    if foreground == "social":
        fore_list = SOCIAL
    if foreground == "solitary":
        fore_list = SOLITARY
    if foreground == "model_d":
        tree_prep = False
    reader = SeqIO.parse("%s/og_cds_%s.1.fas-gb" % ( indir, orthogroup), format = 'fasta')
    seq_dic = {}
    taxa_list = []
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
        taxa_list.append(rec.id[0:4])
    outfile = open("%s/og_cds_%s.afa" % (outdir, orthogroup), 'w')
    outfile.write("%s %s\n" % (len(seq_dic), len(rec.seq)))
    for species, sequence in seq_dic.items():
        outfile.write("%s\n%s\n" % (species[0:4], sequence))
    outfile.close()
    if tree_prep:
        trim_phylo(taxa_list, fore_list, orthogroup, outdir)
    else:
        shutil.copy("/Genomics/kocherlab/berubin/annotation/orthology/model_d.tree", "%s/og_%s.tree" % (outdir, orthogroup))

 
def test_lrt(indir):
    #performs all lrts on the files in a given directory
    #also performs multiple test correction to return adjusted p-value
    p_dic = {}
    p_list = []
    og_list = []
    for og_file in glob("%s/og_*.nul" % (indir)):
        cur_og = int(og_file.split("og_")[1].split(".nul")[0])
        pval = lrt("%s/og_%s.alt" % (indir, cur_og), "%s/og_%s.nul" % (indir, cur_og))
        print "%s: %s" % (cur_og, pval)
        p_list.append(pval)
        og_list.append(cur_og)
    pval_corr = smm.multipletests(p_list, alpha = 0.05, method = 'fdr_i')[1]
    for x in range(pval_corr):
        p_dic[og_list[x]] = pval_corr[x]
    return p_dic

def lrt(alt_file, null_file):
    #lrt for PAML tests
    reader = open(alt_file, 'rU')
    for line in reader:
        if line.startswith("lnL"):
            alt_likely = float(line.split()[4])
    reader = open(null_file, 'rU')
    for line in reader:
        if line.startswith("lnL"):
            nul_likely = float(line.split()[4])
    p = chisqprob(2*(alt_likely - nul_likely), 1)
    return p

def read_frees(indir):
    #reads free ratios files and gets dn/ds ratios
    #can be easily extended to get dn and ds but those are low quality
    og_dnds_dic = {}
    for og_file in glob("%s/og_*.alt" % (indir)):        
        cur_og = int(og_file.split("og_")[1].split(".alt")[0])
        reader = open(og_file, 'rU')
        dnds_tree = False
        dnds_dic = {}
        for line in reader:
            if dnds_tree:
                dnds = PhyloTree(line.strip().replace("#", ":"))
                print dnds
                for leaf in dnds:
                    dnds_dic[leaf.name] = leaf.dist
                dnds_tree = False
            if line.strip() == "w ratios as labels for TreeView:":
                dnds_tree = True
                continue
        og_dnds_dic[cur_og] = dnds_dic

def paml_test(og_list, foreground, test_type, indir, outdir):
    #performs paml test on all OG's in list
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    og_file_list = []
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()
    for cur_og in og_list:
        if test_type == "model_d":
            prep_paml_files(cur_og, indir, outdir, "model_d")
        elif test_type == "free":
            prep_paml_files(cur_og, indir, outdir, "free")
        else:
            prep_paml_files(cur_og, indir, outdir, foreground)
        work_queue.put([cur_og, outdir])
    jobs = []
    for i in range(options.num_threads):
        if test_type == "bs":
            worker = Worker(work_queue, result_queue, paml_tests.branch_site_worker)
        elif test_type == "branch":
            worker = Worker(work_queue, result_queue, paml_tests.branch_worker)
        elif test_type == "free":
              worker = Worker(work_queue, result_queue, paml_tests.free_ratios_worker)          
        elif test_type == "model_d":
              worker = Worker(work_queue, result_queue, paml_tests.branch_site_d_worker)          
        jobs.append(worker)
        worker.start()
    try:
        for j in jobs:
            j.join()
    except KeyboardInterrupt:
        for j in jobs:
            j.terminate()
            j.join()
    
def count_sub_types(seq1, seq2):
    #counts syn and nsyn substitutions between two coding sequences
    #assumes that they are aligned
    empty_chars = ["N", "-", "X"]
    n_same_count = 0
    n_diff_count = 0
    n_total_count = 0
    for x in range(len(seq1)):
        if seq1[x] in empty_chars or seq2[x] in empty_chars:
            continue
        total_count += 1
        if seq1[x] == seq2[x]:
            same_count += 1
        elif seq1[x] != seq2[x]:
            diff_count += 1
    p_seq1 = seq1.replace("-", "N").translate()
    p_seq2 = seq2.replace("-", "N").translate()
    p_same_count = 0
    p_diff_count = 0
    p_total_count = 0
    for x in range(len(seq1)):
        if seq1[x] in empty_chars or seq2[x] in empty_chars:
            continue
        p_total_count += 1
        if seq1[x] == seq2[x]:
            p_same_count += 1
        elif seq1[x] != seq2[x]:
            p_diff_count += 1
    nsyns = p_diff_count
    syns = n_diff_count - nsyns
    return syns, nsyns


def prank_align_worker(og_file, outdir, use_backbone):
    #the worker method for multiprocessing the prank alignments
    cur_og = og_file.split("/")[-1]
    og_num = cur_og.split("_")[2].split(".fa")[0]
    if use_backbone:
        rename_tree(og_file, "%s/og_%s.tree" % (outdir, og_num))
        cmd = ["/Genomics/kocherlab/berubin/local/src/prank/prank", "-d=%s" % og_file, "-o=%s/og_cds_%s" % (outdir, og_num), "-codon", "-t=%s/og_%s.tree" % (outdir,og_num)]
        subprocess.call(cmd)
        gblock("%s/og_cds_%s.1.fas" % (outdir, og_num))
    else:
        cmd = ["/Genomics/kocherlab/berubin/local/src/prank/prank", "-d=%s" % og_file, "-o=%s/og_cds_%s" % (outdir, og_num), "-codon"]
        subprocess.call(cmd)


def gblock(inalignment):
    #run gblocks on given file
    cmd = ["/Genomics/kocherlab/berubin/local/src/Gblocks_0.91b/Gblocks", inalignment, "-t=c", "-b5=h"]
    subprocess.call(cmd)

def prank_align(og_list, indir, outdir, use_backbone): 
    #run prank alignments
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    og_file_list = []
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()
    for og in og_list:
        og_file = "%s/og_cds_%s.fa" % (indir, og)
        work_queue.put([og_file, outdir, use_backbone])
    jobs = []
    for i in range(options.num_threads):
        worker = Worker(work_queue, result_queue, prank_align_worker)
        jobs.append(worker)
        worker.start()
    try:
        for j in jobs:
            j.join()
    except KeyboardInterrupt:
        for j in jobs:
            j.terminate()
            j.join()

def get_cds():
    #get dictionary containing CDS for all genomes
    seq_dic = {}
    for species in ["APUR", "HLIG", "LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LALB", "LFIG", "AAUR"]:
        reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (species, species), format = 'fasta')
        seq_dic[species] = {}
        for rec in reader:
            seq_dic[species][rec.id] = str(rec.seq)
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/annotation/reference_genomes/Dufourea_novaeangliae_v1.1.cds.fa", format = 'fasta')
    species = "Dnov"
    seq_dic[species] = {}
    for rec in reader:
        seq_dic[species][rec.id] = str(rec.seq)
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/data/nmel/Nmel_v1.0.cds.fa", format = 'fasta')
    species = "Nmel"
    seq_dic[species] = {}
    for rec in reader:
        seq_dic[species][rec.id] = str(rec.seq)
    return seq_dic

def write_orthos(ortho_file, seq_dic, min_og_size, paras_allowed, outdir):
    #read/parse orthology file and write files containing all sequences
    #also create an index file that lists the number of taxa in each OG
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    ref_file = open("%s/%s_taxa_%s.index" % (options.base_dir, options.prefix, min_og_size), 'w')
    ref_file.write("#og\tnum_taxa\n")
    counter = 0
    for line in ortho_file:
        if line.startswith("#"):
            continue
        cur_line = line.split()
        if int(cur_line[0]) < min_og_size:
            continue
        if int(cur_line[0]) > min_og_size:
            continue
        if int(cur_line[0]) != int(cur_line[1]):
            if not paras_allowed:
                continue
        ref_file.write("%s\t%s\n" % (counter, cur_line[0]))

        outfile = open("%s/og_cds_%s.fa" % (outdir, counter), 'w')
        for seqs in cur_line[3:]:
            if "*" in seqs:
                continue
            cur_seqs = seqs.split(",")
            for seq in cur_seqs:
                cur_species = seq[0:4]
                outfile.write(">%s\n%s\n" % (seq, seq_dic[cur_species][seq]))
        outfile.close()
        counter += 1
    ref_file.close()

def concatenate_for_raxml(input_dir, outfile):
    #take directory of alignments and concatenate them all into a
    #RAxML formatted fasta file
    full_dic = {}
    for species in SPECIES_LIST:
        full_dic[species] = []
    seq_len = 0
    for alignment in glob("%s/*2.fas" % input_dir):
        reader = SeqIO.parse(alignment, format = 'fasta')
        for rec in reader:
            cur_species = rec.id[0:4]
            full_dic[cur_species].append(str(rec.seq))
        seq_len += len(rec.seq)
    writer = open(outfile, 'w')
    writer.write("%s %s\n" % (len(SPECIES_LIST), seq_len))
    for species, seq_list in full_dic.items():
        writer.write("%s\n%s\n" % (species, "".join(seq_list)))
    writer.close()

def read_ortho_index(min_taxa):
    #get list of all of the OG's with the minumum taxa
    reader = open("%s/%s_taxa_%s.index" % (options.base_dir, options.prefix, options.min_taxa), 'rU')
    og_list = []
    for line in reader:
        if line.startswith("#"):
            continue
        cur_line = line.split()
        if int(cur_line[1]) >= min_taxa:
            og_list.append(int(cur_line[0]))
    return og_list

def limit_list(og_list, lower_bound, higher_bound):
    #limit the number of OG's to be examined
    new_list = []
    for og in og_list:
        if og < lower_bound:
            continue
        if og > higher_bound:
            continue
        else:
            new_list.append(og)
    return new_list

def main():
    if not os.path.isdir(options.base_dir):
        os.mkdir(options.base_dir)
    reader = open(options.ortho_file, 'rU')
    seq_dic = get_cds()
    write_orthos(reader, seq_dic, options.min_taxa, False, "%s/%s_taxa_%s/" % (options.base_dir, options.prefix, options.min_taxa))
    og_list = read_ortho_index(options.min_taxa)[:100]
    use_backbone = False
    prank_align(og_list, "%s/%s_taxa_%s/" % (options.base_dir, options.prefix, options.min_taxa), "%s/%s_prank_no_backbone" % (options.base_dir, options.prefix), use_backbone)
    concatenate_for_raxml("%s/%s_prank_no_backbone" % (options.base_dir, options.prefix), "%s/%s.afa" % (options.base_dir, options.prefix))
    #then run raxml to create a backbone phylogeny
    og_list = read_ortho_index(options.min_taxa)
#    og_list = limit_list(og_list, options.min_og_group, options.max_og_group)
    prank_align(og_list, "%s/%s_taxa_%s/" % (options.base_dir, options.prefix, options.min_taxa), "%s/%s_prank" % (options.base_dir, options.prefix))

    foreground = "social"
    test_type = "model_d"
    test_type = "bs"
#    test_type = "free"
#    paml_test(og_list, foreground, test_type,"%s/%s_prank" % (options.base_dir, options.prefix), "%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type))
#    read_frees("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type))
#    test_lrt("%s/%s_%s_%s" % (options.base_dir, options.prefix, foreground, test_type))

if __name__ == '__main__':
    main()
