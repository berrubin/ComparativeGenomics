import multiprocessing
from multiprocessing import Pool
from potpour import Worker
import copy
import vcf
from Bio import SeqIO
from Bio import Seq
from Gene import Gene
from gff import gffParser
import pickle
import os
import shutil
import subprocess
from glob import glob
from Bio.Phylo.PAML import codeml
#from Bio.Seq import Seq
import ete3
from ete3 import PhyloTree
from Bio.Phylo.PAML.chi2 import cdf_chi2
from scipy.stats import chisqprob
import statsmodels.stats.multitest as smm


def ortho_reader(orthofile):
    #returns dictionary of dictionaries of orthologos genes. 
    #Lower level dictionaries have keys of species codes and 
    #values of gene names. Upper level dictionaries have the
    #OG index (line number in orthofile) as keys.
    reader = open(orthofile, 'rU')
    ortho_dic = {}
    counter = 0
    for line in reader:
        if line.startswith("#"):
            continue
        gene_dic = {}
        cur_line = line.split()
        taxa_count = int(cur_line[0])
        seq_count = int(cur_line[1])
        for gene in cur_line[3:]:
            cur_species = gene[0:4]
            gene_dic[cur_species] = gene
        ortho_dic[counter] = gene_dic
        counter += 1
    return ortho_dic
        
def get_species_data(target_species):
    seq_dic = {}
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/%s/%s_genome_v1.0.fasta" % (target_species, target_species), format = 'fasta')
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    cds_dic = {}
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (target_species, target_species), format = 'fasta')
    for rec in reader:
        cds_dic[rec.id[:-3]] = str(rec.seq)

#    gff_file = gffParser(open("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v1.0_longest_isoform_scaf_5_1M.gff3" % (target_species, target_species), 'rU'))
    gff_file = gffParser(open("/Genomics/kocherlab/berubin/official_release/%s/%s_testset.gff3" % (target_species, target_species), 'rU'))
    gene_dic = gff_file.geneDict()
    gene_objects = {}
    for gene_name in gene_dic.keys():
        gene_objects[gene_name] = Gene(gene_name, gene_dic[gene_name][0], gene_dic[gene_name][1])
    out_dic = {}
#    if data_type == "CDS":
    gene_objects = get_gene_coords(gff_file, gene_dic, gene_objects, seq_dic)
    gene_objects = get_cds_dic(gff_file, gene_dic, gene_objects)
    gene_objects = get_utr_dic(gff_file, gene_dic, "three", gene_objects)
    gene_objects = get_utr_dic(gff_file, gene_dic, "five", gene_objects)
    make_cds_sequences(gene_objects, cds_dic)
    return gene_objects

def make_cds_sequences(gene_objects, cds_dic):
    for gene in gene_objects.values():
        gene.get_cds_sequence(cds_dic)

def get_gene_coords(gff_file, gene_dic, gene_objects, seq_dic):
    for gene_name in gene_dic.keys():
        gene_deets = gff_file.getGene(gene_dic[gene_name][0], gene_name)[0]
        gene_objects[gene_name].set_start(gene_deets["start"])
        gene_objects[gene_name].set_end(gene_deets["end"])
        gene_objects[gene_name].set_sequence(seq_dic[gene_dic[gene_name][0]], 3000)
    return gene_objects

def get_utr_dic(gff_file, gene_dic, prime_end, gene_objects):
    utr_tuples = {}
    for gene_name in gene_dic.keys():
        mrna_list = gff_file.getmRNA(gene_dic[gene_name][0], gene_name) #only one mRNA because working with longest iso
        for mrna_dic in mrna_list:
            if prime_end == "five":
                utr_list = gff_file.getFivePrimeUTR(gene_dic[gene_name][0], mrna_dic["Name"])
            elif prime_end == "three":
                utr_list = gff_file.getThreePrimeUTR(gene_dic[gene_name][0], mrna_dic["Name"])
            tuple_list = []
            for utr_dic in utr_list:
                tuple_list.append((utr_dic["start"], utr_dic["end"]))
#            utr_tuples[mrna_dic["Name"]] = tuple_list
            gene_objects[gene_name].add_utrs(tuple_list, prime_end)
    return gene_objects

def get_introns_dic(gff_file, gene_dic, gene_objects):
    intron_dic = {}
    cds_dic = get_cds_dic(gff_file, gene_dic)
    for gene_name in cds_dic.keys():
        cds_tuples = []
        sorted_tuples = cds_dic[gene_name]
        intron_list = []
        for x in range(len(sorted_tuples)-1):
            cur_intron = (sorted_tuples[x][1], sorted_tuples[x+1][0])
            intron_list.append(cur_intron)
        intron_dic[gene_name] = intron_list
        print intron_list
    return intron_dic

def get_cds_dic(gff_file, gene_dic, gene_objects):
    out_dic = {}
    for gene_name in gene_dic.keys():
        mrna_list = gff_file.getmRNA(gene_dic[gene_name][0], gene_name) #only one mRNA because working with longest iso
        for mrna_dic in mrna_list:
            cds_dic = gff_file.getCDS(gene_dic[gene_name][0], mrna_dic["Name"])
            cds_tuples = []
            for cds in cds_dic:
                cds_tuples.append((cds["start"], cds["end"]))
            sorted_tuples = sorted(cds_tuples, key = lambda tup: tup[0])
#            cds_list = []
#            out_dic[mrna_dic["Name"]] = sorted_tuples
            gene_objects[gene_name].add_cds(sorted_tuples)
    return gene_objects

def gene_vcf_dic(species):
    if os.path.exists("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (species, species)):
        print "Reading %s pickle" % species
        gene_dic = pickle.load(open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (species, species), 'rb'))
    else:
        gene_dic = get_species_data(species)
        reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_testset.vcf.gz" % (species, species))
        for gene_name, gene_object in gene_dic.items():
            gene_object.get_genotypes(reader, 3000)
        pickle.dump(gene_dic, open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (species, species), 'wb'))
    return gene_dic
    

def mk_test(inspecies, outspecies, align_dir):
    ortho_dic = ortho_reader("/Genomics/kocherlab/berubin/annotation/orthology/proteinortho3.proteinortho")
    in_gene_dic = gene_vcf_dic(inspecies)
    out_gene_dic = gene_vcf_dic(outspecies)
#    out_gene_dic = get_species_data(outspecies)
#    reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_testset.vcf.gz" % (outspecies, outspecies))
#    for gene_name, gene_object in out_gene_dic.items():
#        gene_object.get_genotypes(reader, 3000)

    mk_table = open("%s_%s_mk_table.txt" % (inspecies, outspecies), 'w')
    mk_table.write("geneID\tPR\tFR\tPS\tFS\tTsil\tTrepl\tnout\tnpop\n")
    for og_num in ortho_dic.keys():
        if ortho_dic[og_num].keys().count(inspecies) == 1 and ortho_dic[og_num].keys().count(outspecies) == 1:
            if og_num != 2110:
                continue
            reader = SeqIO.parse("%s/og_cds_%s.1.fas" % (align_dir, og_num), format = 'fasta')
            for rec in reader:
                if rec.id[0:4] == inspecies:
                    inseq = str(rec.seq)
                    inseq_name = rec.id[:-3]
                elif rec.id[0:4] == outspecies:
                    outseq = str(rec.seq)
                    outseq_name = rec.id[:-3]
            fix_syn, fix_nsyn = count_sub_types(inseq, outseq)
            in_poly_syn = in_gene_dic[inseq_name].syn_count
            in_poly_nsyn = in_gene_dic[inseq_name].nsyn_count
            in_potent_syn = in_gene_dic[inseq_name].potent_syn
            in_potent_nsyn = in_gene_dic[inseq_name].potent_nsyn
            in_n_size = in_gene_dic[inseq_name].average_n
            mk_table.write("\t".join([str(og_num), str(in_poly_nsyn), str(fix_nsyn), str(in_poly_syn), str(fix_syn), str(in_potent_syn), str(in_potent_nsyn), "1", str(in_n_size)]) + "\n")
    mk_table.close()


def check_for_stop(seq):
    #looks for stop codons in a coding sequence (in that frame)
    x = 0
    while x < len(seq):
        if seq[x:x+3] in STOP_CODONS:
            return x
        x += 3
    return False

def trim_phylo(taxa_list, fore_list, orthogroup, outdir, phylogeny_file):
    #trims taxa and adds foreground tags for PAML analysis tree
#    tree = PhyloTree(options.tree_file)
    tree = PhyloTree(phylogeny_file)
    tree.prune(taxa_list)
    tree.unroot()
    tree_str = tree.write(format = 5)
    for tax in fore_list:
        tree_str = tree_str.replace(tax, "%s#1" % tax)
    outfile = open("%s/og_%s.tree" % (outdir, orthogroup), 'w')
    outfile.write(tree_str)
    outfile.close()

def rename_tree(seqfile, outname, phylogeny_file):
    #trims taxa and renames tips to match sequence names (mostly for prank)
    name_dic = {}
    reader = SeqIO.parse(seqfile, format = 'fasta')
    for rec in reader:
        name_dic[rec.id[0:4]] = rec.id
    tree = PhyloTree(phylogeny_file)
    tree.prune(name_dic.keys())
    tree_str = tree.write(format = 5)
    for k, v in name_dic.items():
        tree_str = tree_str.replace(k, v)
    outfile = open(outname, 'w')
    outfile.write(tree_str)
    outfile.close()

def prep_paml_files(orthogroup, indir, outdir, foreground, phylogeny_file):
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
        trim_phylo(taxa_list, fore_list, orthogroup, outdir, phylogeny_file)
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

def paml_test(og_list, foreground, test_type, indir, outdir, phylogeny_file, num_threads):
    #performs paml test on all OG's in list
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    og_file_list = []
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()
    for cur_og in og_list:
        if test_type == "model_d":
            prep_paml_files(cur_og, indir, outdir, "model_d", phylogeny_file)
        elif test_type == "free":
            prep_paml_files(cur_og, indir, outdir, "free", phylogeny_file)
        else:
            prep_paml_files(cur_og, indir, outdir, foreground, phylogeny_file)
        work_queue.put([cur_og, outdir])
    jobs = []
    for i in range(num_threads):
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
        n_total_count += 1
        if seq1[x] == seq2[x]:
            n_same_count += 1
        elif seq1[x] != seq2[x]:
            n_diff_count += 1
    p_seq1 = str(Seq.Seq(seq1.replace("-", "N")).translate())
    p_seq2 = str(Seq.Seq(seq2.replace("-", "N")).translate())
    p_same_count = 0
    p_diff_count = 0
    p_total_count = 0
    for x in range(len(p_seq1)):
        if p_seq1[x] in empty_chars or p_seq2[x] in empty_chars:
            continue
        p_total_count += 1
        if p_seq1[x] == p_seq2[x]:
            p_same_count += 1
        elif p_seq1[x] != p_seq2[x]:
            p_diff_count += 1
    nsyns = p_diff_count
    syns = n_diff_count - nsyns
    return syns, nsyns


def prank_align_worker(og_file, outdir, use_backbone, phylogeny_file):
    #the worker method for multiprocessing the prank alignments
    cur_og = og_file.split("/")[-1]
    og_num = cur_og.split("_")[2].split(".fa")[0]
    if use_backbone:
        rename_tree(og_file, "%s/og_%s.tree" % (outdir, og_num), phylogeny_file)
        cmd = ["/Genomics/kocherlab/berubin/local/src/prank/prank", "-d=%s" % og_file, "-o=%s/og_cds_%s" % (outdir, og_num), "-codon", "-F", "-t=%s/og_%s.tree" % (outdir,og_num)]
        subprocess.call(cmd)
        gblock("%s/og_cds_%s.1.fas" % (outdir, og_num))
    else:
        cmd = ["/Genomics/kocherlab/berubin/local/src/prank/prank", "-d=%s" % og_file, "-o=%s/og_cds_%s" % (outdir, og_num), "-codon", "-F"]
        subprocess.call(cmd)


def gblock(inalignment):
    #run gblocks on given file
    cmd = ["/Genomics/kocherlab/berubin/local/src/Gblocks_0.91b/Gblocks", inalignment, "-t=c", "-b5=h"]
    subprocess.call(cmd)

def prank_align(og_list, indir, outdir, use_backbone, phylogeny_file, num_threads): 
    #run prank alignments
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    og_file_list = []
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()
    for og in og_list:
        og_file = "%s/og_cds_%s.fa" % (indir, og)
        work_queue.put([og_file, outdir, use_backbone, phylogeny_file])
    jobs = []
    for i in range(num_threads):
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

def write_orthos(ortho_file, seq_dic, paras_allowed, outdir, indexfile):
    #read/parse orthology file and write files containing all sequences
    #also create an index file that lists the number of taxa in each OG
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
#    ref_file = open("%s/%s_ortho.index" % (options.base_dir, options.prefix), 'w')
    ref_file = open(indexfile, 'w')
    ref_file.write("#og\tnum_taxa\tnum_paras\n")
    counter = 0
    reader = open(ortho_file, 'rU')
    for line in reader:
        if line.startswith("#"):
            continue
        cur_line = line.split()
        num_paras = int(cur_line[1]) - int(cur_line[0])
        if int(cur_line[0]) != int(cur_line[1]):
            if not paras_allowed:
                continue
        ref_file.write("%s\t%s\t%s\n" % (counter, cur_line[0], num_paras))

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

def read_ortho_index(index_file, min_taxa, paras_allowed):
    #get list of all of the OG's with the minumum taxa
#    reader = open("%s/%s_ortho.index" % (options.base_dir, options.prefix), 'rU')
    reader = open(index_file, 'rU')
    og_list = []
    for line in reader:
        if line.startswith("#"):
            continue
        cur_line = line.split()
        if not paras_allowed:
            if int(cur_line[2]) > 0:
                continue
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

            
