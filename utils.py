import paml_tests
import os
import sys
import changes
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import MafftCommandline
import multiprocessing
from multiprocessing import Pool
from potpour import Worker
import copy
import vcf
from Bio import SeqIO
from Bio import Seq
from Gene import Gene
from gff import gffParser
import cPickle as pickle
import os
import shutil
import subprocess
from glob import glob
from Bio.Phylo.PAML import codeml
import ete3
from ete3 import PhyloTree
from Bio.Phylo.PAML.chi2 import cdf_chi2
from scipy.stats import chisqprob
import statsmodels.stats.multitest as smm
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from StringIO import StringIO
import datetime

def ortho_reader(orthofile):
    #returns dictionary of dictionaries of orthologos genes. 
    #Lower level dictionaries have keys of species codes and 
    #values of lists of gene names. Upper level dictionaries have the
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
            if cur_species not in gene_dic.keys():
                gene_dic[cur_species] = []
            gene_dic[cur_species].append(gene)
        ortho_dic[counter] = gene_dic
        counter += 1
    return ortho_dic

def get_og_num(query_gene, ortho_dic):
    species = query_gene[0:4]
    query_gene = query_gene + "-RA"
    for og_num, species_dic in ortho_dic.items():
        if species in species_dic.keys():
            if query_gene in species_dic[species]:
                
                return og_num
    return False

def make_og_gene_map(ortho_dic):
    og_map = {}
    for og_num, species_dic in ortho_dic.items():
        for species, genes in species_dic.items():
            for gene in genes:
                og_map[gene[0:-3]] = og_num
    return og_map

def read_species_pickle(target_species):
#    pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle_test" % (target_species, target_species), 'rb')
    pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species), 'rb')
    gene_objects = pickle.load(pickle_file)
    pickle_file.close()
    return gene_objects

def get_species_data(target_species, num_threads):
    #This is a big method. Harvests gene coordinates from GFF3 files
    #and creates Gene objects with all of their characteristics.
#    if os.path.exists("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle_test" % (target_species, target_species)):
    if os.path.exists("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species)):
        print "%s pickle exists" % target_species
        return
#        pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species), 'rb')

        
#        gene_objects = pickle.load(pickle_file)
#        pickle_file.close()
#        return gene_objects
    print "Building %s pickle" % target_species
    official_dir = "/Genomics/kocherlab/berubin/official_release"
    seq_dic = {}
    if target_species == "LALB":
        reader = SeqIO.parse("%s/%s/%s_v3/%s/%s_genome_v3.0.fasta" % (official_dir, target_species, target_species, target_species, target_species), format = 'fasta')
    else:
        reader = SeqIO.parse("%s/%s/%s_genome_v1.0.fasta" % (official_dir, target_species, target_species), format = 'fasta')
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    cds_dic = {}
    if target_species == "LALB":
        reader = SeqIO.parse("%s/%s/%s_v3/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (official_dir, target_species, target_species, target_species, target_species), format = 'fasta')
    else:
        reader = SeqIO.parse("%s/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (official_dir, target_species, target_species), format = 'fasta')

    for rec in reader:
        cds_dic[rec.id[:-3]] = str(rec.seq)

#    gff_file = gffParser(open("%s/%s/%s_testset.gff3" % (official_dir, target_species, target_species), 'rU'))
    if target_species == "LALB":
        gff_file = gffParser(open("%s/%s/%s_v3/%s/%s_OGS_v1.0_longest_isoform.gff3" % (official_dir, target_species, target_species, target_species, target_species), 'rU'))
    else:
        gff_file = gffParser(open("%s/%s/%s_OGS_v1.0_longest_isoform.gff3" % (official_dir, target_species, target_species), 'rU'))
    print "read gff"
#    gff_file = gffParser(open("%s/%s/%s_09600.gff3" % (official_dir, target_species, target_species), 'rU'))
    print "get gene data"
    print str(datetime.datetime.now())
    gene_dic = gff_file.geneDict()
    gene_objects = {}
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    name_list = []
    dic_list = []
    gff_list = []
    seq_list = []
    species_list = []
    for gene_name in gene_dic.keys():
#        if gene_name == "LCAL_00010": #""LMAL_00737" or gene_name == "LLEU_03402":
        work_list.append([gene_name, gene_dic, gff_file, seq_dic, target_species])
    gene_list = pool.map(harvest_worker, work_list)
    pool.close()
    pool.join()
    neighbor_work = []
    for gene in gene_list:
        gene.set_neighbors(gene_list)
    print "gene data gotten"
    print str(datetime.datetime.now())
    print len(gene_list)
    for gene in gene_list:
        gene_objects[gene.name] = gene
    print "Dumping %s pickle" % target_species
#    pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle_test" % (target_species, target_species), 'wb')
    pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species), 'wb')
    pickle.dump(gene_objects, pickle_file)
    pickle_file.close()
    gene_objects = {}

def make_neighborhoods(attribute_list):
    gene = attribute_list[0]
    print "neighbors" + "\t" + gene.name
    gene_list = attribute_list[1]
    gene.set_neighbors(gene_list)

def harvest_worker(attribute_list):#gene_name, gene_dic, gff_file, seq_dic, target_species):
    gene_name = attribute_list[0]
    print gene_name
#    if gene_name not in ["LALB_00502"]:
#        return
    gene_dic = attribute_list[1]
    gff_file = attribute_list[2]
    seq_dic = attribute_list[3]
    target_species = attribute_list[4]
    if target_species == "LALB":
        vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_social_renamed.vcf.gz" % (target_species, target_species))
#        vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_solitary_renamed.vcf.gz" % (target_species, target_species))
    else:
        vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_filtered_miss.vcf.gz" % (target_species, target_species))
#    vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_testset.vcf.gz" % (target_species, target_species))
    cur_gene = Gene(gene_name, gene_dic[gene_name][0], gene_dic[gene_name][1])
    gene_deets = gff_file.getGene(gene_dic[gene_name][0], gene_name)[0]
    cur_gene.set_start(gene_deets["start"])
    cur_gene.set_end(gene_deets["end"])
    cur_gene.set_sequence(seq_dic[gene_dic[gene_name][0]], 5000)
    mrna_list = gff_file.getmRNA(gene_dic[gene_name][0], gene_name) #only one mRNA because working with longest iso
    for mrna_dic in mrna_list:
        cds_dic = gff_file.getCDS(gene_dic[gene_name][0], mrna_dic["Name"])
        cds_tuples = []
        for cds in cds_dic:
            cds_tuples.append((cds["start"], cds["end"]))
        sorted_tuples = sorted(cds_tuples, key = lambda tup: tup[0])
        cur_gene.add_cds(sorted_tuples)
    for mrna_dic in mrna_list:
        utr5_list = gff_file.getFivePrimeUTR(gene_dic[gene_name][0], mrna_dic["Name"])
        tuple_list = []
        for utr_dic in utr5_list:
            tuple_list.append((utr_dic["start"], utr_dic["end"]))
        cur_gene.add_utrs(tuple_list, "five")
        utr3_list = gff_file.getThreePrimeUTR(gene_dic[gene_name][0], mrna_dic["Name"])
        tuple_list = []
        for utr_dic in utr3_list:
            tuple_list.append((utr_dic["start"], utr_dic["end"]))
        cur_gene.add_utrs(tuple_list, "three")
    cur_gene.get_cds_sequence()
    cur_gene.get_flank_sequence(5000)
    cur_gene.get_utr_sequence()
    cur_gene.get_intron_sequence()
    cur_gene.get_genotypes(vcf_reader, 5000)
    return cur_gene

        
def cds_sequence_worker(gene):
    gene.get_cds_sequence()
    gene.get_flank_sequence(5000)
    gene.get_utr_sequence()
    gene.get_intron_sequence()
    return gene    

def get_gene_coords(gff_file, gene_dic, gene_objects, seq_dic):
    #Add coordinate data for all of the Gene objects
    for gene_name in gene_dic.keys():
        gene_deets = gff_file.getGene(gene_dic[gene_name][0], gene_name)[0]
        gene_objects[gene_name].set_start(gene_deets["start"])
        gene_objects[gene_name].set_end(gene_deets["end"])
        gene_objects[gene_name].set_sequence(seq_dic[gene_dic[gene_name][0]], 5000)
    return gene_objects

def get_utr_dic(gff_file, gene_dic, prime_end, gene_objects):
    #Add UTR coordinate data to all of the Gene objects
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
            gene_objects[gene_name].add_utrs(tuple_list, prime_end)
    return gene_objects

def get_cds_dic(gff_file, gene_dic, gene_objects):
    #Create Gene objects
    out_dic = {}
    for gene_name in gene_dic.keys():
        mrna_list = gff_file.getmRNA(gene_dic[gene_name][0], gene_name) #only one mRNA because working with longest iso
        for mrna_dic in mrna_list:
            cds_dic = gff_file.getCDS(gene_dic[gene_name][0], mrna_dic["Name"])
            cds_tuples = []
            for cds in cds_dic:
                cds_tuples.append((cds["start"], cds["end"]))
            sorted_tuples = sorted(cds_tuples, key = lambda tup: tup[0])
            gene_objects[gene_name].add_cds(sorted_tuples)
    return gene_objects


def check_og_complete(og_num, inspecies, outspecies, ortho_dic):
    #make sure that a particular OG has the required species
    #represented a single time. This is needed to check on 
    #neighboring genes before trying to use them
    if inspecies not in ortho_dic[og_num] or outspecies not in ortho_dic[og_num]:
        return False
    if ortho_dic[og_num][inspecies][0].count(inspecies) > 1 or ortho_dic[og_num][outspecies][0].count(outspecies) > 1:
        return False
    return True

def hka_test(inspecies, outspecies, seq_type, ortho_dic, out_path, num_threads, align_dir):
    #do all possible iterations of the HKA test for all genes between
    #two species.
    get_species_data(inspecies, num_threads)
#    sys.exit()
    get_species_data(outspecies, num_threads)
    if not os.path.isdir("%s/hka_tests/" % (out_path)):
        os.mkdir("%s/hka_tests/" % (out_path))
    og_map = open("%s/%s_%s_og_map.txt" % (out_path, inspecies, outspecies), 'w')
    hka_table = open("%s/hka_tests/%s_%s_%s_hka_direct.txt" % (out_path, inspecies, outspecies, seq_type), 'w')
    hka_table.write("Title: ingroup: %s, outgroup: %s\n" % (inspecies, outspecies))
    hka_line_list = []
    numloci = 0
    og_list = []
    for og_num in ortho_dic.keys():
        if inspecies not in ortho_dic[og_num] or outspecies not in ortho_dic[og_num]:
            continue
        if ortho_dic[og_num][inspecies][0].count(inspecies) > 1 or ortho_dic[og_num][outspecies][0].count(outspecies) > 1:
            continue
        if len(ortho_dic[og_num][inspecies]) == 1 and len(ortho_dic[og_num][outspecies]):
            og_list.append(og_num)
#    og_list = [2110]
#    og_list = [870]
#    og_list = og_list[0:50]
    work_list = []
    in_gene_dic = read_species_pickle(inspecies)
    print "%s pickle read" % inspecies
    out_gene_dic = read_species_pickle(outspecies)
    print "%s pickle read" % outspecies
    og_map_dic = make_og_gene_map(ortho_dic)
    for og in og_list:
        inspecies_gene = ortho_dic[og][inspecies][0][:-3]
        outspecies_gene = ortho_dic[og][outspecies][0][:-3]
        in_gene = in_gene_dic[inspecies_gene]
        out_gene = out_gene_dic[outspecies_gene]
        neighbor_dic = {}
        for neighbor_name in in_gene.neighbors:
            try:
                neighbor_og = og_map_dic[neighbor_name]
            except KeyError:
                print "%s not in orthology dictionary" % neighbor_name
                continue
            if not neighbor_og:
                continue
            if check_og_complete(neighbor_og, inspecies, outspecies, ortho_dic):
                outneighbor_name = ortho_dic[neighbor_og][outspecies][0][:-3]
                neighbor_dic[neighbor_og] = (in_gene_dic[neighbor_name], out_gene_dic[outneighbor_name])
        if len(neighbor_dic) < 5:
            continue
        og_map.write("OG_%s\t%s\t%s\n" % (og, inspecies_gene, outspecies_gene))
#        hka_line_list.append(hka_worker([inspecies, outspecies, seq_type, og, in_gene, out_gene, out_path, align_dir, neighbor_dic]))

        work_list.append([inspecies, outspecies, seq_type, og, in_gene, out_gene, out_path, align_dir, neighbor_dic])
    og_map.close()
    pool = multiprocessing.Pool(processes = num_threads)
    hka_lines = pool.map(hka_worker, work_list)
    pool.close()
    pool.join()
    full_length_lines = []
    jody_lines = []
    hkadirect_lines = []
    direct_pairwise_lines = []
    for line in hka_lines:
        if "too_short" not in line and "too_different" not in line and "HKA_failed" not in line[0] and "insufficient_sampling" not in line and "no_alignment" not in line:
            jody_lines.append(line[0])
    jody_correction(jody_lines, "%s/hka_tests/%s_%s_%s_hka_results_correct.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_results_correct_faster.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_results_correct_slower.txt" % (out_path, inspecies, outspecies, seq_type))
    for line in hka_lines:
        if "too_short" not in line and "too_different" not in line and "insufficient_sampling" not in line and "no_alignment" not in line:
            if len(line[1]) > 0:
                direct_pairwise_lines.append(line[2])
    pairwise_direct_correction(direct_pairwise_lines, "%s/hka_tests/%s_%s_%s_hka_direct_pairwise_results_correct.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_direct_pairwise_results_correct_faster.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_direct_pairwiseresults_correct_slower.txt" % (out_path, inspecies, outspecies, seq_type))

    for line in hka_lines:
        if "too_short" not in line and "too_different" not in line and "insufficient_sampling" not in line and "no_alignment" not in line:
            if len(line[1]) > 0:
                hkadirect_lines.append(line[1])

#this is all for running HKAdirect
    hka_table.write("nloci %s\n" % (len(hkadirect_lines)))
    hka_table.write("IDlocus\tnsam\tSegSites\tDivergence\tlength_pol\tlength_div\tfactor_chrn\n")
    print len(hkadirect_lines)
    for line in hkadirect_lines:
        hka_table.write(line)
    hka_table.close()
    cmd = ["/Genomics/kocherlab/berubin/local/src/HKAdirect/bin/HKAdirect", "%s/hka_tests/%s_%s_%s_hka_direct.txt" % (out_path, inspecies, outspecies, seq_type)]
    with open("%s/hka_tests/%s_%s_%s_hka_direct_results.txt" % (out_path, inspecies, outspecies, seq_type), 'w') as outfile:
        subprocess.call(cmd, stdout = outfile)
    outfile.close()
    hka_direct_correction("%s/hka_tests/%s_%s_%s_hka_direct_results.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_direct_results_correct.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_direct_results_correct_faster.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_direct_results_correct_slower.txt" % (out_path, inspecies, outspecies, seq_type))
    in_gene_dic = {}
    out_gene_dic = {}

def hka_direct_correction(infile, outfile, fasterfile, slowerfile):
    #read in HKAdirect output, correct the p-values and make
    #nice output
    reader = open(infile, 'rU')
    start_reading = False
    chi_list = []
    p_list = []
    slower_p_list = []
    faster_p_list = []
    for line in reader:
        if line.startswith("#IDloc"):
            start_reading = True
            continue
        if not start_reading:
            continue
        cur_line = line.split()
        if len(cur_line) == 0:
            break
        p = chisqprob(float(cur_line[-1]), 1)
        cur_line.append(p)
        chi_list.append(cur_line)
        p_list.append(p)
    pval_corr = smm.multipletests(p_list, alpha = 0.05, method = 'fdr_i')[1]
    for x in range(len(pval_corr)):
        chi_list[x].append(pval_corr[x])
    chi_list.sort(key = lambda x: x[-1])
    outfile = open(outfile, 'w')
    outfile.write("IDloci\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texpHKA_S\texpVar_S\texpHKA_div\texpVar_D\texpHKA_theta\tpartialHKA\tp_val\tbh5_p\n")
    for locus in chi_list:
        outfile.write("\t".join([str(i) for i in locus]) + "\n")
    outfile.close()
    faster = open(fasterfile, 'w')
    faster.write("IDloci\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texpHKA_S\texpVar_S\texpHKA_div\texpVar_D\texpHKA_theta\tpartialHKA\tp_val\tbh5_p\n")
    slower = open(slowerfile, 'w')
    slower.write("IDloci\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texpHKA_S\texpVar_S\texpHKA_div\texpVar_D\texpHKA_theta\tpartialHKA\tp_val\tbh5_p\n")
    for locus in chi_list:
        if float(locus[-1]) < 0.01:
            obs_s = float(locus[2])
            exp_s = float(locus[7])
            if obs_s > exp_s:
                faster.write("\t".join([str(i) for i in locus]) + "\n")
            else:
                slower.write("\t".join([str(i) for i in locus]) + "\n")
    faster.close()
    slower.close()

def pairwise_direct_correction(line_list, outfile, fasterfile, slowerfile):
    #take in alist of results from pairwise HKAdirect, correct the p-values
    #and make pretty output
    chi_list = []
    p_list = []
    for line in line_list:
        p = float(line[-1])
        if p < 0:
            continue
        chi_list.append(line)
        p_list.append(p)
    pval_corr = smm.multipletests(p_list, alpha = 0.05, method = 'fdr_i')[1]
    for x in range(len(pval_corr)):
        chi_list[x].append(pval_corr[x])
    chi_list.sort(key = lambda x: x[-1])
    outfile = open(outfile, 'w')
    outfile.write("OG\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texp_S\texpVar_S\texp_div\texpVar_D\texp_theta\tpartialHKA\tpval\tbh5_p\n")
    for locus in chi_list:
        outfile.write("\t".join([str(i) for i in locus]) + "\n")
    outfile.close()
    faster = open(fasterfile, 'w')
    faster.write("OG\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texp_S\texpVar_S\texp_div\texpVar_D\texp_theta\tpartialHKA\tpval\tbh5_p\n")
    slower = open(slowerfile, 'w')
    slower.write("OG\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texp_S\texpVar_S\texp_div\texpVar_D\texp_theta\tpartialHKA\tpval\tbh5_p\n")
    for locus in chi_list:
        if float(locus[-1]) < 0.01:
            obs_s = float(locus[2])
            exp_s = float(locus[7])
            if obs_s > exp_s:
                faster.write("\t".join([str(i) for i in locus]) + "\n")
            else:
                slower.write("\t".join([str(i) for i in locus]) + "\n")
    faster.close()
    slower.close()
    

def jody_correction(line_list, outfile, fasterfile, slowerfile):
    #take in a list of results from Hey's HKA and correct the 
    #p-values and write nice output files
    chi_list = []
    p_list = []
    for line in line_list:
        p = float(line.strip().split()[-1])
        chi_list.append(line.split())
        p_list.append(p)
    pval_corr = smm.multipletests(p_list, alpha = 0.05, method = 'fdr_i')[1]
    for x in range(len(pval_corr)):
        chi_list[x].append(pval_corr[x])
    chi_list.sort(key = lambda x: x[-1])
    outfile = open(outfile, 'w')
    outfile.write("OG\tobs_div\texp_div\tobs_poly\texp_poly\tp_val\tbh5_p\n")
    for locus in chi_list:
        outfile.write("\t".join([str(i) for i in locus]) + "\n")
    outfile.close()
    faster = open(fasterfile, 'w')
    faster.write("OG\tobs_div\texp_div\tobs_poly\texp_poly\tp_val\tbh5_p\n")
    slower = open(slowerfile, 'w')
    slower.write("OG\tobs_div\texp_div\tobs_poly\texp_poly\tp_val\tbh5_p\n")
    for locus in chi_list:
        if float(locus[-1]) < 0.01:
            obs_s = float(locus[1])
            exp_s = float(locus[2])
            if obs_s > exp_s:
                faster.write("\t".join([str(i) for i in locus]) + "\n")
            else:
                slower.write("\t".join([str(i) for i in locus]) + "\n")
    faster.close()
    slower.close()


def hka_worker(attribute_list):
    #workhorse for running HKA tests. takes a list of parameters 
    #to make it easier to parallelize. 
    inspecies = attribute_list[0]
    outspecies = attribute_list[1]
    seq_type = attribute_list[2]
    og_num = attribute_list[3]
    in_gene = attribute_list[4]
    out_gene = attribute_list[5]
    out_path = attribute_list[6]
    align_dir = attribute_list[7]
    neighbor_dic = attribute_list[8]                                 
    print og_num
    if not os.path.exists("%s/og_cds_%s.1.fas" % (align_dir, og_num)):
        return "OG_%s\tinsufficient_sampling\t" % og_num
    inpoly, outpoly, inseq, outseq, insample, outsample = gather_hka_data(in_gene, out_gene, seq_type)
    if len(inseq) < 100 or len(outseq) < 100:
        return "OG_%s\ttoo_short\n" % (og_num)
    if in_gene.strand == -1:
        inseq = str(Seq.Seq(inseq).reverse_complement())
    if out_gene.strand == -1:
        outseq = str(Seq.Seq(outseq).reverse_complement())
    if len(outseq) < len(inseq):
        inseq = inseq[0:len(outseq)]
    if len(inseq) < len(outseq):
        outseq = outseq[0:len(inseq)]
    pairwise_data = mafft_pairwise_diff_count(inseq, outseq, inspecies, outspecies, og_num, in_gene, out_gene, out_path, seq_type)
    if pairwise_data == "no_alignment":
        return "OG_%s\ttoo_different\n" % (og_num)
    average_diff = pairwise_data[0]
    align_len = pairwise_data[1]
    inpoly = pairwise_data[2]
    outpoly = pairwise_data[3]
    inlen = pairwise_data[4]
    outlen = pairwise_data[5]
    if align_len < 500 or inlen < 500:
        return "OG_%s\ttoo_different\n" % (og_num)
    if 1.0 * average_diff / align_len > 0.10:
        return "OG_%s\ttoo_different\t%s\n" % (og_num, 1.0 * average_diff / align_len)
    if 1.0 * align_len / inlen < 0.90:
        return "OG_%s\ttoo_different\t%s\n" % (og_num, 1.0 * align_len / inlen)
    if insample < 5 or outsample < 5:
        return "OG_%s\tinsufficient_sampling\t" % og_num
    #this is the line for Jody Hey's HKA
#    return "OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, 1.0, len(inseq), len(outseq), align_len, insample*2, outsample*2, inpoly, outpoly, average_diff)
    if not os.path.exists("%s/hka_tests/hkadirect_pairwise_%s_%s" % (out_path, inspecies, outspecies)):
        os.mkdir("%s/hka_tests/hkadirect_pairwise_%s_%s" % (out_path, inspecies, outspecies))
        os.mkdir("%s/hka_tests/hey_%s_%s" % (out_path, inspecies, outspecies))
    directfile = open("%s/hka_tests/hkadirect_pairwise_%s_%s/OG_%s_%s_%s_%s_hka_direct_pairwise.txt" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies, seq_type), 'w')
    directfile.write("Title: ingroup: %s, outgroup: %s\n" % (inspecies, outspecies))
    directfile.write("nloci 2\n")# % (len(hkadirect_lines)))
    directfile.write("IDlocus\tnsam\tSegSites\tDivergence\tlength_pol\tlength_div\tfactor_chrn\n")

    outfile = open("%s/hka_tests/hey_%s_%s/OG_%s_%s_%s_%s_hka.txt" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies, seq_type), 'w')
    neighbor_stats_list = [] 
    direct_stats_list = []
    for neighbor_og, neighbor in neighbor_dic.items():
        neighbor_stats = fourfold_degenerate_string([neighbor[0], neighbor[1], align_dir, neighbor_og, inspecies, outspecies])

        if "insufficient_sampling" not in neighbor_stats:
            neighbor_stats_list.append(neighbor_stats[1])
            direct_stats_list.append(neighbor_stats[0])
    if len(direct_stats_list) < 2:
        return "OG_%s\tinsufficient_sampling\t" % og_num
    direct_background = sum_fourf_list(direct_stats_list)
    if inspecies == "LALB":
        insample = insample
    else:
        insample = insample * 2
    if outspecies == "LALB":
        outsample = outsample
    else:
        outsample = outsample * 2
    outfile.write("HKA table\n%s\n%s %s\n" % (len(neighbor_stats_list) + 1, inspecies, outspecies))
    outfile.write("%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, seq_type[0:2], 1.0, inlen, outlen, align_len, int(insample), int(outsample), inpoly, outpoly, average_diff))
    for neighbor in neighbor_stats_list:
        outfile.write(neighbor)
    outfile.close()
    switch_dir = os.getcwd()
    os.chdir("%s/hka_tests/hey_%s_%s" % (out_path, inspecies, outspecies))
    cmd = ["/Genomics/kocherlab/berubin/local/src/HKA/hka", "-DOG_%s_%s_%s_%s_hka.txt" % (og_num, inspecies, outspecies, seq_type), "-R%s/hka_tests/hey_%s_%s/OG_%s_%s_%s_%s_hka_results" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies, seq_type), "-S100"]
    FNULL = open(os.devnull, 'w')
    subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
    os.chdir(switch_dir)
    p_val, obs_poly, exp_poly = parse_hka_result(inspecies, outspecies, seq_type, "%s/hka_tests/hey_%s_%s/OG_%s_%s_%s_%s_hka_results.hka" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies, seq_type))
    directfile.write("OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, insample, inpoly, average_diff, inlen, align_len, 1))
    directfile.write(direct_background)
    directfile.close()
    cmd = ["/Genomics/kocherlab/berubin/local/src/HKAdirect/bin/HKAdirect", "%s/hka_tests/hkadirect_pairwise_%s_%s/OG_%s_%s_%s_%s_hka_direct_pairwise.txt" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies, seq_type)]
    with open("%s/hka_tests/hkadirect_pairwise_%s_%s/OG_%s_%s_%s_%s_hka_direct_pairwise_results.txt" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies, seq_type), 'w') as outfile:
        subprocess.call(cmd, stdout = outfile)
    outfile.close()
    directp = get_direct_p("%s/hka_tests/hkadirect_pairwise_%s_%s/OG_%s_%s_%s_%s_hka_direct_pairwise_results.txt" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies, seq_type))
    if p_val == "HKA_failed":
        return "HKA_failed","OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, insample, inpoly, average_diff, inlen, align_len, 1), directp
    #the first string here is the result from this individual run of Hey's HKA
    #the second string is the input line for HKAdirect
    #the third is the output of pairwise hkadirect
    return "OG_%s\t%s\t%s\t%s\n" % (og_num, obs_poly, exp_poly, p_val), "OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, insample, inpoly, average_diff, inlen, align_len, 1), directp

def get_direct_p(direct_resultsfile):
    #gets the output stats from an HKAdirect output file
    reader = open(direct_resultsfile, 'rU')
    for line in reader:
        if line.startswith("OG_"):
            data_line = line.split()
        if line.startswith("Significance"):
            pval = float(line.split("=")[-1].strip())
            data_line.append(pval)
            break
    return data_line

def fourfold_degenerate_string(attribute_list):
    #returns two tuples of data about the fourfold degenerate
    #sites in a particular gene alignment. the first two is a set of
    #data that will be used for HKAdirect. The second tuple is just
    #a single string formatted for use by Hey's HKA
    ingene = attribute_list[0]
    outgene = attribute_list[1] 
    align_dir = attribute_list[2]
    og = attribute_list[3] 
    inspecies = attribute_list[4]
    outspecies = attribute_list[5]
    if inspecies == "LALB":
        insample = ingene.average_n
    else:
        insample = ingene.average_n * 2
    if outspecies == "LALB":
        outsample = outgene.average_n
    else:
        outsample = outgene.average_n * 2
    if not os.path.exists("%s/og_cds_%s.1.fas" % (align_dir, og)):
        return "OG_%s\tinsufficient_sampling\t" % og
    if insample < 5:
        return "insufficient_sampling"
    inseq, outseq = get_prank_aligned_seqs(align_dir, og, inspecies, outspecies)
    fix_fourfold = count_fourfold(inseq, outseq, ingene, outgene)
    in_poly_fourf = ingene.fourfold
    out_poly_fourf = outgene.fourfold
    in_potent_fourfold = ingene.potent_fourfold
    out_potent_fourfold = outgene.potent_fourfold
    align_potent_fourf = potential_aligned_fourfold_sites(inseq, outseq)
    
    return (fix_fourfold, in_poly_fourf, out_poly_fourf, in_potent_fourfold, out_potent_fourfold, align_potent_fourf, int(insample), int(outsample)), ("%s_4d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og, 1.0, in_potent_fourfold, out_potent_fourfold, align_potent_fourf, int(insample), int(outsample), in_poly_fourf, out_poly_fourf, fix_fourfold))
    
def sum_fourf_list(fourf_list):
    #sums all of the variables for fourfold degenerate sites
    #across a bunch of genes and makes them into a string
    #for HKAdirect
    fix_fourfold = in_poly_fourf = out_poly_fourf = in_potent_fourfold = out_potent_fourfold = align_potent_fourf = insample = outsample = 0
    x = 0
    num_loci = 0.0
    for fourf_stats in fourf_list:
        if "insufficient" in fourf_stats:
            continue
        fix_fourfold += fourf_stats[0]
        in_poly_fourf += fourf_stats[1]
        out_poly_fourf += fourf_stats[2]
        in_potent_fourfold += fourf_stats[3]
        out_potent_fourfold += fourf_stats[4]
        align_potent_fourf += fourf_stats[5]
        insample += fourf_stats[6]
        outsample += fourf_stats[7]
        num_loci += 1
    return "4d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (round(insample / num_loci), in_poly_fourf, fix_fourfold, in_potent_fourfold, align_potent_fourf, 1)

def parse_hka_result(inspecies, outspecies, seq_type, hkafile):
    #this is for reading the result files from Jody Hey's HKA
    reader = open(hkafile, 'rU')
    poly_data = True
    fixed_data = False
    line = reader.readline()
    p_val = -1
    chi_val = -1
    flank_div = []
    in_flank_poly = []
    while "program by Jody Hey" not in line:
        line = reader.readline()
        if "LOCI AND LENGTHS" in line:
            line = reader.readline()
            cur_line = line.split()
            if int(cur_line[2]) == 0 and int(cur_line[3]) == 0 and int(cur_line[4]) == 0:
                return "HKA_failed", "NA", "NA", "NA", "NA"
            if "program by Jody Hey" in line:
                return "HKA_failed", "NA", "NA", "NA", "NA"
        if line.startswith("*** LOCUS") and poly_data:
            line = reader.readline()
            line = reader.readline()
            in_flank_poly = line.split()[2:4]
            line = reader.readline()
            out_flank_poly = line.split()[2:4]
            line = reader.readline()
            line = reader.readline()
            line = reader.readline()
            line = reader.readline()
            in_syn_poly = line.split()[2:4]
            line = reader.readline()
            out_syn_poly = line.split()[2:4]
            poly_data = False
            fixed_data = True
        #THIS NEXT BIT IS ONLY FOR TWO LOCI
        if line.startswith("*** LOCUS") and fixed_data:
            line = reader.readline()
            line = reader.readline()
            flank_div = line.split()[0:2]
            reader.readline()
            reader.readline()
            syn_div = reader.readline().split()[0:2]
            fixed_data = False
        if "SUM OF DEVIATIONS:" in line:
            chi_val = float(line.strip().split(" ")[-1])
        if line.startswith("Degrees of Freedom"):
            if "Degrees of Freedom:  " in line:
                dof = int(line.strip().split(" ")[4])
                break
            else:
                dof = int(line.strip().split(" ")[3])
                break
    if chi_val == -1:
        return "HKA_failed", "NA", "NA"
    p_val = chisqprob(chi_val, dof)
    return p_val, float(in_flank_poly[0]), float(in_flank_poly[1])
        
def get_prank_aligned_seqs(align_dir, og_num, inspecies, outspecies):
    #get the aligned sequences for two species in a particular OG
    reader = SeqIO.parse("%s/og_cds_%s.1.fas" % (align_dir, og_num), format = 'fasta')
    for rec in reader:
        if rec.id[0:4] == inspecies:
            inseq = str(rec.seq)
            inseq_name = rec.id[:-3]
        elif rec.id[0:4] == outspecies:
            outseq = str(rec.seq)
            outseq_name = rec.id[:-3]
    return inseq, outseq

def align_len(seq1, seq2):
    missing = ["N", "-"]
    total_len = 0
    for x in range(len(seq1)):
        if seq1[x] not in missing and seq2[x] not in missing:
            total_len += 1
    return total_len

def potential_aligned_sites(inseq, outseq):
    #the number of synonymous and nonsynonymous sites in two aligned
    #sequences. remember that these are counted in kind of a 
    #complicated way because some sites have potential to be both 
    #synonymous and nonsynonymous
    x = 0
    syn_sites = 0
    nsyn_sites = 0
    outsyn_sites = 0
    potent_dic = changes.potent_dic()
    while x < len(inseq):
        ambig_codon = False
        incodon = inseq[x:x+3].upper()
        outcodon = outseq[x:x+3].upper()
        for ambig in ["N", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "-"]:
            if ambig in incodon or ambig in outcodon:
                ambig_codon = True
        if ambig_codon:
            x = x + 3
            ambig_codon = False
            continue
        nsyn_sites += potent_dic["N"][incodon]
        syn_sites += potent_dic["S"][incodon]
        outsyn_sites += potent_dic["S"][outcodon]
        x = x + 3
    if outsyn_sites <= syn_sites:
        syn_sites = outsyn_sites
    return syn_sites, nsyn_sites

def potential_aligned_fourfold_sites(inseq, outseq):
    #the number of sites that are fourfold degenerate in two aligned
    #sequences.
    x = 0
    fourf_sites = 0
    fourfold_list = changes.fourfold_codons()
    incount = 0
    while x < len(inseq):
        ambig_codon = False
        incodon = str(inseq[x:x+3].upper())
        outcodon = str(outseq[x:x+3].upper())
        for ambig in ["N", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "-"]:
            if ambig in incodon or ambig in outcodon:
                ambig_codon = True
        if ambig_codon:
            x = x + 3
            ambig_codon = False
            continue
        if incodon in fourfold_list:
            incount += 1
            if str(incodon)[0] == str(outcodon)[0] and str(incodon)[1] == str(outcodon)[1]:
                fourf_sites += 1
        x = x + 3
    return fourf_sites


def gather_hka_data(inspecies_gene, outspecies_gene, seq_type):
    #Gather together all of the necessary HKA test data from different
    #types of sequence.
    if seq_type == "flank":
        inpoly = inspecies_gene.flank_subs
        outpoly = outspecies_gene.flank_subs
        #again, trying to make it faster and less memory
#        inseq = "".join(in_gene_dic[inspecies_gene].flank_dic.values())
#        outseq = "".join(out_gene_dic[outspecies_gene].flank_dic.values())
        inseq = inspecies_gene.flank_seq
        outseq = outspecies_gene.flank_seq       
        insample = inspecies_gene.average_n
        outsample = outspecies_gene.average_n
    elif seq_type == "intron":
        inpoly = inspecies_gene.intron_subs
        outpoly = outspecies_gene.intron_subs
        inseq = inspecies_gene.intron_seq
        outseq = outspecies_gene.intron_seq
        insample = inspecies_gene.average_n
        outsample = outspecies_gene.average_n
    elif seq_type == "three_prime_utr":
        inpoly = inspecies_gene.utr3_subs
        outpoly = outspecies_gene.utr3_subs
        inseq = "".join(inspecies_gene.utr3_dic.values())
        outseq = "".join(outspecies_gene.utr3_dic.values())
        insample = inspecies_gene.average_n
        outsample = outspecies_gene.average_n
    elif seq_type == "five_prime_utr":
        inpoly = inspecies_gene.utr5_subs
        outpoly = outspecies_gene.utr5_subs
        inseq = "".join(inspecies_gene.utr5_dic.values())
        outseq = "".join(outspecies_gene.utr5_dic.values())
        insample = inspecies_gene.average_n
        outsample = outspecies_gene.average_n

    return inpoly, outpoly, inseq, outseq, insample, outsample    

def muscle_pairwise_diff_count(seq1, seq2, inspecies, outspecies, og_num):
    #Align two sequences using muscle and return the number of 
    #differences between them. Not up to date. Use MAFFT instead.
    handle = StringIO()
    rec1 = SeqRecord(Seq.Seq(seq1), id = "inseq")
    rec2 = SeqRecord(Seq.Seq(seq2), id = "outseq")
    recs = [rec1, rec2]
    SeqIO.write(recs, handle, "fasta")
    data = handle.getvalue()
    muscle_cline = MuscleCommandline()
    stdout, stderr = muscle_cline(stdin = data)
    align = AlignIO.read(StringIO(stdout), "fasta")
    align_dic = {}
    outfile = open("/Genomics/kocherlab/berubin/annotation/orthology/muscle_files/OG_%s_%s_%s.afa" % (og_num, inspecies, outspecies), 'w')
    for rec in align:
        align_dic[rec.id] = str(rec.seq)
        outfile.write(">%s\n%s\n" % (rec.id, str(rec.seq)))
    outfile.close()
    counter = 0
    missing_data = ["N", "-"]
    indel_count = 0
    for x in range(len(align_dic["inseq"])):
        if align_dic["inseq"][x] in missing_data or align_dic["outseq"][x] in missing_data:
            indel_count += 1
            continue
        if align_dic["inseq"][x] != align_dic["outseq"][x]:
            counter += 1
    return counter, len(align_dic["inseq"]) - indel_count

def mafft_pairwise_diff_count(seq1, seq2, inspecies, outspecies, og_num, in_gene, out_gene, out_path, seq_type):
    #Align two sequences using mafft and return the number of 
    #differences between them. Since the number of polymorphisms
    #counted in these regions is dependent on where the alignable
    #regions are, this also provides the number of polymorphisms
    if not os.path.exists("%s/conservation_files" % (out_path)):
        os.mkdir("%s/conservation_files" % (out_path))

    if not os.path.exists("%s/conservation_files/%s_%s/" % (out_path, inspecies, outspecies)):
        os.mkdir("%s/conservation_files/%s_%s" % (out_path, inspecies, outspecies))
    intuple, outtuple = conserved_noncoding(seq1, seq2, og_num, inspecies, outspecies, "%s/conservation_files/%s_%s" % (out_path, inspecies, outspecies))
    if intuple[1] - intuple[0] < 500 or outtuple[1] - outtuple[0] < 500:
        return "no_alignment"
    if abs((intuple[1] - intuple[0]) - (outtuple[1] - outtuple[0])) > 500:
        return "no_alignment"
    if not os.path.exists("%s/mafft_files" % (out_path)):
        os.mkdir("%s/mafft_files" % (out_path))

    if not os.path.exists("%s/mafft_files/%s_%s" % (out_path, inspecies, outspecies)):
        os.mkdir("%s/mafft_files/%s_%s" % (out_path, inspecies, outspecies))
    seqfile = open("%s/mafft_files/%s_%s/OG_%s_%s_%s.fa" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies), 'w')
    seqfile.write(">inseq\n%s\n" % (seq1[intuple[0]:intuple[1]]))
    seqfile.write(">outseq\n%s\n" % (seq2[outtuple[0]:outtuple[1]]))

    seqfile.close()
    FNULL = open(os.devnull, 'w')
    with open("%s/mafft_files/%s_%s/OG_%s_%s_%s.afa" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies), 'w') as outfile:
        subprocess.call(["linsi", "%s/mafft_files/%s_%s/OG_%s_%s_%s.fa" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies)], stdout = outfile, stderr=FNULL)

    outfile.close()
    align_dic = {}
    reader = SeqIO.parse("%s/mafft_files/%s_%s/OG_%s_%s_%s.afa" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies), format = 'fasta')
    for rec in reader:
        align_dic[rec.id] = str(rec.seq)
    average_diff = 0
    missing_data = ["N", "-", "n"]
    indel_count = 0
    if seq_type == "flank":
        if in_gene.strand == -1:
            in_gene_alt_dic = in_gene.check_polymorphisms_fixed(align_dic, in_gene.end, "inseq")
        else:
            in_gene_alt_dic = in_gene.check_polymorphisms_fixed(align_dic, in_gene.flank_end - intuple[1], "inseq")
        if out_gene.strand == -1:
            out_gene_alt_dic = out_gene.check_polymorphisms_fixed(align_dic, out_gene.flank_end - outtuple[1], "outseq")
        else:
            out_gene_alt_dic = out_gene.check_polymorphisms_fixed(align_dic, out_gene.flank_start + outtuple[0], "outseq")
    elif seq_type == "intron":
        in_gene_alt_dic = in_gene.coding_fixed_align(align_dic["inseq"])
        out_gene_alt_dic = out_gene.coding_fixed_align(align_dic["outseq"])

    for x in range(len(align_dic["inseq"])):
        if align_dic["inseq"][x] in missing_data or align_dic["outseq"][x] in missing_data:
            indel_count += 1
            continue
        insite_list = [align_dic["inseq"][x].lower()]

        outsite_list = [align_dic["outseq"][x].lower()]
        if x in in_gene_alt_dic.keys():
            insite_list.append(str(in_gene_alt_dic[x]).lower())
        if x in out_gene_alt_dic.keys():
            outsite_list.append(str(out_gene_alt_dic[x]).lower())
        overlap = False
        for nuc in insite_list:
            if nuc in outsite_list:
                overlap = True
        if not overlap:
            average_diff += 1
    return average_diff, len(align_dic["inseq"]) - indel_count, len(in_gene_alt_dic), len(out_gene_alt_dic), intuple[1] - intuple[0], outtuple[1] - outtuple[0]

def mk_test(inspecies, outspecies, ortho_dic, align_dir, out_path, num_threads):
    get_species_data(inspecies, num_threads)
#    sys.exit()
    get_species_data(outspecies, num_threads)
    in_gene_dic = read_species_pickle(inspecies)
    print "%s pickle read" % inspecies
    out_gene_dic = read_species_pickle(outspecies)
    print "%s pickle read" % outspecies
    if not os.path.exists("%s/mk_tests/" % (out_path)):
        os.mkdir("%s/mk_tests/" % (out_path))
    mk_table = open("%s/mk_tests/%s_%s_mk_snipre.txt" % (out_path, inspecies, outspecies), 'w')
    mk_table.write("geneID\tPR\tFR\tPS\tFS\tTsil\tTrepl\tnout\tnpop\n")
    
    for og_num in ortho_dic.keys():
#        print "OG %s" % og_num
#        if og_num != 7:
#            continue
        if inspecies not in ortho_dic[og_num] or outspecies not in ortho_dic[og_num]:
            continue
        if ortho_dic[og_num][inspecies][0].count(inspecies) > 1 or ortho_dic[og_num][outspecies][0].count(outspecies) > 1:
            continue
        if len(ortho_dic[og_num][inspecies]) == 1 and len(ortho_dic[og_num][outspecies]) == 1:

            inspecies_gene = ortho_dic[og_num][inspecies][0][:-3]
            outspecies_gene = ortho_dic[og_num][outspecies][0][:-3]
            inseq, outseq = get_prank_aligned_seqs(align_dir, og_num, inspecies, outspecies)
            if inseq.count("-") / (len(inseq) * 1.0) > 0.05 or outseq.count("-") / (len(outseq) * 1.0) > 0.05:
                continue

#            if og_num != 2110:
#                continue
#            fix_syn, fix_nsyn = count_sub_types(inseq, outseq)
            in_gene = in_gene_dic[inspecies_gene]
            out_gene = out_gene_dic[outspecies_gene]

            fix_syn, fix_nsyn = fixed_sub_types(inspecies, outspecies, og_num, align_dir, in_gene, out_gene)
            in_poly_syn = in_gene.syn_count #in_gene_dic[inseq_name].syn_count
            in_poly_nsyn = in_gene.nsyn_count #in_gene_dic[inseq_name].nsyn_count
            in_potent_syn = in_gene.potent_syn #in_gene_dic[inseq_name].potent_syn
            in_potent_nsyn = in_gene.potent_nsyn #in_gene_dic[inseq_name].potent_nsyn
            in_n_size = in_gene.average_n #in_gene_dic[inseq_name].average_n * 2
            if inspecies == "LALB":
                insample = in_n_size
            else:
                insample = in_n_size * 2
#            if (fix_syn + fix_nsyn) / (in_potent_syn + in_poten_nsyn) > 0.05:
#                continue
            if insample == 0:
                continue
            mk_table.write("\t".join([str(og_num), str(in_poly_nsyn), str(fix_nsyn), str(in_poly_syn), str(fix_syn), str(in_potent_syn), str(in_potent_nsyn), "1", str(insample)]) + "\n")
    mk_table.close()

def fixed_sub_types(inspecies, outspecies, og_num, align_dir, in_gene, out_gene):
    inseq, outseq = get_prank_aligned_seqs(align_dir, og_num, inspecies, outspecies)
    in_gene_alt_dic = in_gene.coding_fixed_align(inseq)
    out_gene_alt_dic = out_gene.coding_fixed_align(outseq)

    nsyn_count = 0
    syn_count = 0
    same_codons = 0
    x = 0
    while x < len(inseq):
        #ignore codons with more than one divergent site (Bin's thing)
        diff_count = 0
        missing_data = False
        for cod_index in range(3):
            if inseq[x+cod_index] in ["N", "-"] or outseq[x+cod_index] in ["N", "-"]:
                x = x + 3
                missing_data = True
                break
        if missing_data:
            continue

        incodon_list = []
        outcodon_list = []
        for cod_index in range(3):
            incodon_list.append([inseq[x+cod_index].upper()])
            outcodon_list.append([outseq[x+cod_index].upper()])
            if x+cod_index in in_gene_alt_dic.keys():
                incodon_list[cod_index].append(str(in_gene_alt_dic[x+cod_index]).upper())
            if x+cod_index in out_gene_alt_dic.keys():
                outcodon_list[cod_index].append(str(out_gene_alt_dic[x+cod_index]).upper())
#        print incodon_list
#        print x
#        print outcodon_list
        overlap = 0
        overlap_incodon = list(inseq[x:x+3])
        overlap_outcodon = list(outseq[x:x+3])
        for cur_site in range(len(incodon_list)):
            for nuc in incodon_list[cur_site]:
                if nuc in outcodon_list[cur_site]:
                    overlap_incodon[cur_site] = nuc
                    overlap_outcodon[cur_site] = nuc
                    #if there is overlap at a site then just got to next site
                    #don't need to keep checking current site
                    break
        incodon = "".join(overlap_incodon)
        outcodon = "".join(overlap_outcodon)
#        print incodon
#        print outcodon
        if incodon == outcodon:
            same_codons += 1
        else:
            num_diffs = 0
            for cod_index in range(3):
                if incodon[cod_index] != outcodon[cod_index]:
                    num_diffs += 1
            if num_diffs == 1:
                inaa = str(Seq.Seq(incodon).translate())
                outaa = str(Seq.Seq(outcodon).translate())
                if inaa == outaa:
                    syn_count += 1
                else:
                    nsyn_count += 1
#        print syn_count
#        print nsyn_count
        x = x + 3
#    print in_gene_alt_dic
#    print out_gene_alt_dic

    return syn_count, nsyn_count        


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
    SOCIAL = ["HLIG","LMAL", "LMAR", "LPAU", "LZEP", "AAUR", "LCAL", "LALB", "LZEP", "HRUB"]
    REV_SOLITARY = ["APUR", "LLEU", "LOEN", "LVIE", "LFIG"]
    if foreground == "social":
        fore_list = SOCIAL
    if foreground == "solitary":
        fore_list = REV_SOLITARY
    if foreground == "model_d":
        tree_prep = False
    if foreground == "ancestral":
        reader = SeqIO.parse("%s/og_cds_%s.1.fas" % ( indir, orthogroup), format = 'fasta')
    else:
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
        elif test_type == "ancestral":
            cur_out_dir = "%s/OG_%s" % (outdir, cur_og)
#            if not os.path.exists("%s/OG_%s" % (outdir, cur_og)):
#                os.mkdir("%s/OG_%s" % (outdir, cur_og))
#            prep_paml_files(cur_og, indir, "%s/OG_%s" % (outdir, cur_og), "ancestral", phylogeny_file)
            prep_paml_files(cur_og, indir, outdir, "ancestral", phylogeny_file)
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
        elif test_type == "ancestral":
            worker = Worker(work_queue, result_queue, paml_tests.ancestor_reconstruction)          
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

def count_fourfold(seq1, seq2, in_gene, out_gene):
    #counts the number of differences at fourfold degenerate sites
    #in two sequences. Takes into account the fourfold degenerate
    #polymorphisms in both sequences so as not to overcount divergence
    fourfold_list = changes.fourfold_codons()
    empty_chars = ["N", "-", "X"]
    x = 0
    fourfold_count = 0
    in_gene_alt_dic = in_gene.coding_fixed_align(seq1)
    out_gene_alt_dic = out_gene.coding_fixed_align(seq2)
    while x < len(seq1):
        if seq1[x:x+3] in fourfold_list:
            if seq1[x:x+2] == seq2[x:x+2]:
                if seq1[x+2] != seq2[x+2]:
                    in_alt_list = [seq1[x+2]]
                    out_alt_list = [seq2[x+2]]
                    if x+2 in in_gene_alt_dic.keys():
                        in_alt_list.append(in_gene_alt_dic[x+2])
                    if x+2 in out_gene_alt_dic.keys():
                        out_alt_list.append(out_gene_alt_dic[x+2])
                    overlap = False
                    for nuc in in_alt_list:
                        if nuc in out_alt_list:
                            overlap = True
                    if not overlap:
                        fourfold_count += 1
        x = x + 3
    return fourfold_count

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
    for species in ["APUR", "HLIG", "LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LFIG", "AAUR", "AVIR", "HRUB"]:
        reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (species, species), format = 'fasta')
        seq_dic[species] = {}
        for rec in reader:
            seq_dic[species][rec.id] = str(rec.seq)
    species = "LALB"
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/LALB/LALB_v3/LALB/LALB_OGS_v1.0_longest_isoform.cds.fasta", format = 'fasta')
    seq_dic["LALB"] = {}
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
    SPECIES_LIST = ["APUR", "HLIG", "LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LFIG", "AAUR", "AVIR", "HRUB", "LALB", "Nmel", "Dnov"]
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

            
def conserved_noncoding(seq1, seq2, og_num, inspecies, outspecies, working_dir):
    #this method is for identifying the alignable parts of noncoding
    #sequence
    outfile = open("%s/OG_%s_%s_%s_in.fa" % (working_dir, og_num, inspecies, outspecies), 'w')
    outfile.write(">%s\n%s\n" % (inspecies, seq1))
    outfile.close()
    query_seq = open("%s/OG_%s_%s_%s_out.fa" % (working_dir, og_num, inspecies, outspecies), 'w')
    query_seq.write(">%s\n%s\n" % (outspecies, seq2))
    query_seq.close()
    cmd = ["/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/makeblastdb", "-in", "%s/OG_%s_%s_%s_in.fa" % (working_dir, og_num, inspecies, outspecies), "-out", "%s/OG_%s_%s_%s_in_db" % (working_dir, og_num, inspecies, outspecies), "-dbtype", "nucl"]
    FNULL = open(os.devnull, 'w')
    subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
    cmd = ["/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/blastn", "-query", "%s/OG_%s_%s_%s_out.fa" % (working_dir, og_num, inspecies, outspecies), "-db", "%s/OG_%s_%s_%s_in_db" % (working_dir, og_num, inspecies, outspecies), "-outfmt", "6", "-out", "%s/OG_%s_%s_%s.txt" % (working_dir, og_num, inspecies, outspecies)]
    subprocess.call(cmd)
    longest_intuple, longest_outtuple = parse_conserved_noncoding("%s/OG_%s_%s_%s.txt" % (working_dir, og_num, inspecies, outspecies))
    return longest_intuple, longest_outtuple

def parse_conserved_noncoding(blastfile):
    if os.stat(blastfile).st_size == 0:
        #if there are no hits then there is nothing alignable
        return (0, 0), (0, 0)
    reader = open(blastfile, 'rU')
    coord_pairs = {}
    for line in reader:
        cur_line = line.split()
        if float(cur_line[2]) < 90:
            continue
        out_tuple = (int(cur_line[6]), int(cur_line[7]))
        in_tuple = (int(cur_line[8]), int(cur_line[9]))
        coord_pairs[in_tuple] = out_tuple
    if len(coord_pairs.keys()) == 0:
        return (0, 0), (0, 0)
    in_tuple_list = coord_pairs.keys()
    for in_coords in coord_pairs.keys():
        in_tuple_list = overlap(in_coords, in_tuple_list)
    out_tuple_list = coord_pairs.values()
    for out_coords in coord_pairs.values():
        out_tuple_list = overlap(out_coords, out_tuple_list)
    longest_intuple = in_tuple_list[0]
    for in_tuple in in_tuple_list:
        if in_tuple[1]-in_tuple[0] > longest_intuple[1] - longest_intuple[0]:
            longest_tuple = in_tuple
    longest_outtuple = out_tuple_list[0]
    for out_tuple in out_tuple_list:
        if out_tuple[1]-out_tuple[0] > longest_outtuple[1] - longest_outtuple[0]:
            longest_outtuple = out_tuple
    return longest_intuple, longest_outtuple

         
def overlap(mytuple, tuplelist):
    included = False
    i = 0
    while i < len(tuplelist):
        cur_tuple = tuplelist[i]
        if mytuple[0] <= cur_tuple[0] and mytuple[1] >= cur_tuple[0]-50:
            if mytuple[1] >= cur_tuple[1]:
                new_tuple = mytuple
                mytuple = new_tuple
                del tuplelist[i]
                i = -1
            else:
                new_tuple = (mytuple[0], cur_tuple[1])
                mytuple = new_tuple
                del tuplelist[i]
                i = -1
        elif mytuple[1] >= cur_tuple[1] and mytuple[0] <= cur_tuple[1] + 50:
            if mytuple[0] >= cur_tuple[0]:
                new_tuple = (cur_tuple[0], mytuple[1])
                mytuple = new_tuple
                del tuplelist[i]
                i = -1   
        elif mytuple[0] >= cur_tuple[0] and mytuple[1] <= cur_tuple[1]:
            included = True
            break
        i += 1
    if not included:
        tuplelist.append(mytuple)
    return tuplelist
   
def read_ancestral_seq(infile):
    reader = open(infile, 'rU')
    tree_line = False
    ancestral_seqs = False
    node_index_dic = {}
    node_seqs = {}
    for line in reader:
        if line.startswith("tree with node labels for Rod Page's TreeView"):
            tree_line = True
            continue
        if tree_line:
            tree_string = line.strip().replace(" ", "")
            tree = PhyloTree(tree_string, format = 8)
            for node in tree.traverse():
                if not node.is_leaf():
                    children = []
                    for child in node.traverse():
                        if child.is_leaf():
                            children.append(child.name.split("_")[1])
                    node_index_dic[node.name] = children
            tree_line = False
        if line.startswith("List of extant and reconstructed sequences"):
            ancestral_seqs = True
            continue
        if ancestral_seqs:
            if line.startswith("node #"):
 #               print line
                node_name = line.split()[1][1:]
                cur_seq = "".join(line.split()[2:])
                node_seqs[node_name] = cur_seq
#                print node_name
 #               print cur_seq
        if line.startswith("Overall accuracy of the"):
            break
    return node_index_dic, node_seqs

def compile_ancestrals(og_list, indir, aligndir, outdir):
    
    species_syn_dists = {}
    species_nsyn_dists = {}
    species_fourf_dists = {}
    for og_num in og_list:
        print og_num
        node_index_dic, node_seqs = read_ancestral_seq("%s/og_%s_working/rst" % (indir, og_num))
        inseq_dic = {}
        reader = SeqIO.parse("%s/og_cds_%s.1.fas" % (aligndir, og_num), format = 'fasta')
        for rec in reader:
            inseq_dic[rec.id[0:-5]] = str(rec.seq)
        if len(inseq_dic) < 4:
            continue
        nearest_anc = get_nearest_anc(node_index_dic, inseq_dic)
        
#        print inseq_dic.keys()
#        print node_seqs
        for species, seq in inseq_dic.items():
            syn_pos, nsyn_pos, fourf_pos = sub_positions(seq, node_seqs[nearest_anc[species]])
#            syn_pos, nsyn_pos, fourf_pos = sub_positions(seq, inseq_dic["ECOL"])
            if species not in species_syn_dists.keys():
                species_syn_dists[species] = []
                species_nsyn_dists[species] = []
                species_fourf_dists[species] = []
                
            syn_dists = calc_distances(syn_pos)
#            print "syn dists"
#            print syn_pos
#            print syn_dists
            nsyn_dists = calc_distances(nsyn_pos)
            fourf_dists = calc_distances(fourf_pos)
#            print "nsyn dists"
#            print nsyn_pos
#            print nsyn_dists
            species_syn_dists[species] += syn_dists
            species_nsyn_dists[species] += nsyn_dists
            species_fourf_dists[species] += fourf_dists
#            print species_syn_dists
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for species in species_syn_dists.keys():
        outfile = open("%s/%s_nsyn_distances.txt" % (outdir, species), 'w')
        for distance in species_nsyn_dists[species]:
            outfile.write(str(distance) + "\n")
        outfile.close()
        outfile = open("%s/%s_syn_distances.txt" % (outdir, species), 'w')
        for distance in species_syn_dists[species]:
            outfile.write(str(distance) + "\n")
        outfile.close()
        outfile = open("%s/%s_fourf_distances.txt" % (outdir, species), 'w')
        for distance in species_fourf_dists[species]:
            outfile.write(str(distance) + "\n")
        outfile.close()


#    print species_nsyn_dists["SAL1"]
#    print species_syn_dists["SAL1"]
#    print len(species_syn_dists["SAL1"])
 #   print len(species_nsyn_dists["SAL1"])
    

def calc_distances(pos_list):
    dist_list = []
    x = 0
    used_pos = []
    for position in pos_list:
        
        for other_pos in pos_list:

#            print pos_list
            if position == other_pos or other_pos in used_pos:
                continue
            dist_list.append(abs(position - other_pos))
        used_pos.append(position)
    return dist_list

def sub_positions(inseq, ancestor):
    fourfold_list = changes.fourfold_codons()
    empty_chars = ["N", "-", "X"]
    insyns = []
    innsyns = []
    fourf = []
    x = 0
#    print inseq
#    print ancestor
    while x < len(inseq):
        incodon = inseq[x:x+3]
        anccodon = ancestor[x:x+3]
        if incodon == anccodon or "N" in incodon or "-" in incodon:
            x = x + 3
            continue
        inp = str(Seq.Seq(incodon.replace("-", "N")).translate())
        ancp = str(Seq.Seq(anccodon.replace("-", "N")).translate())
        if incodon[0] != anccodon[0]:
            change_index = x
        elif incodon[1] != anccodon[1]:
            change_index = x + 1
        elif incodon[2] != anccodon[2]:
            change_index = x + 2
        if inp == ancp:
            insyns.append(change_index)
        else:
            innsyns.append(change_index)
        if incodon in fourfold_list and anccodon in fourfold_list:
            fourf.append(change_index)
        x = x + 3
    return insyns, innsyns, fourf

def get_nearest_anc(index_dic, inseq_dic):
    nearest_anc = {}
    for inname in inseq_dic.keys():
        fewest_children = len(inseq_dic)
        smallest_node = -1
        for index, children in index_dic.items():
            if inname in children:
                if len(children) <= fewest_children:
                    fewest_children = len(children)
                    smallest_node = index
        nearest_anc[inname] = smallest_node
    return nearest_anc
