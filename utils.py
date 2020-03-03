import re
import json
from Bio.Phylo.PAML import baseml
import spotProblematicSeqsModules
import copy
import random
import numpy
import paml_tests
import os
import sys
import changes
from Bio.Align.Applications import MuscleCommandline
import multiprocessing
from multiprocessing import Pool
from potpour import Worker
import copy
import vcf
from Bio import SeqIO
from Bio import Seq
from Gene import Gene
from NCAR import NCAR
from gff import gffParser
import cPickle as pickle
import shutil
import subprocess
from glob import glob
from Bio.Phylo.PAML import codeml
import ete3
from ete3 import PhyloTree
from Bio.Phylo.PAML.chi2 import cdf_chi2
from scipy.stats import chisqprob
import statsmodels.stats.multitest as smm # this is required for a lot of stuff but is taking forever to load
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from StringIO import StringIO
import datetime
import scipy.stats

STOP_CODONS = ["TGA", "TAA", "TAG"]

def read_exec_paths(paths_file):
    paths_dic = {}
    reader = open(paths_file, 'rU')
    for line in reader:
        cur_line = line.split()
        paths_dic[cur_line[0]] = cur_line[1]
    return paths_dic

def read_params(params_file):
    reader = open(params_file, 'rU')
    cds_dic = {}
    for line in reader:
        cur_line = line.split()
        cds_dic[cur_line[0]] = cur_line[1]
    return cds_dic

def read_orthofile(orthoformat, ortho_file):
    if orthoformat == "orthofinder":
        ortho_dic = orthofinder_reader(ortho_file)
    elif orthoformat == "orthomcl":
        ortho_dic = orthomcl_reader(ortho_file)
    elif orthoformat == "proteinortho":
        ortho_dic = ortho_reader(ortho_file)
    else:
        print "Unknown orthogroup file format. Please specify orthofinder, orthomcl, or proteinortho format using the -w option."
        sys.exit()
    return ortho_dic

def ortho_reader(orthofile):
    #returns dictionary of dictionaries of orthologos genes. 
    #Lower level dictionaries have keys of species codes and 
    #values of lists of gene names. Upper level dictionaries have the
    #OG index (line number in orthofile) as keys.
    #species with paralogs are not included at all.
    reader = open(orthofile, 'rU')
    ortho_dic = {}
    counter = 10000
    for line in reader:
        if line.startswith("#"):
            continue
        gene_dic = {}
        cur_line = line.split()
        taxa_count = int(cur_line[0])
        seq_count = int(cur_line[1])
        for gene in cur_line[3:]:
            if gene == "*":
                continue
            cur_species = gene[0:4]
            if cur_species not in gene_dic.keys():
                gene_dic[cur_species] = []
                
            gene_dic[cur_species].append(gene)
        ortho_dic[counter] = gene_dic
        counter += 1
    return ortho_dic

def orthomcl_reader(orthofile):
    #returns dictionary of dictionaries of orthologous genes inferred by orthomcl. 
    #Lower level dictionaries have keys of species codes and 
    #values of lists of gene names. Upper level dictionaries have the
    #OG index (line number in orthofile) as keys.
    #species with paralogs are not included at all.
    reader = open(orthofile, 'rU')
    ortho_dic = {}
    for line in reader:
        gene_dic = {}
        cur_line = line.split()
        cur_og = int(cur_line[0].split("_")[1][:-1])
        for gene in cur_line[1:]:
            cur_species = gene.split("|")[1][0:4]
            if cur_species not in gene_dic.keys():
                gene_dic[cur_species] = []
                
            gene_dic[cur_species].append(gene.split("|")[1])
            gene_dic[cur_species] = [",".join(gene_dic[cur_species])]
        ortho_dic[cur_og] = gene_dic
    return ortho_dic

def orthofinder_reader(orthofile):
    #returns dictionary of dictionaries of orthologous genes inferred by OrthoFinder. 
    #Lower level dictionaries have keys of species codes and 
    #values of lists of gene names. Upper level dictionaries have the
    #OG index (line number in orthofile) as keys.
    #species with paralogs are not included at all.
    reader = open(orthofile, 'rU')
    ortho_dic = {}
    for line in reader:
        gene_dic = {}
        cur_line = line.split()
        cur_og = cur_line[0][4:-1]
        cur_og = int(cur_og) + 10000 #this is so that og numbers start at 10000. we don't want og numbers with different numbers of digits
        if len(cur_line) == 2:
            continue #skip orthogroups with a single sequence
        for gene in cur_line[1:]:
            cur_species = gene[0:4]
            if cur_species not in gene_dic.keys():
                gene_dic[cur_species] = []
                
            gene_dic[cur_species].append(gene)
            gene_dic[cur_species] = [",".join(gene_dic[cur_species])]
        ortho_dic[cur_og] = gene_dic
    return ortho_dic

def ncar_ortho_dic(ncar_file, min_taxa):
    reader = open(ncar_file, 'rU')
    ortho_dic = {}
    for line in reader:
        gene_dic = {}
        cur_line = line.split()
        cur_nmel = cur_line[0]
        cur_og = cur_line[1]
        if int(cur_line[2]) < min_taxa:
            continue
        seq_list = cur_line[3].split(",")
        for seq in seq_list:
            cur_species = seq[0:4]
            gene_dic[cur_species] = seq
        ortho_dic[cur_og] = gene_dic
    return ortho_dic

def write_ncars(ortho_dic, seq_dic, outdir, min_taxa, filtered_ortho_file):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    new_orthos = open(filtered_ortho_file, 'w')
    for ncar_num, ncar_seqs in ortho_dic.items():
        outfile = open("%s/ncar_%s.fasta" % (outdir, ncar_num), 'w')
        new_seq_list = []
        good_count = 0
        for species, seq in ncar_seqs.items():
            if len(seq_dic[species][seq]) < 1000:
                if len(seq_dic[species][seq].replace("N","").replace("n","")) >= 250:
                    good_count += 1
        if good_count < min_taxa:
            continue
        for species, seq in ncar_seqs.items():
            if len(seq_dic[species][seq]) < 1000:
                if len(seq_dic[species][seq].replace("N","").replace("n","")) >= 250:
                    outfile.write(">%s\n%s\n" % (seq, seq_dic[species][seq]))
                    new_seq_list.append(seq)
        new_orthos.write("NA\t%s\t%s\t%s\n" % (ncar_num, len(new_seq_list), ",".join(new_seq_list)))
        outfile.close()
    new_orthos.close()

def write_ncar_cnees(ortho_dic, seq_dic, outdir, min_taxa, filtered_ortho_file):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    new_orthos = open(filtered_ortho_file, 'w')
    for ncar_num, ncar_seqs in ortho_dic.items():
        outfile = open("%s/ncar_%s.fasta" % (outdir, ncar_num), 'w')
        new_seq_list = []
        good_count = 0
        for species, seq in ncar_seqs.items():
            if len(seq_dic[species][seq].replace("N","").replace("n","")) >= 100:
                    good_count += 1
        if good_count < min_taxa:
            continue
        for species, seq in ncar_seqs.items():
            if len(seq_dic[species][seq].replace("N","").replace("n","")) >= 100:
                    outfile.write(">%s\n%s\n" % (seq, seq_dic[species][seq]))
                    new_seq_list.append(seq)
        new_orthos.write("NA\t%s\t%s\t%s\n" % (ncar_num, len(new_seq_list), ",".join(new_seq_list)))
        outfile.close()
    new_orthos.close()


def get_og_num(query_gene, ortho_dic):
    species = query_gene[0:4]
#    query_gene = query_gene + "-RA"
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

def read_species_pickle(target_species, pickle_dir, gene_or_ncar):
#    pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle_test" % (target_species, target_species), 'rb')
#    pickle_file = open("/scratch/tmp/berubin/resequencing/hic/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species), 'rb')
#    pickle_file = open("%s/testset_%s_%s_vcf_dic.pickle" % (pickle_dir, target_species, gene_or_ncar), 'rb')
    pickle_file = open("%s/%s_%s_vcf_dic.pickle" % (pickle_dir, target_species, gene_or_ncar), 'rb')
    gene_objects = pickle.load(pickle_file)
    pickle_file.close()
    return gene_objects

def harvest_ncar_worker(attribute_list):
    cur_ncar = attribute_list[0]
#    if gene_name not in ["LZEP_01797"]:
#        return
    scaf_seq = attribute_list[1]
    target_species = attribute_list[2]
    vcf_reader = vcf.Reader(filename = "/Genomics/kocherlab/berubin/Dropbox/data/freebayes_v2.1/%s_filtered_miss.vcf.gz" % (target_species))
#    vcf_reader = vcf.Reader(filename = "/Genomics/kocherlab/berubin/Dropbox/data/freebayes_v2.1/%s_chr_7_11143675_11696910.vcf" % target_species)
    cur_ncar.set_sequence(scaf_seq, 0)
    try:
        reader = vcf_reader.fetch(cur_ncar.scaf, cur_ncar.start, cur_ncar.end)
        cur_ncar.get_genotypes(vcf_reader, 0)
    except ValueError:
        print "can't get that bit of vcf: %s" % cur_ncar.name
        cur_ncar.alts = {}
        cur_ncar.refs = {}
        cur_ncar.called_counts = {}
    return cur_ncar


def get_species_ncar_data(target_species, num_threads, basedir, ncar_gff, genome_fasta, ortho_dic, gene_dic, pickle_dir):
#    if os.path.exists("/scratch/tmp/berubin/resequencing/hic/%s/genotyping/%s_ncar_vcf_dic.pickle" % (target_species, target_species)):
#    if os.path.exists("%s/%s_cnee_vcf_dic.pickle" % (pickle_dir, target_species)):
#    if os.path.exists("%s/%s_ncar_vcf_dic.pickle" % (pickle_dir, target_species)):
    if os.path.exists("%s/%s_cactusncar_vcf_dic.pickle" % (pickle_dir, target_species)):
        print "%s NCAR pickle exists" % target_species
        return
    print "Building %s NCAR pickle" % target_species
    seq_dic = {}
    reader = SeqIO.parse(genome_fasta, format = 'fasta')
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    print "genome read"
    all_ncars_list = []
    for og_num, ncar_dic in ortho_dic.items():
#        all_ncars_list = all_ncars_list + ncar_dic[target_species]
        if not ncar_dic.get(target_species, None) == None:
            all_ncars_list.append(ncar_dic[target_species])
    gff_file = open(ncar_gff, 'rU')
    print "read gff"
    ncar_list = []
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []

    for line in gff_file:
        cur_line = line.split()
        cur_scaf = cur_line[0]
#        if cur_scaf != "NMEL_chr_1":# and cur_scaf != "LVIE_chr_4":
#            continue
        cur_start = int(cur_line[3])
        cur_end = int(cur_line[4])
#        if cur_start < 11143675 or cur_end > 11696910:
#            continue
        cur_strand = cur_line[6]
        cur_name = cur_line[8].split(";")[0].split("ID=")[1]
        if cur_name not in all_ncars_list:
            if cur_name == "LZEP:LZEP_chr_4:13156535:13156885:+":
                print "lost to all_ncars_list"
            print "not in orthodic: %s" % cur_name
            continue
#####TURN THIS ON IF DOING NCARS
#        if cur_end - cur_start > 1000:
#            if cur_name == "LZEP:LZEP_chr_4:13156535:13156885:+":
#                print "lost to length"
#            print "too long: %s" % cur_name
#            continue
        cur_ncar = NCAR(cur_name, cur_scaf, cur_start, cur_end, cur_strand)
#        if cur_name != "LZEP:LZEP_chr_4:13156535:13156885:+":
#            continue
        work_list.append([cur_ncar, seq_dic[cur_scaf], target_species])
#        ncar= harvest_ncar_worker([cur_ncar, seq_dic[cur_scaf], target_species])
    ncar_list = pool.map(harvest_ncar_worker, work_list)
    pool.close()
    pool.join()
    neighbor_work = []
    for ncar in ncar_list:
        ncar.set_neighbors(gene_dic.values())
    print "ncar data gotten"
    print str(datetime.datetime.now())
    print len(ncar_list)
    all_samples = []
    ncar_objects = {}
    for ncar in ncar_list:
        ncar_objects[ncar.name] = ncar
        for sample_dic in ncar.sample_gts.values():
            for sample in sample_dic.keys():
                all_samples.append(sample)
            all_samples = list(set(all_samples))
    print "Dumping %s pickle" % target_species
#    pickle_file = open("%s/%s_cnee_vcf_dic.pickle" % (pickle_dir, target_species), 'wb')
    pickle_file = open("%s/%s_cactusncar_vcf_dic.pickle" % (pickle_dir, target_species), 'wb')
    pickle.dump(ncar_objects, pickle_file)
    pickle_file.close()
    sys.exit()

def get_species_data(target_species, num_threads, basedir, ingenome, coding_gff, pickle_dir):
    #This is a big method. Harvests gene coordinates from GFF3 files
    #and creates Gene objects with all of their characteristics.
#    if os.path.exists("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle_test" % (target_species, target_species)):
#    if os.path.exists("/scratch/tmp/berubin/resequencing/hic/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species)):
    if os.path.exists("%s/%s_gene_vcf_dic.pickle" % (pickle_dir, target_species)):
        print "%s pickle exists" % target_species
        return
#        pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species), 'rb')

        
#        gene_objects = pickle.load(pickle_file)
#        pickle_file.close()
#        return gene_objects
    print "Building %s pickle" % target_species
#    official_dir = "/Genomics/kocherlab/berubin/official_release_v2.1"
    seq_dic = {}
#    if target_species == "LALB":
#        reader = SeqIO.parse("%s/%s/%s_v3/%s/%s_genome_v3.0.fasta" % (official_dir, target_species, target_species, target_species, target_species), format = 'fasta')
#    if target_species in ["BIMP", "PMIB"]:
#        reader = SeqIO.parse("%s/%s/%s_genome_v2.0.fasta" % (official_dir, target_species, target_species), format = 'fasta')
#    else:
#        reader = SeqIO.parse("%s/%s/%s_genome_v2.1.fasta" % (official_dir, target_species, target_species), format = 'fasta')
    reader = SeqIO.parse(ingenome, format = 'fasta')
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    cds_dic = {}
#    if target_species == "LALB":
#        reader = SeqIO.parse("%s/%s/%s_v3/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (official_dir, target_species, target_species, target_species, target_species), format = 'fasta')
#    else:
#    reader = SeqIO.parse("%s/%s/%s_OGS_v2.1_longest_isoform.cds.fasta" % (official_dir, target_species, target_species), format = 'fasta')

#    for rec in reader:
#        cds_dic[rec.id[:-3]] = str(rec.seq)

#    gff_file = gffParser(open("%s/%s/%s_testset.gff3" % (official_dir, target_species, target_species), 'rU'))
#    if target_species == "LALB":
#        gff_file = gffParser(open("%s/%s/%s_v3/%s/%s_OGS_v1.0_longest_isoform.gff3" % (official_dir, target_species, target_species, target_species, target_species), 'rU'))
#    else:
#    gff_file = gffParser(open("%s/%s/%s_OGS_v2.1_longest_isoform.gff3" % (official_dir, target_species, target_species), 'rU'))
    gff_file = gffParser(open(coding_gff, 'rU'))
    print "read gff"
#    gff_file = gffParser(open("%s/%s/%s_09600.gff3" % (official_dir, target_species, target_species), 'rU'))
    print "get gene data"
    print str(datetime.datetime.now())
    gene_dic = gff_file.geneDict()
    gene_objects = {}
#    work_queue = multiprocessing.Queue()
#    result_queue = multiprocessing.Queue()
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    name_list = []
    dic_list = []
    gff_list = []
    seq_list = []
    species_list = []
    for gene_name in gene_dic.keys():
#        if gene_name == "APUR_11398": #"LZEP_05457":""LMAL_00737" or gene_name == "LLEU_03402":
#        if gene_dic[gene_name][0] == "LZEP_chr_7" or gene_dic[gene_name][0] == "LVIE_chr_4": 
        if True:
            work_list.append([gene_name, gene_dic, gff_file, seq_dic, target_species])
#            harvest_worker([gene_name, gene_dic, gff_file, seq_dic, target_species])
    gene_list = pool.map(harvest_worker, work_list)
    pool.close()
    pool.join()
    neighbor_work = []
    for gene in gene_list:
        gene.set_neighbors(gene_list)
    print "gene data gotten"
    print str(datetime.datetime.now())
    print len(gene_list)
    all_samples = []
    for gene in gene_list:
        gene_objects[gene.name] = gene
        for sample_dic in gene.sample_gts.values():
            for sample in sample_dic.keys():
                all_samples.append(sample)
            all_samples = list(set(all_samples))
    print "Dumping %s pickle" % target_species
#    pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle_test" % (target_species, target_species), 'wb')
#    pickle_file = open("/scratch/tmp/berubin/resequencing/hic/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species), 'wb')
    pickle_file = open("%s/%s_gene_vcf_dic.pickle" % (pickle_dir, target_species), 'wb')
    pickle.dump(gene_objects, pickle_file)
    pickle_file.close()
    if not os.path.exists("%s/pnps_files/" % (basedir)):
        os.mkdir("%s/pnps_files/" % (basedir))
    pnps_file = open("%s/pnps_files/%s_pnps.txt" % (basedir, target_species), 'w')
    outline = "Gene"
    print all_samples
    for sample in all_samples:
        outline = outline + "\t" + sample + "_pn" + "\t" + sample + "_ps" + "\t" + sample + "_ngenod" + "\t" + sample + "_sgenod"
    outline = outline + "\t" + "#nsynsites\t#synsites\tpotennsyn\tpotensyn"
    pnps_file.write(outline + "\n")
    for gene in gene_list:
        outline = gene.name
        for sample in all_samples:
            if sample in gene.sample_nsyn_counts.keys():
                outline = outline + "\t" + str(gene.sample_nsyn_counts[sample])
            else:
                outline = outline + "\t" + "0"
            if sample in gene.sample_syn_counts.keys():
                outline = outline + "\t" + str(gene.sample_syn_counts[sample])
            else:
                outline = outline + "\t" + "0"
            if sample in gene.sample_nsyn_genod.keys():
                outline = outline + "\t" + str(gene.sample_nsyn_genod[sample])
            else:
                outline = outline + "\t" + "0"
            if sample in gene.sample_syn_genod.keys():
                outline = outline + "\t" + str(gene.sample_syn_genod[sample])
            else:
                outline = outline + "\t" + "0"
        outline = outline + "\t" + str(len(gene.nsyn_coords)) + "\t" + str(len(gene.syn_coords)) + "\t" + str(round(gene.potent_nsyn, 3)) + "\t" + str(round(gene.potent_syn, 3))
        pnps_file.write(outline + "\n")
    pnps_file.close()

    pnps_file = open("%s/pnps_files/%s_pnps_high.txt" % (basedir, target_species), 'w')
    outline = "Gene"
    print all_samples
    for sample in all_samples:
        outline = outline + "\t" + sample + "_pn" + "\t" + sample + "_ps" + "\t" + sample + "_ngenod" + "\t" + sample + "_sgenod"
    outline = outline + "\t" + "#nsynsites\t#synsites\tpotennsyn\tpotensyn"
    pnps_file.write(outline + "\n")
    for gene in gene_list:
        outline = gene.name
        for sample in all_samples:
            if sample in gene.sample_hi_counts.keys():
                outline = outline + "\t" + str(gene.sample_hi_counts[sample])
            else:
                outline = outline + "\t" + "0"
            if sample in gene.sample_syn_counts.keys():
                outline = outline + "\t" + str(gene.sample_syn_counts[sample])
            else:
                outline = outline + "\t" + "0"
            if sample in gene.sample_hi_genod.keys():
                outline = outline + "\t" + str(gene.sample_hi_genod[sample])
            else:
                outline = outline + "\t" + "0"
            if sample in gene.sample_syn_genod.keys():
                outline = outline + "\t" + str(gene.sample_syn_genod[sample])
            else:
                outline = outline + "\t" + "0"
        outline = outline + "\t" + str(len(gene.hi_impact_coords)) + "\t" + str(len(gene.syn_coords)) + "\t" + str(round(gene.potent_nsyn, 3)) + "\t" + str(round(gene.potent_syn, 3))
        pnps_file.write(outline + "\n")
    pnps_file.close()

    pnps_file = open("%s/pnps_files/%s_pnps_low.txt" % (basedir, target_species), 'w')
    outline = "Gene"
    print all_samples
    for sample in all_samples:
        outline = outline + "\t" + sample + "_pn" + "\t" + sample + "_ps" + "\t" + sample + "_ngenod" + "\t" + sample + "_sgenod"
    outline = outline + "\t" + "#nsynsites\t#synsites\tpotennsyn\tpotensyn"
    pnps_file.write(outline + "\n")
    for gene in gene_list:
        outline = gene.name
        for sample in all_samples:
            if sample in gene.sample_low_counts.keys():
                outline = outline + "\t" + str(gene.sample_low_counts[sample])
            else:
                outline = outline + "\t" + "0"
            if sample in gene.sample_syn_counts.keys():
                outline = outline + "\t" + str(gene.sample_syn_counts[sample])
            else:
                outline = outline + "\t" + "0"
            if sample in gene.sample_low_genod.keys():
                outline = outline + "\t" + str(gene.sample_low_genod[sample])
            else:
                outline = outline + "\t" + "0"
            if sample in gene.sample_syn_genod.keys():
                outline = outline + "\t" + str(gene.sample_syn_genod[sample])
            else:
                outline = outline + "\t" + "0"
        outline = outline + "\t" + str(len(gene.low_impact_coords)) + "\t" + str(len(gene.syn_coords)) + "\t" + str(round(gene.potent_nsyn, 3)) + "\t" + str(round(gene.potent_syn, 3))
        pnps_file.write(outline + "\n")
    pnps_file.close()


    gene_objects = {}
    sys.exit()

def make_neighborhoods(attribute_list):
    gene = attribute_list[0]
    print "neighbors" + "\t" + gene.name
    gene_list = attribute_list[1]
    gene.set_neighbors(gene_list)

def harvest_worker(attribute_list):#gene_name, gene_dic, gff_file, seq_dic, target_species):
    gene_name = attribute_list[0]
    print gene_name
#    if gene_name not in ["LZEP_01797"]:
#        return
    gene_dic = attribute_list[1]
    gff_file = attribute_list[2]
    seq_dic = attribute_list[3]
    target_species = attribute_list[4]
    if target_species == "LALB":
#        vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_social_renamed.vcf.gz" % (target_species, target_species))
        vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_solitary_renamed.vcf.gz" % (target_species, target_species))
    elif target_species in ["BIMP", "PMIB"]:
        vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/bombus.bmel.vcf.gz" % (target_species))
    else:
#        vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_filtered_miss.vcf.gz" % (target_species, target_species))
        vcf_reader = vcf.Reader(filename = "/Genomics/kocherlab/berubin/Dropbox/data/freebayes_v2.1/%s_filtered_miss.vcf.gz" % (target_species))
#        vcf_reader = vcf.Reader(filename = "/Genomics/kocherlab/berubin/Dropbox/data/freebayes_v2.1/%s_nab2.vcf.gz" % (target_species))
#        vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_testset.vcf.gz" % (target_species, target_species))
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
    cur_gene.get_first_intron_sequence()
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
#    print inspecies
#    print og_num
#    print ortho_dic[og_num]
    if inspecies not in ortho_dic[og_num].keys() or outspecies not in ortho_dic[og_num].keys():
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
        if len(ortho_dic[og_num][inspecies]) == 1 and len(ortho_dic[og_num][outspecies]) == 1:
            og_list.append(og_num)
#    og_list = [2110]
#    og_list = [870]
#    og_list = og_list[0:50]
    print "starting ogs: %s" % len(og_list)
    print "orthodic ogs: %s" % len(ortho_dic.keys())
    work_list = []
    in_gene_dic = read_species_pickle(inspecies)
    print "%s pickle read" % inspecies
    out_gene_dic = read_species_pickle(outspecies)
    print "%s pickle read" % outspecies
    og_map_dic = make_og_gene_map(ortho_dic)
    not_enough_neighbors = 0
#    og_list = [9860]
    for og in og_list:
        inspecies_gene = ortho_dic[og][inspecies][0][:-3]
        outspecies_gene = ortho_dic[og][outspecies][0][:-3]
        ##THIS IS JUST FOR TESTSET
#        if inspecies_gene not in in_gene_dic.keys() or outspecies_gene not in out_gene_dic.keys():
#            continue
        print og
        in_gene = in_gene_dic[inspecies_gene]
        out_gene = out_gene_dic[outspecies_gene]
        neighbor_dic = {}
        in_gene.neighbors.append(in_gene.name)
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
                ###THIS IS ALSO JUST FOR TESTSET
#                if neighbor_name not in in_gene_dic.keys() or outneighbor_name not in out_gene_dic.keys():
#                    continue
                neighbor_dic[neighbor_og] = (in_gene_dic[neighbor_name], out_gene_dic[outneighbor_name])
        if len(neighbor_dic) < 4:
            not_enough_neighbors += 1
            continue
        og_map.write("OG_%s\t%s\t%s\n" % (og, inspecies_gene, outspecies_gene))
        og_map.flush()
#        hka_line_list.append(hka_worker([inspecies, outspecies, seq_type, og, in_gene, out_gene, out_path, align_dir, neighbor_dic]))

        work_list.append([inspecies, outspecies, seq_type, og, in_gene, out_gene, out_path, align_dir, neighbor_dic])
    print "not enough neighbors: %s" % not_enough_neighbors
    print len(work_list)
    og_map.close()
    pool = multiprocessing.Pool(processes = num_threads)
    hka_lines = pool.map(hka_worker, work_list)
    pool.close()
    pool.join()
    full_length_lines = []
    jody_lines = []
    hkadirect_lines = []
    direct_pairwise_lines = []
    too_short = 0
    too_different_1 = 0
    too_different_2 = 0
    too_different_3 = 0
    insuf_sam = 0
    no_align_1 = 0
    no_align_2 = 0
    no_align_3 = 0
    no_align_4 = 0
    print "starting lines: %s" % len(hka_lines)
    for line in hka_lines:
        if "too_short" in line:
            too_short += 1
        if "too_different_1" in line:
            too_different_1 += 1

        if "too_different_2" in line:
            too_different_2 += 1
        if "too_different_3" in line:
            too_different_3 += 1

        if "insufficient_sampling" in line:
            insuf_sam +=1
        if "no_alignment_1" in line:
            no_align_1 += 1
        if "no_alignment_2" in line:
            no_align_2 += 1
        if "no_alignment_3" in line:
            no_align_3 += 1
        if "no_alignment_4" in line:
            no_align_4 += 1

        if "too_short" not in line and "too_different" not in line and "HKA_failed" not in line[0] and "insufficient_sampling" not in line and "no_alignment" not in line:
            jody_lines.append(line[0])
    print "too short: %s" % too_short
    print "too diff 1: %s" % too_different_1
    print "too diff 2: %s" % too_different_2
    print "too diff 3: %s" % too_different_3
    print "insuff: %s" % insuf_sam
    print "no align 1: %s" % no_align_1
    print "no align 2: %s" % no_align_2
    print "no align 3: %s" % no_align_3
    print "no align 4: %s" % no_align_4

#    jody_correction(jody_lines, "%s/hka_tests/%s_%s_%s_hka_results_correct.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_results_correct_faster.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_results_correct_slower.txt" % (out_path, inspecies, outspecies, seq_type))
    for line in hka_lines:
        if "too_short" not in line and "too_different" not in line and "insufficient_sampling" not in line and "no_alignment" not in line:
            if len(line[1]) > 0:
                if "FAILURE" not in line[2]:
                    direct_pairwise_lines.append(line[2])
    pairwise_direct_correction(direct_pairwise_lines, "%s/hka_tests/%s_%s_%s_hka_direct_pairwise_results_correct.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_direct_pairwise_results_correct_faster.txt" % (out_path, inspecies, outspecies, seq_type), "%s/hka_tests/%s_%s_%s_hka_direct_pairwise_results_correct_slower.txt" % (out_path, inspecies, outspecies, seq_type))

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
    pval_corr = smm.multipletests(p_list, alpha = 0.1, method = 'fdr_bh')[1]
    for x in range(len(pval_corr)):
        chi_list[x].append(pval_corr[x])
    chi_list.sort(key = lambda x: x[-2])
    outfile = open(outfile, 'w')
    outfile.write("IDloci\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texpHKA_S\texpVar_S\texpHKA_div\texpVar_D\texpHKA_theta\tpartialHKA\tp_val\tbh5_p\n")
    for locus in chi_list:
        outfile.write("\t".join([str(i) for i in locus]) + "\n")
    outfile.close()
    faster = open(fasterfile, 'w')
    faster.write("IDloci\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texpHKA_S\texpVar_S\texpHKA_div\texpVar_D\texpHKA_theta\tpartialHKA\tp_val\tbh5_p\tobs_exp\n")
    slower = open(slowerfile, 'w')
    slower.write("IDloci\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texpHKA_S\texpVar_S\texpHKA_div\texpVar_D\texpHKA_theta\tpartialHKA\tp_val\tbh5_p\tobs_exp\n")
    for locus in chi_list:
#        if float(locus[-1]) < 0.01:
        obs_s = float(locus[2])
        exp_s = float(locus[7])
        if obs_s > exp_s:
            faster.write("\t".join([str(i) for i in locus]) + "\t" + str(obs_s / exp_s) + "\n")
        else:
            slower.write("\t".join([str(i) for i in locus]) + "\t" + str(obs_s / exp_s) + "\n")
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
    pval_corr = smm.multipletests(p_list, alpha = 0.1, method = 'fdr_bh')[1]
    for x in range(len(pval_corr)):
        chi_list[x].append(pval_corr[x])
    chi_list.sort(key = lambda x: x[-2])
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
#        if float(locus[-1]) < 0.01:
        obs_s = float(locus[2])
        exp_s = float(locus[7])
        if obs_s > exp_s:
            faster.write("\t".join([str(i) for i in locus]) + "\t" + str(obs_s / exp_s) + "\n")
        else:
            slower.write("\t".join([str(i) for i in locus]) + "\t" + str(obs_s / exp_s) + "\n")
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
    pval_corr = smm.multipletests(p_list, alpha = 0.05, method = 'fdr_bh')[1]
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
    if not os.path.exists("%s/og_cds_%s.afa" % (align_dir, og_num)):
        return "OG_%s\tinsufficient_sampling\t" % og_num
    inpoly, outpoly, inseq, outseq, insample, outsample = gather_hka_data(in_gene, out_gene, seq_type)
    if len(inseq) < 50 or len(outseq) < 50:
        return "OG_%s\ttoo_short\n" % (og_num)
    if in_gene.strand == -1:
        inseq = str(Seq.Seq(inseq).reverse_complement())
    if out_gene.strand == -1:
        outseq = str(Seq.Seq(outseq).reverse_complement())
#    if len(outseq) < len(inseq):
#        inseq = inseq[0:len(outseq)]
#    if len(inseq) < len(outseq):
#        outseq = outseq[0:len(inseq)]
    pairwise_data = fsa_pairwise_diff_count(inseq, outseq, inspecies, outspecies, og_num, in_gene, out_gene, out_path, seq_type)
    if pairwise_data == "no_alignment_1":
        return "OG_%s\tno_alignment_1\n" % (og_num)
    if pairwise_data == "no_alignment_2":
        return "OG_%s\tno_alignment_2\n" % (og_num)
    if pairwise_data == "no_alignment_3":
        return "OG_%s\tno_alignment_3\n" % (og_num)
    if pairwise_data == "no_alignment_4":
        return "OG_%s\tno_alignment_4\n" % (og_num)

    average_diff = pairwise_data[0]
    align_len = pairwise_data[1]
    inpoly = pairwise_data[2]
    outpoly = pairwise_data[3]
    inlen = pairwise_data[4]
    outlen = pairwise_data[5]
#    if align_len < 500 or inlen < 500:
    if align_len < 50 or inlen < 50:
        return "OG_%s\ttoo_different_1\n" % (og_num)
    if 1.0 * average_diff / align_len > 0.10:
        return "OG_%s\ttoo_different_2\t%s\n" % (og_num, 1.0 * average_diff / align_len)
    if 1.0 * align_len / inlen < 0.40:
        return "OG_%s\ttoo_different_3\t%s\n" % (og_num, 1.0 * align_len / inlen)
    if insample < 5:# or outsample < 5:
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
#    FNULL = open(os.devnull, 'w')
#    subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
    os.chdir(switch_dir)
#    p_val, obs_poly, exp_poly = parse_hka_result(inspecies, outspecies, seq_type, "%s/hka_tests/hey_%s_%s/OG_%s_%s_%s_%s_hka_results.hka" % (out_path, inspecies, outspecies, og_num, inspecies, outspecies, seq_type))
    p_val = obs_poly = exp_poly = -1
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
        if "ERROR: Results are not available. Time obtained is negative." in line:
            return "FAILURE"
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

def ingene_fourfold_degenerate_string(attribute_list):
    #returns two tuples of data about the fourfold degenerate
    #sites in a particular gene alignment. the first two is a set of
    #data that will be used for HKAdirect. The second tuple is just
    #a single string formatted for use by Hey's HKA
    ingene = attribute_list[0]
    inspecies = attribute_list[1]
    inseq = attribute_list[2]
    outseq = attribute_list[3]
    fix_fourfold = count_ingene_fourfold(inseq, outseq, ingene)
    in_poly_fourf = ingene.fourfold
    in_potent_fourfold = ingene.potent_fourfold
    align_potent_fourf = potential_aligned_fourfold_sites(inseq, outseq)
    
    return fix_fourfold, in_poly_fourf, in_potent_fourfold, align_potent_fourf

def sum_ingene_fourfs(fourf_list):
    #sums all of the variables for fourfold degenerate sites
    #across a bunch of genes and makes them into a string
    #for HKAdirect
    fix_fourfold = in_poly_fourf = in_potent_fourfold = align_potent_fourf = 0
    x = 0
    num_loci = 0.0
    for fourf_stats in fourf_list:
        fix_fourfold += fourf_stats[0]
        in_poly_fourf += fourf_stats[1]
        in_potent_fourfold += fourf_stats[2]
        align_potent_fourf += fourf_stats[3]
        num_loci += 1
    return float(in_poly_fourf), float(fix_fourfold), float(in_potent_fourfold), float(align_potent_fourf)

    
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

def get_fsa_aligned_seqs(align_dir, og_num, inspecies, outspecies, gene_or_ncar):
    #get the aligned sequences for two species in a particular OG
    if gene_or_ncar == "gene":
        reader = SeqIO.parse("%s/og_cds_%s.afa" % (align_dir, og_num), format = 'fasta')
    elif gene_or_ncar == "ncar":
        reader = SeqIO.parse("%s/ncar_%s.afa" % (align_dir, og_num), format = 'fasta')
    for rec in reader:
        if rec.id[0:4] == inspecies:
            inseq = str(rec.seq).upper()
            inseq_name = rec.id[:-3]
        elif rec.id[0:4] == outspecies:
            outseq = str(rec.seq).upper()
            outseq_name = rec.id[:-3]
    return inseq, outseq

def get_fsa_aligned_dic(align_dir, og_num, gene_or_ncar):
    #get the aligned sequences for two species in a particular OG
    if gene_or_ncar == "ncar":
        reader = SeqIO.parse("%s/ncar_%s.afa" % (align_dir, og_num), format = 'fasta')
    elif gene_or_ncar == "gene":
        reader = SeqIO.parse("%s/og_cds_%s.afa" % (align_dir, og_num), format = 'fasta')
    seq_dic = {}
    for rec in reader:
        seq_dic[rec.id[0:4]] = str(rec.seq).upper()
    return seq_dic


def get_prank_aligned_dic(align_dir, og_num):
    #get the aligned sequences for two species in a particular OG
    reader = SeqIO.parse("%s/og_cds_%s.1.fas" % (align_dir, og_num), format = 'fasta')
    seq_dic = {}
    for rec in reader:
        seq_dic[rec.id[0:4]] = str(rec.seq)
    return seq_dic

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
    elif seq_type == "first_intron":
        inpoly = inspecies_gene.first_intron_subs
        outpoly = outspecies_gene.first_intron_subs
        inseq = inspecies_gene.first_intron_seq
        outseq = outspecies_gene.first_intron_seq
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
    if intuple[0] == 0 and intuple[1] == 0 and outtuple[1] ==0 and outtuple[0] == 0:
        return "no_alignment_3"

    if intuple[0] == -1 and intuple[1] == -1 and outtuple[1] == -1 and outtuple[0] == -1:
        return "no_alignment_4"
    if (intuple[1] - intuple[0]) < 2000 or (outtuple[1] - outtuple[0]) < 2000:
#        print "too different alignment: %s" % og_num
        return "no_alignment_1"
    if abs((intuple[1] - intuple[0]) - (outtuple[1] - outtuple[0])) > 500:
        return "no_alignment_2"
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

def fsa_pairwise_diff_count(seq1, seq2, inspecies, outspecies, og_num, in_gene, out_gene, out_path, seq_type):
    #Align two sequences using mafft and return the number of 
    #differences between them. Since the number of polymorphisms
    #counted in these regions is dependent on where the alignable
    #regions are, this also provides the number of polymorphisms
    if not os.path.exists("%s/conservation_files_%s" % (out_path, seq_type)):
        os.mkdir("%s/conservation_files_%s" % (out_path, seq_type))

    if not os.path.exists("%s/conservation_files_%s/%s_%s/" % (out_path, seq_type, inspecies, outspecies)):
        os.mkdir("%s/conservation_files_%s/%s_%s" % (out_path, seq_type, inspecies, outspecies))
    intuple, outtuple = conserved_noncoding(seq1, seq2, og_num, inspecies, outspecies, "%s/conservation_files_%s/%s_%s" % (out_path, seq_type, inspecies, outspecies))
    if intuple[0] == 0 and intuple[1] == 0 and outtuple[1] ==0 and outtuple[0] == 0:
        return "no_alignment_3"

    if intuple[0] == -1 and intuple[1] == -1 and outtuple[1] == -1 and outtuple[0] == -1:
        return "no_alignment_4"
    if (intuple[1] - intuple[0]) < 50 or (outtuple[1] - outtuple[0]) < 50:
#        print "too different alignment: %s" % og_num
        return "no_alignment_1"
#    if abs((intuple[1] - intuple[0]) - (outtuple[1] - outtuple[0])) > 500:
#        return "no_alignment_2"
    if not os.path.exists("%s/fsa_files_%s" % (out_path, seq_type)):
        os.mkdir("%s/fsa_files_%s" % (out_path, seq_type))

    if not os.path.exists("%s/fsa_files_%s/%s_%s" % (out_path, seq_type, inspecies, outspecies)):
        os.mkdir("%s/fsa_files_%s/%s_%s" % (out_path, seq_type, inspecies, outspecies))
    seqfile = open("%s/fsa_files_%s/%s_%s/OG_%s_%s_%s.fa" % (out_path, seq_type, inspecies, outspecies, og_num, inspecies, outspecies), 'w')
#    seqfile.write(">inseq\n%s\n" % (seq1[intuple[0]:intuple[1]]))
#    seqfile.write(">outseq\n%s\n" % (seq2[outtuple[0]:outtuple[1]]))
    seqfile.write(">inseq\n%s\n" % (seq1))
    seqfile.write(">outseq\n%s\n" % (seq2))
    intuple = (1, len(seq1))
    outtuple = (1, len(seq2))

    seqfile.close()
    cmd = ["/Genomics/kocherlab/berubin/local/src/fsa-1.15.9/bin/fsa", "%s/fsa_files_%s/%s_%s/OG_%s_%s_%s.fa" % (out_path, seq_type, inspecies, outspecies, og_num, inspecies, outspecies)]
    FNULL = open(os.devnull, 'w')
    with open("%s/fsa_files_%s/%s_%s/OG_%s_%s_%s.afa" % (out_path, seq_type, inspecies, outspecies, og_num, inspecies, outspecies), 'w') as outfile:
        subprocess.call(cmd, stdout = outfile, stderr = FNULL)
    outfile.close()
    align_dic = {}
    reader = SeqIO.parse("%s/fsa_files_%s/%s_%s/OG_%s_%s_%s.afa" % (out_path, seq_type, inspecies, outspecies, og_num, inspecies, outspecies), format = 'fasta')
    for rec in reader:
        align_dic[rec.id] = str(rec.seq)
    average_diff = 0
    missing_data = ["N", "-", "n"]
    indel_count = 0

    if seq_type == "flank":
        if in_gene.strand == -1:
############GGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAHHHHHHHHHH
#            in_gene_alt_dic = in_gene.check_polymorphisms_fixed(align_dic, in_gene.end, "inseq")
            in_gene_alt_dic = in_gene.check_polymorphisms_fixed(align_dic, in_gene.flank_end - intuple[1], "inseq")
        else:
#            in_gene_alt_dic = in_gene.check_polymorphisms_fixed(align_dic, in_gene.flank_end - intuple[1], "inseq")
            in_gene_alt_dic = in_gene.check_polymorphisms_fixed(align_dic, in_gene.flank_start + intuple[0], "inseq")
        if out_gene.strand == -1:
            out_gene_alt_dic = out_gene.check_polymorphisms_fixed(align_dic, out_gene.flank_end - outtuple[1], "outseq")
        else:
            out_gene_alt_dic = out_gene.check_polymorphisms_fixed(align_dic, out_gene.flank_start + outtuple[0], "outseq")

    elif seq_type == "intron":
        in_gene_alt_dic = in_gene.coding_fixed_align(align_dic["inseq"])
        out_gene_alt_dic = out_gene.coding_fixed_align(align_dic["outseq"])
    elif seq_type == "first_intron":
        if in_gene.strand == -1:
            in_gene_alt_dic = in_gene.check_polymorphisms_fixed(align_dic, in_gene.cds_coords[-1][0] - intuple[1], "inseq")
        else:
            in_gene_alt_dic = in_gene.check_polymorphisms_fixed(align_dic, in_gene.cds_coords[0][1] + intuple[0], "inseq")
        if out_gene.strand == -1:
            out_gene_alt_dic = out_gene.check_polymorphisms_fixed(align_dic, out_gene.cds_coords[-1][0] - outtuple[1], "outseq")
        else:
            out_gene_alt_dic = out_gene.check_polymorphisms_fixed(align_dic, out_gene.cds_coords[0][1] + outtuple[0], "outseq")
        print in_gene.name
        for k in sorted(in_gene_alt_dic.keys()):
            print "%s:%s" % (k, in_gene_alt_dic[k])
        print in_gene.name
        for k in sorted(out_gene_alt_dic.keys()):
            print "%s:%s" % (k, out_gene_alt_dic[k])


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


def includes_paralogs(og_dic):
    for k, v in og_dic.items():
        if v[0].count(",") > 0:
            return True

def prop_shared_sequence(inseq, outseq):
    #counts the number of sites in two aligned sequences where
    #both have bases rather than gaps. returns the proportion
    #of the total bases (not including gaps) in the sequence
    #with a lower proportion shared
    bases = ["A", "C", "T", "G"]
    shared_count = 0.0
    for x in range(len(inseq)):
        if inseq[x].upper() in bases and outseq[x].upper() in bases:
            shared_count += 1
    inlen = len(inseq) - inseq.count("-")
#    outlen = len(outseq) - outseq.count("-")
    inshared = shared_count / inlen
#    outshared = shared_count / outlen
#    if inshared > outshared:
#        return outshared
#    else:
    return inshared

def ncar_pairwise_diff_count(inspecies, seq1, seq2, in_ncar):
    in_ncar_alt_dic = in_ncar.noncoding_fixed_align(seq1)    
    x = 0
    same_base = 0
    diff_base = 0
    while x < len(seq1):
        if seq1[x] in ["N", "-", "n"] or seq2[x] in ["N","-", "n"]:
            x += 1
            continue
        seq1_bases = [seq1[x].upper()]
        if in_ncar_alt_dic.get(x, None) != None:
            seq1_bases.append(str(in_ncar_alt_dic[x]).upper())
        if seq2[x].upper() in seq1_bases:
            same_base += 1
        else:
            diff_base += 1
        x += 1
    return same_base, diff_base
        
        
        
def ncar_mk_test(inspecies, outspecies, ortho_dic, align_dir, anc_dir, basedir, num_threads, min_taxa, ingff, ingenome, coding_gff, pickle_dir, coding_ortho_dic, coding_anc_dir, coding_align_dir):
    if not os.path.exists("%s/ncar_mk_tests/" % (basedir)):
        os.mkdir("%s/ncar_mk_tests/" % (basedir))
    get_species_data(inspecies, num_threads, basedir, ingenome, coding_gff, pickle_dir)
    in_gene_dic = read_species_pickle(inspecies, pickle_dir, "gene")
    print "%s gene pickle read" % inspecies
    get_species_ncar_data(inspecies, num_threads, basedir, ingff, ingenome, ortho_dic, in_gene_dic, pickle_dir)
#    in_ncar_dic = read_species_pickle(inspecies, pickle_dir, "cnee")
#    in_ncar_dic = read_species_pickle(inspecies, pickle_dir, "ncar")
    in_ncar_dic = read_species_pickle(inspecies, pickle_dir, "cactusncar")
    print "%s ncar pickle read" % inspecies
    outfile = open("%s/ncar_mk_tests/%s_%s_mk.txt" % (basedir, inspecies, outspecies), 'w')
    outfile.write("NCAR\t#poly\t#fixed\t#neib_4d_poly\t#neib_4d_fixed\tncar_sites\t#neib_4d\tNCAR_N\talpha\todds_rat\tpval\n")
    neighbor_file = open("%s/ncar_mk_tests/%s_%s_mk_neighbors.txt" % (basedir, inspecies, outspecies), 'w')
    neighbor_file.write("NCAR\tOG\tfix_fourfold\tin_poly_fourf\tin_potent_fourfold\talign_potent_fourf\n")
    neighbor_fasta = open("%s/ncar_mk_tests/%s_%s_mk_neighbors.fasta" % (basedir, inspecies, outspecies), 'w')
    ncar_fasta = open("%s/ncar_mk_tests/%s_%s_mk_ncar.fasta" % (basedir, inspecies, outspecies), 'w')
    snipre_table = open("%s/ncar_mk_tests/%s_%s_mk_snipre.txt" % (basedir, inspecies, outspecies), 'w')
    snipre_table.write("geneID\tPR\tFR\tPS\tFS\tTsil\tTrepl\tnout\tnpop\n")

    og_map_dic = make_og_gene_map(coding_ortho_dic)
    not_enough_neighbors = 0
    for ncar_num, species_gene_dic in ortho_dic.items():
#        print ncar_num
#        print species_gene_dic
        if len(species_gene_dic.keys()) < min_taxa:
            continue
        inspecies_gene = ortho_dic[ncar_num][inspecies]
#        if "NMEL_chr_1:" not in inspecies_gene:# and "LVIE_chr_4" not in inspecies_gene:
#            continue
        print inspecies_gene
        
        outspecies_gene = ortho_dic[ncar_num][outspecies]
        inseq, outseq = get_fsa_aligned_seqs(align_dir, ncar_num, inspecies, outspecies, "ncar")
        node_index_dic, node_seqs = read_ancestral_seq("%s/ncar_%s_working/rst" % (anc_dir, ncar_num))
        seq_dic = get_fsa_aligned_dic(align_dir, ncar_num, "ncar")
        nearest_anc = get_nearest_anc(node_index_dic, seq_dic)
        inseq = seq_dic[inspecies].upper()
        if nearest_anc[inspecies] == -1:
            continue
        outseq = node_seqs[nearest_anc[inspecies]].upper()
        ncar_fasta.write(">%s_%s\n%s\n>%s_%s\n%s\n" % (ncar_num, inspecies, inseq, ncar_num, outspecies, outseq))
        if prop_shared_sequence(inseq, outseq) < 0.9:
            continue
        in_ncar = in_ncar_dic[inspecies_gene]
        same_base, diff_base = ncar_pairwise_diff_count(inspecies, inseq, outseq, in_ncar)
        in_poly = in_ncar.poly_count
        neighbor_dic = {}
        for neighbor_name in in_ncar.neighbors:
            if in_gene_dic.get(neighbor_name, None) == None:
#                print "%s not in gene dictionary" % neighbor_name
                continue
            if og_map_dic.get(neighbor_name, None) == None:
#                print "%s not in orthology dictionary" % neighbor_name
                continue
            neighbor_og = og_map_dic[neighbor_name]
            if outspecies == "DNOV":
                coding_outspecies = "Dnov"
            else:
                coding_outspecies = outspecies
            if check_og_complete(neighbor_og, inspecies, coding_outspecies, coding_ortho_dic):
                if coding_outspecies == "Dnov":
                    outneighbor_name = coding_ortho_dic[neighbor_og][coding_outspecies][0]
                else:
                    outneighbor_name = coding_ortho_dic[neighbor_og][coding_outspecies][0][:-3]
                neighbor_dic[neighbor_og] = in_gene_dic[neighbor_name]
        if len(neighbor_dic) < 4:
            not_enough_neighbors += 1
            continue

        neighbor_stats_list = []
        direct_stats_list = []
        av_ns = [in_ncar.average_n]
        for neighbor_og, neighbor in neighbor_dic.items():
#            print neighbor_og
            av_ns.append(neighbor.cds_average_n)
            node_index_dic, node_seqs = read_ancestral_seq("%s/og_%s_working/rst" % (coding_anc_dir, neighbor_og))
            seq_dic = get_fsa_aligned_dic(coding_align_dir, neighbor_og, "gene")
            nearest_anc = get_nearest_anc(node_index_dic, seq_dic)
#            print nearest_anc
            outseq = node_seqs[nearest_anc[inspecies]]
            inseq = seq_dic[inspecies]
            neighbor_stats = ingene_fourfold_degenerate_string([neighbor, inspecies, inseq, outseq])
            neighbor_fasta.write(">%s_%s_%s\n%s\n>%s_%s_%s\n%s\n" % (ncar_num, neighbor_og, inspecies, inseq, ncar_num, neighbor_og, outspecies, outseq))
            neighbor_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (ncar_num, neighbor_og, neighbor_stats[0], neighbor_stats[1], neighbor_stats[2], neighbor_stats[3]))
            neighbor_stats_list.append(neighbor_stats)
        in_poly_fourf, fix_fourfold, in_potent_fourfold, align_potent_fourf = sum_ingene_fourfs(neighbor_stats_list)
#        print diff_base
#        print in_poly_fourf
        if diff_base == 0 or in_poly_fourf == 0 or fix_fourfold == 0:
            alpha = "NA"
            neutrality_index = "NA"
            odds_rat = "NA"
            pval = "NA"
        else:
            alpha = 1 - ((fix_fourfold*in_poly) / (diff_base*in_poly_fourf))
            neutrality_index = (in_poly / in_poly_fourf) / (diff_base / fix_fourfold) #the neutrality index is the same as the odds ratio
        odds_rat, pval = scipy.stats.fisher_exact([[fix_fourfold, in_poly_fourf],[diff_base, in_poly]])
        average_n = round(numpy.mean(av_ns))
        outfile.write("NCAR_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ncar_num, in_poly, diff_base, in_poly_fourf, fix_fourfold, in_ncar.potential_poly_sites(), align_potent_fourf, in_ncar.average_n, alpha, odds_rat, pval))
        if in_ncar.average_n > 0:
            snipre_table.write("NCAR_%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\t%s\n" % (ncar_num, in_poly, diff_base, in_poly_fourf, fix_fourfold, in_ncar.potential_poly_sites(), align_potent_fourf, in_ncar.average_n * 2))
    print "not enough neighbors: %s" % not_enough_neighbors
    outfile.close()

def mk_test(inspecies, outspecies, ortho_dic, align_dir, anc_dir, basedir, num_threads, min_taxa, ingenome, coding_gff, pickle_dir):
    if not os.path.exists("%s/coding_mk_tests/" % (basedir)):
        os.mkdir("%s/coding_mk_tests/" % (basedir))
    get_species_data(inspecies, num_threads, basedir, ingenome, coding_gff, pickle_dir)
#    sys.exit()
####    get_species_data(outspecies, num_threads)
    in_gene_dic = read_species_pickle(inspecies, pickle_dir, "gene")
    print "%s gene pickle read" % inspecies
###    out_gene_dic = read_species_pickle(outspecies)
    outfile = open("%s/coding_mk_tests/%s_%s_mk.txt" % (basedir, inspecies, outspecies), 'w')
    outfile.write("OG\t#pn\t#dn\t#ps\t#ds\tpotential_nsyn\tpotential_syn\tmean_N\talpha\todds_rat\tpval\n")
    snipre_table = open("%s/coding_mk_tests/%s_%s_mk_snipre.txt" % (basedir, inspecies, outspecies), 'w')
    snipre_table.write("geneID\tPR\tFR\tPS\tFS\tTsil\tTrepl\tnout\tnpop\n")
    pinpis_table = open("%s/coding_mk_tests/%s_pinpis.txt" % (basedir, inspecies), 'w')
    prop_shared = 0
    starting = 0
    nearest = 0
    for og_num in ortho_dic.keys():
#        print "OG %s" % og_num
#        if og_num != 7:
#            continue
        if len(ortho_dic[og_num].keys()) < min_taxa:
            continue
        if inspecies not in ortho_dic[og_num] or outspecies not in ortho_dic[og_num]:
            continue
        if ortho_dic[og_num][inspecies][0].count(inspecies) > 1 or ortho_dic[og_num][outspecies][0].count(outspecies) > 1:
            continue
#        if includes_paralogs(ortho_dic[og_num]):
#            continue
        if len(ortho_dic[og_num][inspecies]) == 1 and len(ortho_dic[og_num][outspecies]) == 1:
            starting += 1
            inspecies_gene = ortho_dic[og_num][inspecies][0][:-3]
            outspecies_gene = ortho_dic[og_num][outspecies][0][:-3]
            inseq, outseq = get_fsa_aligned_seqs(align_dir, og_num, inspecies, outspecies, "gene")

            ######Turn this on to do ancestral
            node_index_dic, node_seqs = read_ancestral_seq("%s/og_%s_working/rst" % (anc_dir, og_num))
            seq_dic = get_fsa_aligned_dic(align_dir, og_num, "gene")
            nearest_anc = get_nearest_anc(node_index_dic, seq_dic)
            inseq = seq_dic[inspecies]
            if nearest_anc[inspecies] == -1:
                nearest += 1
                continue
            outseq = node_seqs[nearest_anc[inspecies]]
#            if inseq.count("-") / (len(inseq) * 1.0) > 0.05 or outseq.count("-") / (len(outseq) * 1.0) > 0.05:
#                continue
            if prop_shared_sequence(inseq, outseq) < 0.9:
                prop_shared += 1
                continue

#            if og_num != 2110:
#                continue
            fix_syn, fix_nsyn = count_sub_types(inseq, outseq)
            in_gene = in_gene_dic[inspecies_gene]
#####            out_gene = out_gene_dic[outspecies_gene]

#            fix_syn, fix_nsyn = fixed_sub_types(inspecies, outspecies, og_num, align_dir, inseq, outseq, in_gene)
            in_poly_syn = float(in_gene.syn_count) #in_gene_dic[inseq_name].syn_count
            in_poly_nsyn = float(in_gene.nsyn_count) #in_gene_dic[inseq_name].nsyn_count
            if in_poly_syn == -1 or in_poly_nsyn == -1:
                continue
            in_potent_syn = float(in_gene.potent_syn) #in_gene_dic[inseq_name].potent_syn
            in_potent_nsyn = float(in_gene.potent_nsyn) #in_gene_dic[inseq_name].potent_nsyn
            in_n_size = float(in_gene.average_n) #in_gene_dic[inseq_name].average_n * 2
            if inspecies == "LALB":
                insample = in_n_size
            else:
                insample = in_n_size * 2
#            if (fix_syn + fix_nsyn) / (in_potent_syn + in_poten_nsyn) > 0.05:
#                continue
#            if insample < 5:
#                continue
            pis = in_poly_syn / float(in_potent_syn)
            pin = in_poly_nsyn / float(in_potent_nsyn)
            if pis > 0 and pin > 0:
                pinpis = pin / pis
#                pinpis_table.write("%s\t%s\n" % ("\t".join(str(og_num)), str(pinpis)))
                pinpis_table.write("%s\t%s\n" % (str(og_num), str(pinpis)))
            if fix_nsyn == 0 or in_poly_syn == 0 or fix_syn == 0:
                alpha = "NA"
                neutrality_index = "NA"
                odds_rat = "NA"
                pval = "NA"
            else:
                alpha = 1 - ((fix_syn*in_poly_nsyn) / (fix_nsyn*in_poly_syn))
                neutrality_index = (in_poly_nsyn / in_poly_syn) / (fix_nsyn / fix_syn) #the neutrality index is the same as the odds ratio
            odds_rat, pval = scipy.stats.fisher_exact([[fix_syn, in_poly_syn],[fix_nsyn, in_poly_nsyn]])
            outfile.write("OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, in_poly_nsyn, fix_nsyn, in_poly_syn, fix_syn, in_potent_nsyn, in_potent_syn, in_gene.average_n, alpha, odds_rat, pval))
            snipre_table.write("\t".join([str(og_num), str(in_poly_nsyn), str(fix_nsyn), str(in_poly_syn), str(fix_syn), str(in_potent_syn), str(in_potent_nsyn), "1", str(insample)]) + "\n")
    snipre_table.close()
    pinpis_table.close()
    outfile.close()
    print nearest
    print starting
    print prop_shared

def fixed_sub_types(inspecies, outspecies, og_num, align_dir, inseq, outseq, in_gene):
#    inseq, outseq = get_prank_aligned_seqs(align_dir, og_num, inspecies, outspecies)
    in_gene_alt_dic = in_gene.coding_fixed_align(inseq)
#    out_gene_alt_dic = out_gene.coding_fixed_align(outseq)

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
#            if x+cod_index in out_gene_alt_dic.keys():
#                outcodon_list[cod_index].append(str(out_gene_alt_dic[x+cod_index]).upper())
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

def pairs_coding_div(inspecies, outspecies, ortho_dic, align_dir, basedir, num_threads, min_taxa):
    if not os.path.exists("%s/paired_cds_div/" % (basedir)):
        os.mkdir("%s/paired_cds_div/" % (basedir)) 
    outfile = open("%s/paired_cds_div/%s_%s_cds_div.txt" % (basedir, inspecies, outspecies), 'w')
    total_sames = 0
    total_diffs = 0
    for og_num in ortho_dic.keys():
        if len(ortho_dic[og_num].keys()) < min_taxa:
            continue
        if inspecies not in ortho_dic[og_num] or outspecies not in ortho_dic[og_num]:
            continue
        if ortho_dic[og_num][inspecies][0].count(inspecies) > 1 or ortho_dic[og_num][outspecies][0].count(outspecies) > 1:
            continue
        if not len(ortho_dic[og_num][inspecies]) == 1 and len(ortho_dic[og_num][outspecies]) == 1:
            continue
        inspecies_gene = ortho_dic[og_num][inspecies][0][:-3]
        outspecies_gene = ortho_dic[og_num][outspecies][0][:-3]
        inseq, outseq = get_fsa_aligned_seqs(align_dir, og_num, inspecies, outspecies, "gene")
        seq_dic = get_fsa_aligned_dic(align_dir, og_num, "gene")

        if prop_shared_sequence(inseq, outseq) < 0.9 or prop_shared_sequence(outseq, inseq) < 0.9:
            continue
        inseq = inseq.upper()
        outseq = outseq.upper()
        x = 0
        same_count = 0
        diff_count = 0
        
        while x < len(inseq):
            if inseq[x] in ["-", "N"] or outseq[x] in ["-", "N"]:
                x += 1
                continue
            if inseq[x] == outseq[x]:
                same_count += 1
                total_sames += 1
            elif inseq[x] != outseq[x]:
                diff_count += 1
                total_diffs += 1
            x += 1
        outfile.write("OG_%s\t%s\t%s\t%s\n" % (og_num, same_count, diff_count, float(diff_count) / (same_count + diff_count)))
    outfile.close()
    print "%s\t%s\t%s\t%s\t%s" % (inspecies, outspecies, total_sames, total_diffs, float(total_diffs) / (total_sames + total_diffs))

def gather_fixed_v_shared(inspecies, outspecies, ortho_dic, align_dir, basedir, num_threads, min_taxa, pickle_dir):
    if not os.path.exists("%s/fixed_v_shared/" % (basedir)):
        os.mkdir("%s/fixed_v_shared/" % (basedir))

    in_gene_dic = read_species_pickle(inspecies, pickle_dir, "gene")
    print "%s gene pickle read" % inspecies
    out_gene_dic = read_species_pickle(outspecies, pickle_dir.replace(inspecies, outspecies), "gene")
    print "%s gene pickle read" % outspecies
    outfile = open("%s/fixed_v_shared/%s_%s_polys.txt" % (basedir, inspecies, outspecies), 'w')
    outfile.write("OG\tps\tpn\tshare_s\tshare_n\n")
    prop_shared = 0
    starting = 0
    nearest = 0
    total_syn = 0
    total_nsyn = 0
    total_shared_syn = 0
    total_shared_nsyn = 0
    for og_num in ortho_dic.keys():
#        print "OG %s" % og_num
#        if og_num != 7:
#            continue
        if len(ortho_dic[og_num].keys()) < min_taxa:
            continue
        if inspecies not in ortho_dic[og_num] or outspecies not in ortho_dic[og_num]:
            continue
        if ortho_dic[og_num][inspecies][0].count(inspecies) > 1 or ortho_dic[og_num][outspecies][0].count(outspecies) > 1:
            continue
#        if includes_paralogs(ortho_dic[og_num]):
#            continue
        if len(ortho_dic[og_num][inspecies]) == 1 and len(ortho_dic[og_num][outspecies]) == 1:
            starting += 1
            inspecies_gene = ortho_dic[og_num][inspecies][0][:-3]
            outspecies_gene = ortho_dic[og_num][outspecies][0][:-3]
            inseq, outseq = get_fsa_aligned_seqs(align_dir, og_num, inspecies, outspecies, "gene")

            seq_dic = get_fsa_aligned_dic(align_dir, og_num, "gene")
            inseq = seq_dic[inspecies]

            if prop_shared_sequence(inseq, outseq) < 0.9 or prop_shared_sequence(outseq, inseq) < 0.9:
                prop_shared += 1
                continue
            in_gene = in_gene_dic[inspecies_gene]
            out_gene = out_gene_dic[outspecies_gene]

            if in_gene.cds_average_n < 4 or out_gene.cds_average_n < 4:
                continue
#            fix_syn, fix_nsyn = count_sub_types(inseq, outseq)
            poly_syn, poly_nsyn, shared_syns, shared_nsyns = fixed_v_shared_subs(inspecies, outspecies, og_num, align_dir, inseq, outseq, in_gene, out_gene)
            total_nsyn += poly_nsyn
            total_syn += poly_syn
            total_shared_syn += shared_syns
            total_shared_nsyn += shared_nsyns

#             in_poly_syn = float(in_gene.syn_count) 
#             in_poly_nsyn = float(in_gene.nsyn_count) 
#             if in_poly_syn == -1 or in_poly_nsyn == -1:
#                 continue
#             in_potent_syn = float(in_gene.potent_syn)
#             in_potent_nsyn = float(in_gene.potent_nsyn)
#             in_n_size = float(in_gene.average_n) 
#             if inspecies == "LALB":
#                 insample = in_n_size
#             else:
#                 insample = in_n_size * 2
# #            if (fix_syn + fix_nsyn) / (in_potent_syn + in_poten_nsyn) > 0.05:
# #                continue
# #            if insample < 5:
# #                continue
#             pis = in_poly_syn / float(in_potent_syn)
#             pin = in_poly_nsyn / float(in_potent_nsyn)
#             if pis > 0 and pin > 0:
#                 pinpis = pin / pis
# #                pinpis_table.write("%s\t%s\n" % ("\t".join(str(og_num)), str(pinpis)))
#                 pinpis_table.write("%s\t%s\n" % (str(og_num),
#                                                 str(pinpis)))
            outfile.write("OG_%s\t%s\t%s\t%s\t%s\n" % (og_num, poly_syn, poly_nsyn, shared_syns, shared_nsyns))#, in_poly_syn, in_poly_nsyn))
#    pinpis_table.close()
    outfile.close()
    print "%s\t%s\t%s\t%s\t%s\t%s" % (inspecies, outspecies, total_syn, total_nsyn, total_shared_syn, total_shared_nsyn)
#    print nearest
#    print starting
#    print prop_shared


def fixed_v_shared_subs(inspecies, outspecies, og_num, align_dir, inseq, outseq, in_gene, out_gene):
    inseq, outseq = get_fsa_aligned_seqs(align_dir, og_num, inspecies, outspecies, "gene")
    in_gene_alt_dic = in_gene.coding_fixed_align(inseq)
    out_gene_alt_dic = out_gene.coding_fixed_align(outseq)
    temp_dic = {}
    for k, v in in_gene_alt_dic.items():
        if len(v) == 1:
            temp_dic[k] = v
    in_gene_alt_dic = temp_dic
    temp_dic = {}
    for k, v in out_gene_alt_dic.items():
        if len(v) == 1:
            temp_dic[k] = v
    out_gene_alt_dic = temp_dic
#    print in_gene_alt_dic
#    print out_gene_alt_dic
    nsyn_count = 0
    syn_count = 0
    same_codons = 0
    x = 0
    shared_polys = 0
    shared_syns = 0
    shared_nsyns = 0
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
        has_variant = False
        if in_gene.strand == 1:
            ingap_count = inseq[0:x].count("-")
        else:
            ingap_count = -1*inseq[x:].count("-")
        if out_gene.strand == 1:
            outgap_count = outseq[0:x].count("-")
        else:
            outgap_count = -1*outseq[x:].count("-")
        for cod_index in range(3):
            incodon_list.append([inseq[x+cod_index].upper()])
            outcodon_list.append([outseq[x+cod_index].upper()])
            if x+cod_index-ingap_count in in_gene_alt_dic.keys():
                incodon_list[cod_index].append(str(in_gene_alt_dic[x+cod_index-ingap_count]).upper())
                has_variant = True
            if x+cod_index-outgap_count in out_gene_alt_dic.keys():
                outcodon_list[cod_index].append(str(out_gene_alt_dic[x+cod_index-outgap_count]).upper())
        if not has_variant:
            x = x+3
            continue
        overlap = 0
        overlap_incodon = list(inseq[x:x+3])
        overlap_outcodon = list(outseq[x:x+3])
        variant_sites = 0
        
        #remove codons with multiple variants
        for cur_site in incodon_list:
            if len(cur_site) > 1:
                variant_sites += 1
        if variant_sites > 1:
            x = x+3
            continue

        variant_sites = 0
        for cur_site in outcodon_list:
            if len(cur_site) > 1:
                variant_sites += 1
        if variant_sites > 1:
            x = x+3
            continue

        #is incodon variant syn or nonsyn
        ref_allele = []
        alt_allele = []
        for cur_site in incodon_list:
            if len(cur_site) > 1:
                ref_allele.append(cur_site[0])
                alt_allele.append(cur_site[1])
            else:
                ref_allele.append(cur_site[0])
                alt_allele.append(cur_site[0])
        ref_allele = "".join(ref_allele)
        alt_allele = "".join(alt_allele)
        is_syn = True
        if str(Seq.Seq(ref_allele).translate()) != str(Seq.Seq(alt_allele).translate()):
            is_syn = False
            
                                   
        if is_syn:
            syn_count += 1
        else:
            nsyn_count += 1
        # print x
        # print is_syn
        # print incodon_list
        # print outcodon_list


        for cur_site in range(len(incodon_list)):
            poly_overlap = 0
            for nuc in incodon_list[cur_site]:
                if nuc in outcodon_list[cur_site]:
                    poly_overlap += 1
                    
            if poly_overlap > 1:
                shared_polys += 1
                if is_syn:
                    shared_syns += 1
                else:
                    shared_nsyns += 1
#                 print in_gene.name
#                 print out_gene.name
# #            if len(incodon_list[cur_site]) > 1 or len(outcodon_list[cur_site]) > 1:

#                 print x
#                 print poly_overlap
#                 print incodon_list
#                 print outcodon_list
                
#                    overlap_incodon[cur_site] = nuc
#                    overlap_outcodon[cur_site] = nuc
                    #if there is overlap at a site then just got to next site
                    #don't need to keep checking current site
#                    break

        incodon = "".join(overlap_incodon)
        outcodon = "".join(overlap_outcodon)
#        print incodon
#        print outcodon
        # if incodon == outcodon:
        #     same_codons += 1
        # else:
        #     num_diffs = 0
        #     for cod_index in range(3):
        #         if incodon[cod_index] != outcodon[cod_index]:
        #             num_diffs += 1
        #     if num_diffs == 1:
        #         inaa = str(Seq.Seq(incodon).translate())
        #         outaa = str(Seq.Seq(outcodon).translate())
        #         if inaa == outaa:
        #             syn_count += 1
        #         else:
        #             nsyn_count += 1
#        print syn_count
#        print nsyn_count
        x = x + 3
#    print in_gene_alt_dic
#    print out_gene_alt_dic

    return syn_count, nsyn_count, shared_syns, shared_nsyns


def remove_stops(seq):
    x = 0
    while x < len(seq):
        if seq[x:x+3] in STOP_CODONS:
            seq = seq[:x] + "NNN" + seq[x+3:]
        x += 3
    return seq

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

def trim_phylo_hyphy(taxa_list, fore_list, orthogroup, outdir, phylogeny_file):
    #trims taxa and adds foreground tags for HYPHY analysis tree
    tree = PhyloTree(phylogeny_file)
    tree.prune(taxa_list)
    tree.unroot()
    tree_str = tree.write(format = 5)
    for leaf in tree.get_leaves():
        cur_tax = leaf.name
        if cur_tax in fore_list:
            tree_str = tree_str.replace(cur_tax, "%s{foreground}" % cur_tax)
    outfile = open("%s/og_%s.tree" % (outdir, orthogroup), 'w')
    outfile.write(tree_str)
    outfile.close()

def trim_phylo_alldaughters_hyphy(taxa_list, fore_list, orthogroup, outdir, phylogeny_file):
    #trims taxa and adds foreground tags for HYPHY analysis tree
    tree = PhyloTree(phylogeny_file)
    tree.prune(taxa_list)
    tree.unroot()
    for node in tree.traverse("postorder"):
        good_daughters = 0
        if node.is_leaf():
            continue
        for daughter in node.get_leaves():
            if daughter.name in fore_list:
                good_daughters += 1
        if good_daughters == len(node.get_leaves()):
            node.name = "{foreground}"
    tree_str = tree.write(format = 3)
    tree_str = tree_str.replace("NoName", "")
    for leaf in tree.get_leaves():
        cur_tax = leaf.name
        if cur_tax in fore_list:
            tree_str = tree_str.replace(cur_tax, "%s{foreground}" % cur_tax)
    
    outfile = open("%s/og_%s.tree" % (outdir, orthogroup), 'w')
    outfile.write(tree_str)
    outfile.close()


def trim_phylo_ancestral_hyphy(taxa_list, orthogroup, outdir, phylogeny_file):
    #trims taxa and adds foreground tags for HYPHY analysis tree
    tree = PhyloTree(phylogeny_file, format = 1)
    tree.prune(taxa_list, preserve_branch_length = True)
    tree.unroot()
    tree_str = tree.write(format = 3)
    tree_str = tree_str.replace("NoName", "")
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

def aligned_og_completeness(og_list, align_dir, min_taxa):
    new_og_list = []
    for og in og_list:
        reader = SeqIO.parse("%s/og_cds_%s.afa" % (align_dir, og), format = 'fasta')
        seq_dic = {}
        for rec in reader:
            seq_dic[rec.id] = str(rec.seq)
        if len(seq_dic) >= min_taxa:
            new_og_list.append(og)
    return new_og_list

def prep_paml_files(orthogroup, indir, outdir, foreground, phylogeny_file, test_type, min_taxa, use_gblocks, exclude_taxa):
    #formats fasta and tree files for PAML analysis
    print orthogroup
    print foreground
    tree_prep = True
    terminals = False
    fore_list = []
    SOCIAL = ["HLIG","LMAL", "LMAR", "LPAU", "AAUR", "LZEP"]
    REV_SOLITARY = ["LLEU", "LOEN", "LVIE", "LFIG", "APUR", "HQUA"]
    if foreground in ["HLIG","LMAL", "LMAR", "LPAU", "AAUR", "LZEP", "LLEU", "LOEN", "LVIE", "LFIG", "APUR", "FNIG", "BNIG", "HQUA", "HRUB", "AVIR"]:
        fore_list = [foreground]
        terminals = True
    if foreground == "social":
        fore_list = SOCIAL
    elif foreground == "solitary":
        fore_list = REV_SOLITARY
    elif foreground == "model_d":
        tree_prep = False
    elif foreground == "lasioglossum" or foreground == "augochlorine" or foreground == "halictus" or foreground == "lasihali" or foreground == "lasiaugo":
        tree_prep = False
    elif "clade" in foreground:
        tree_prep = False
    elif isinstance(foreground, list):
        fore_list = foreground
    if use_gblocks:
        if foreground == "ancestral":
            reader = SeqIO.parse("%s/og_cds_%s.afa-gb" % ( indir, orthogroup), format = 'fasta')
        elif foreground == "yn":
            reader = SeqIO.parse("%s/og_cds_%s.afa-gb" % ( indir, orthogroup), format = 'fasta')
            tree_prep = False
        elif foreground == "aaml_blengths":
            reader = SeqIO.parse("%s/og_cds_%s.afa-gb" % ( indir, orthogroup), format = 'fasta')
        else:
            reader = SeqIO.parse("%s/og_cds_%s.afa-gb" % ( indir, orthogroup), format = 'fasta')
    else:
        if test_type == "ancestral":
            if foreground == "ncar":
                reader = SeqIO.parse("%s/ncar_%s.afa" % ( indir, orthogroup), format = 'fasta')
            else:
                reader = SeqIO.parse("%s/og_cds_%s.afa" % ( indir, orthogroup), format = 'fasta')
        elif foreground == "yn":
            reader = SeqIO.parse("%s/og_cds_%s.afa" % ( indir, orthogroup), format = 'fasta')
            tree_prep = False
        elif foreground == "aaml_blengths":
            reader = SeqIO.parse("%s/og_cds_%s.afa" % ( indir, orthogroup), format = 'fasta')
        else:
            reader = SeqIO.parse("%s/og_cds_%s.afa" % ( indir, orthogroup), format = 'fasta')

    seq_dic = {}
    taxa_list = []
    for rec in reader:
        if rec.id[0:4] not in exclude_taxa:
 #           if rec.id[0:4] == "LMAR":
            if foreground == "ncar":
                seq_dic[rec.id] = str(rec.seq).upper()
            else:
                seq_dic[rec.id] = mask_selenocysteine(str(rec.seq).upper())
            taxa_list.append(rec.id[0:4])
#    if len(seq_dic) < min_taxa:
#        print orthogroup
#        return "too short"
    if terminals:
        if foreground not in taxa_list:
            return False
    if foreground == "ncar":
        outfile = open("%s/ncar_%s.afa" % (outdir, orthogroup), 'w')
    else:
        outfile = open("%s/og_cds_%s.afa" % (outdir, orthogroup), 'w')
        for species, sequence in seq_dic.items():
            if len(sequence) % 3 > 0:
                outfile.close()
                print "Not divisible by 3 %s" % orthogroup
                return False
            if str(Seq.Seq(sequence.replace("-", "N")).translate()).count("*") > 1:
                outfile.close()
                print "Too many stop codons %s" % orthogroup
                return False
        
            
    if foreground == "aaml_blengths":
        outfile.write("%s %s\n" % (len(seq_dic), len(rec.seq) / 3))
    elif test_type not in ["RELAX", "aBSREL"]:
        outfile.write("%s %s\n" % (len(seq_dic), len(rec.seq)))
    for species, sequence in seq_dic.items():
        if foreground == "aaml_blengths":
            outfile.write("%s\n%s\n" % (species[0:4], str(Seq.Seq(sequence.replace("-", "N")).translate())))
        elif foreground == "free" or foreground == "yn":
            if len(sequence) % 3 > 0:
                print species
                print str(Seq.Seq(sequence.replace("-", "N")).translate())
            if str(Seq.Seq(sequence.replace("-", "N")).translate()).count("*") > 0:
                outfile.close()
                print "Too many stop codons in %s %s" % (species, orthogroup)
                return False
            else:
                outfile.write("%s\n%s\n" % (species[0:4], sequence))
        elif test_type in ["RELAX", "aBSREL"]:
            outfile.write(">%s\n%s\n" % (species[0:4], sequence))
        else:
            outfile.write("%s\n%s\n" % (species[0:4], sequence))
    outfile.close()
    if tree_prep:
        if test_type == "RELAX":
            print foreground
            if foreground == "INTREE":
#                shutil.copyfile(phylogeny_file, "%s/og_%s.tree" % (outdir, orthogroup))
                trim_phylo_ancestral_hyphy(taxa_list, orthogroup, outdir, phylogeny_file)
            elif foreground[0] == "DAUGHTERS":
                print "DAUGHT"
                trim_phylo_alldaughters_hyphy(taxa_list, fore_list[1:], orthogroup, outdir, phylogeny_file)
            else:
                trim_phylo_hyphy(taxa_list, fore_list, orthogroup, outdir, phylogeny_file)
        else:
            trim_phylo(taxa_list, fore_list, orthogroup, outdir, phylogeny_file)
    elif "lasihali" in foreground:
        lasihali_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file)
    elif "augochlorine" in foreground:
        augo_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file)
    elif "halictus" in foreground:
        hali_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file)
    elif "lasioglossum" in foreground:
        lasi_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file)
    elif "lasiaugo" in foreground:
        lasiaugo_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file)

#        shutil.copy("/Genomics/kocherlab/berubin/annotation/orthology/halictids_rough/lasioglossum.tree", "%s/og_%s.tree" % (outdir, orthogroup))
    elif foreground != "yn":
        shutil.copy("/Genomics/kocherlab/berubin/annotation/orthology/model_d.tree", "%s/og_%s.tree" % (outdir, orthogroup))
#    for seq1 in seq_dic.values():
#        for seq2 in seq_dic.values():
#            if align_len(seq1,seq2) < 100:
#                return False
    return True

def mask_selenocysteine(inseq):
    index = 0
    new_seq = []
    while index < len(inseq):
#        print Seq.Seq(inseq[index:index+3].replace("-","N")).translate()
#        print Seq.Seq(inseq[index:index+3].replace("-","N"))
        if inseq[index:index+3] in STOP_CODONS:
            new_seq.append("NNN")
        elif "N" in inseq[index:index+3]:
            new_seq.append("NNN")
        elif str(Seq.Seq(inseq[index:index+3].replace("-","N")).translate()) in ["J", "Z", "B"]:

            new_seq.append("NNN")
        else:
            new_seq.append(inseq[index:index+3])
        index += 3
    return "".join(new_seq)

def terminal_test_overlap(indir, prefix, test_type, target_taxa, soc_list, sol_list, pairs_list):
    sig_dic = {}
    og_list = []
    for taxon in target_taxa:
        reader = open("%s/%s_%s_%s.lrt" % (indir, prefix, taxon, test_type), 'rU')
        sig_dic[taxon] = {}
        for line in reader:
            if line.startswith("OG"):
                continue
            cur_line = line.strip().split()
            cur_og = cur_line[0]
            cur_p = float(cur_line[1])
            sig_dic[taxon][cur_og] = cur_p
            if cur_p < 0.05:
                og_list.append(cur_og)
    og_list = list(set(og_list))
    outfile = open("%s/%s_%s_soc_sol_counts.txt" % (indir, prefix, test_type), 'w')
    outfile.write("OG\tsocs\tsols\tsoc_pairs\tsol_pairs\t%s\n" % "\t".join(target_taxa))

    for og in og_list:
#        if og != "3499":
#            continue
        soc_count = 0
        sol_count = 0
        soc_pairs_count = 0
        sol_pairs_count = 0
        for soc in soc_list:
            if og in sig_dic[soc].keys():
                if sig_dic[soc][og] < 0.05:
                    soc_count += 1
        for sol in sol_list:
            if og in sig_dic[sol].keys():
                if sig_dic[sol][og] < 0.05:
                    sol_count += 1

        for pair in pairs_list:
#            print og
#            print pair
#            print sig_dic[pair[0]][og]

            if og in sig_dic[pair[0]].keys() and og in sig_dic[pair[1]].keys():
                if sig_dic[pair[0]][og] < 0.05 and sig_dic[pair[1]][og] > 0.05:
                    soc_pairs_count += 1
#            if og in sig_dic[pair[1]].keys() and og in sig_dic[pair[0]].keys():           
                elif sig_dic[pair[1]][og] < 0.05 and sig_dic[pair[0]][og] > 0.05:
                    sol_pairs_count += 1
        out_ps = []
        for taxon in target_taxa:
            if og in sig_dic[taxon].keys():
                out_ps.append(str(round(sig_dic[taxon][og], 5)))
            else:
                out_ps.append("NA")
        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (og, soc_count, sol_count, soc_pairs_count, sol_pairs_count, "\t".join(out_ps)))
    outfile.close()
        
    

def test_lrt_branch(indir, outpath, database_file, ortho_dic, go_dir):
    #performs all lrts on the files in a given directory
    #also performs multiple test correction to return adjusted p-value
    p_dic = {}
    p_list = []
    og_list = []
    pos_change = []
    neg_change = []
    for og_file in glob("%s/og_*.nul" % (indir)):
        cur_og = int(og_file.split("og_")[1].split(".nul")[0])
        pval = lrt("%s/og_%s.alt" % (indir, cur_og), "%s/og_%s.nul" % (indir, cur_og))
        if pval == "bad_len":
            continue
#        print "%s: %s" % (cur_og, pval)
        p_list.append(pval)
        og_list.append(cur_og)
        reader = open("%s/og_%s.alt" % (indir, cur_og), 'rU')
        for line in reader:
            if line.startswith("w (dN/dS) for branches:"):
                cur_line = line.strip().split()
                backw = float(cur_line[-2])
                forew = float(cur_line[-1])
                if backw < forew:
                    pos_change.append(cur_og)
                else:
                    neg_change.append(cur_og)
    pval_corr = smm.multipletests(p_list, alpha = 0.05, method = 'fdr_i')[1]
    fastfile = open("%s_faster.lrt" % outpath, 'w')
    slowfile = open("%s_slower.lrt" % outpath, 'w')
    outfile = open("%s.lrt" % outpath, 'w')
    outfile.write("OG\tpval_corr\n")
    fastfile.write("OG\tpval_corr\n")
    slowfile.write("OG\tpval_corr\n")
    for x in range(len(pval_corr)):
        p_dic[og_list[x]] = pval_corr[x]

    sorted_pvals = sorted(p_dic.items(), key = lambda x: x[1])
    sig_ogs_fast = []
    sig_ogs_slow = []
    for pval in sorted_pvals:
        if float(pval[1]) < 0.05:
            if pval[0] in pos_change:
                sig_ogs_fast.append(pval[0])
                fastfile.write("%s\t%s\n" % (pval[0], pval[1]))
            else:
                sig_ogs_slow.append(pval[0])
                slowfile.write("%s\t%s\n" % (pval[0], pval[1]))
        outfile.write("%s\t%s\n" % (pval[0], pval[1]))
    fastfile.close()
    slowfile.close()
    outfile.close()
    run_termfinder(sig_ogs_fast, og_list, database_file, ortho_dic, "%s_faster" % go_dir)
    run_termfinder(sig_ogs_slow, og_list, database_file, ortho_dic, "%s_slower" % go_dir)
    return p_dic
 

def test_lrt(indir, outpath, database_file, ortho_dic, go_dir):
    #performs all lrts on the files in a given directory
    #also performs multiple test correction to return adjusted p-value
    p_dic = {}
    p_list = []
    og_list = []
    for og_file in glob("%s/og_*.nul" % (indir)):
        cur_og = int(og_file.split("og_")[1].split(".nul")[0])
        pval = lrt("%s/og_%s.alt" % (indir, cur_og), "%s/og_%s.nul" % (indir, cur_og))
        if pval == "bad_len":
            continue
#        print "%s: %s" % (cur_og, pval)
        p_list.append(pval)
        og_list.append(cur_og)
    pval_corr = smm.multipletests(p_list, alpha = 0.05, method = 'fdr_i')[1]
    outfile = open("%s.lrt" % outpath, 'w')
    outfile.write("OG\tpval_corr\n")
    for x in range(len(pval_corr)):
        p_dic[og_list[x]] = pval_corr[x]

    sorted_pvals = sorted(p_dic.items(), key = lambda x: x[1])
    sig_ogs = []
    for pval in sorted_pvals:
        if float(pval[1]) < 0.05:
            sig_ogs.append(pval[0])
        outfile.write("%s\t%s\n" % (pval[0], pval[1]))
    outfile.close()
    run_termfinder(sig_ogs, og_list, database_file, ortho_dic, go_dir)
    return p_dic

def lrt(alt_file, null_file):
    #lrt for PAML tests
    reader = open(alt_file, 'rU')
    firstline = True
    for line in reader:
        if firstline:
            seq_len = int(line.strip().split()[1])
            firstline = False
            if seq_len < 300:
                return "bad_len"
        if line.startswith("lnL"):
            alt_likely = float(line.split()[4])
    reader = open(null_file, 'rU')
    for line in reader:
        if line.startswith("lnL"):
            nul_likely = float(line.split()[4])
    p = chisqprob(2*(alt_likely - nul_likely), 1)
    return p

def compile_gene_trees(og_list, indir, ref_tree, outdir):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    og_bl_dic = {}
    ref_tr = PhyloTree(ref_tree)
    bad_count = 0
    for og in og_list:
        print og
        if not os.path.exists("%s/RAxML_bestTree.og_%s.tree" % (indir, og)):
            continue
        tree = PhyloTree("%s/RAxML_bestTree.og_%s.tree" % (indir, og))
#        tree.set_outgroup("AMEL")
        cur_bl_dic = {}
        cur_leaves = []
        for leaf in tree:
            if leaf.name == "BTER":
                continue
            cur_leaves.append(leaf.name)
        tree.prune(cur_leaves, preserve_branch_length = True)
        for leaf in tree:
            cur_bl_dic[leaf.name] = leaf.dist
#        tree.prune(cur_bl_dic.keys(), preserve_branch_length = True)
        time_tree = PhyloTree(ref_tree)
        time_tree.prune(cur_bl_dic.keys(), preserve_branch_length = True)
        time_dic = standardize_to_time(cur_bl_dic, time_tree)
#        time_dic = cur_bl_dic


        if tree.robinson_foulds(ref_tr)[0] > 0:
            print og
            print tree
            print tree.robinson_foulds(ref_tr)
            bad_count += 1
            continue
        for leaf in tree:
            if leaf.name == "BTER":
                continue
            if leaf.name not in og_bl_dic.keys():
                
                og_bl_dic[leaf.name] = []
            og_bl_dic[leaf.name].append(time_dic[leaf.name])
#            og_bl_dic[leaf.name].append(leaf.dist)
    outfile = open("%s/bls.txt" % (outdir), 'w')
    for species in og_bl_dic.keys():
        outfile.write("%s\t%s\t%s\n" % (species, numpy.median(og_bl_dic[species]), numpy.std(og_bl_dic[species])))
    print bad_count
    outfile.close()

def standardize_to_time(target_dic, ref_tree):
    standard_dic = {}
    for k, v in target_dic.items():
        if k in ["Dnov", "Nmel"]:
            continue
        standard_dic[k] = v / ref_tree.get_leaves_by_name(k)[0].dist
    return standard_dic

def read_frees(indir, outdir, database_file, go_dir, get_dn_ds, time_tree, og_list):
    #reads free ratios files and gets dn/ds ratios
    #can be easily extended to get dn and ds but those are low quality
    soc_sol_pairs = {"LMAR":"LFIG", "LZEP":"LVIE", "LPAU":"LOEN", "AAUR":"APUR"}
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    og_dnds_dic = {}
    og_ds_dic = {}
    og_dn_dic = {}
    species_list = []
    for og_file in glob("%s/og_*.alt" % (indir)):        
        cur_og = int(og_file.split("og_")[1].split(".alt")[0])
        if cur_og not in og_list:
            continue
        reader = open(og_file, 'rU')
        dnds_tree = False
        ds_tree = False
        dn_tree = False
        ds_dic = {}
        dn_dic = {}
        dnds_dic = {}
        first_line = True
        for line in reader:
            if first_line:
                seq_len = int(line.strip().split()[1])
                first_line = False
            if dnds_tree:
                dnds = PhyloTree(line.strip().replace("#", ":"))
                for leaf in dnds:
                    if leaf.dist < 4 and leaf.dist > 0.0001:
                        if ds_dic[leaf.name] > 0.001 and ds_dic[leaf.name] < 10:
                            if dn_dic[leaf.name] > 0.001 and dn_dic[leaf.name] < 10:
                                dnds_dic[leaf.name] = leaf.dist
                    if leaf.name not in species_list:
                        species_list.append(leaf.name)
                dnds_tree = False
            if ds_tree:
                ds = PhyloTree(line.strip())
                for leaf in ds:
                    ds_dic[leaf.name] = leaf.dist
                ds_tree = False
            if dn_tree:
                dn = PhyloTree(line.strip())
                for leaf in dn:
                    dn_dic[leaf.name] = leaf.dist
                dn_tree = False
            if line.strip() == "dS tree:":
                ds_tree = True
                continue
            if line.strip() == "dN tree:":
                dn_tree = True
                continue
            if line.strip() == "w ratios as labels for TreeView:":
                dnds_tree = True
                continue
        if seq_len < 300:
            continue
        og_dnds_dic[cur_og] = dnds_dic
        if get_dn_ds:
            taxa_pres = dn_dic.keys()
            if len(taxa_pres) == 0:
                continue
            if "Dnov" in taxa_pres:
                taxa_pres.remove("Dnov")
            if "Nmel" in taxa_pres:
                taxa_pres.remove("Nmel")
            tree = PhyloTree(time_tree)
            tree.prune(taxa_pres, preserve_branch_length = True)
            dn_time_dic = standardize_to_time(dn_dic, tree)
            taxa_pres = ds_dic.keys()
            if "Dnov" in taxa_pres:
                taxa_pres.remove("Dnov")
            if "Nmel" in taxa_pres:
                taxa_pres.remove("Nmel")
            tree = PhyloTree(time_tree)
            tree.prune(taxa_pres, preserve_branch_length = True)
            ds_time_dic = standardize_to_time(ds_dic, tree)
            og_ds_dic[cur_og] = ds_time_dic
            og_dn_dic[cur_og] = dn_time_dic
            og_ds_dic[cur_og] = ds_dic
            og_dn_dic[cur_og] = dn_dic
    sum_file = open("%s/dnds_summary.txt" % outdir, 'w')
    for species in species_list:
        outfile = open("%s/%s_free_dnds.txt" % (outdir, species), 'w')
        dnds_list = []
        for og in og_dnds_dic.keys():
            if species in og_dnds_dic[og].keys():
                outfile.write("%s\t%s\n" % (og, og_dnds_dic[og][species]))
                dnds_list.append(og_dnds_dic[og][species])
            else:
                outfile.write("%s\tNA\n" % (og))
            
        outfile.close()
        sum_file.write("%s\t%s\t%s\n" % (species, numpy.mean(dnds_list), numpy.var(dnds_list)))
    if get_dn_ds:
        sum_file = open("%s/dn_ds_summary.txt" % outdir, 'w')
        for species in species_list:
            dn_list = []
            ds_list = []
            dnfile = open("%s/%s_free_dn.txt" % (outdir, species), 'w')
            dsfile = open("%s/%s_free_ds.txt" % (outdir, species), 'w')
            for og in og_dn_dic.keys():
                if species in og_dn_dic[og].keys():
                    if og_dn_dic[og][species] < 1:
                        dnfile.write("%s\t%s\n" % (og, og_dn_dic[og][species]))
                        dn_list.append(og_dn_dic[og][species])
                    else:
                        dnfile.write("%s\tNA\n" % (og))
                    if og_ds_dic[og][species] < 1:
                        dsfile.write("%s\t%s\n" % (og, og_ds_dic[og][species]))
                        ds_list.append(og_ds_dic[og][species])
                    else:
                        dsfile.write("%s\tNA\n" % (og))
                else:
                    dnfile.write("%s\tNA\n" % (og))
                    dsfile.write("%s\tNA\n" % (og))
            sum_file.write("%s\t%s\t%s\t%s\t%s\n" % (species, numpy.mean(dn_list), numpy.var(dn_list), numpy.mean(ds_list), numpy.var(ds_list)))
            dnfile.close()
            dsfile.close()
        sum_file.close()

    soc_larger_file = open("%s/soc_larger.txt" % outdir, 'w')
    sol_larger_file = open("%s/sol_larger.txt" % outdir, 'w')
    soc_list = []
    sol_list = []
    background_ogs = []
    for og in og_dnds_dic.keys():
        if og not in og_list:
            continue
        sol_larger = True
        soc_larger = True
        background_og = True
        for soc, sol in soc_sol_pairs.items():
            if soc in og_dnds_dic[og].keys() and sol in og_dnds_dic[og].keys():
                if og_dnds_dic[og][soc] <= og_dnds_dic[og][sol]:
                    soc_larger = False
                if og_dnds_dic[og][sol] <= og_dnds_dic[og][soc]:
                    sol_larger = False
            else:
                soc_larger = False
                sol_larger = False
                background_og = False
        if soc_larger:
            soc_larger_file.write("%s\n" % og)
            soc_list.append(og)
        if sol_larger:
            sol_larger_file.write("%s\n" % og)
            sol_list.append(og)
        if background_og:
            background_ogs.append(og)
    print len(background_ogs)
    soc_larger_file.close()
    sol_larger_file.close()
#    run_termfinder(soc_list, background_ogs, database_file, ortho_dic, "%s_soc" % go_dir)
#    run_termfinder(sol_list, background_ogs, database_file, ortho_dic, "%s_sol" % go_dir)
    sum_file.close()

def read_aaml_phylos(og_list, indir, outdir, outfile, min_taxa):#, full_tree):
    #reads free ratios files and gets dn/ds ratios
    #can be easily extended to get dn and ds but those are low quality
    soc_sol_pairs = {"LMAR":"LFIG", "LZEP":"LVIE", "LPAU":"LOEN", "AAUR":"APUR"}
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    outfile = open("%s/%s" % (outdir, outfile), 'w')
    og_phylo_dic = {}
    for cur_og in og_list:
        if not os.path.exists("%s/og_%s.alt" % (indir, cur_og)):
            print "%s/og_%s.alt is missing" % (indir, cur_og)
            continue
        og_file =open("%s/og_%s.alt" % (indir, cur_og), 'rU')
        aa_tree = False
        first_line = True
        line = og_file.readline()
        while True:
            if first_line:
#                print cur_og
#                print line
                seq_len = int(line.strip().split()[1])
                first_line = False
            if aa_tree:
                aa = PhyloTree(line.strip().replace("#", ":"))
#                if len(aa.get_leaves()) < min_taxa:
#                    break
                og_phylo_dic[cur_og] = aa

                outfile.write("OG_%s\t%s\n" % (cur_og, aa.write(format = 5)))
                aa_tree = False
                break
            if line.strip().startswith("tree length = "):
                og_file.readline()
                og_file.readline()
                og_file.readline()
                line = og_file.readline()
                aa_tree = True
                continue
#            if seq_len < 100:
#                break
            line = og_file.readline()
            if not line:
                break
#    marine_test(og_phylo_dic, PhyloTree(full_tree))
#    blen_dists = mean_blengths(og_phylo_dic, PhyloTree(full_tree))
#    for node, dist in blen_dists.items():
#        outfile = open("%s/%s_blens.txt" % (outdir, node), 'w')
#        for blen in dist:
#            outfile.write("%s\n" % (blen))
#        outfile.close()
    outfile.close()
    return og_phylo_dic
        
def aaml_time_phylos(og_list, indir, outdir, calib_tree, fore_list):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    outfile = open("%s/aaml_time_calibrated.txt" % outdir, 'w')
    og_phylo_dic = read_aaml_phylos(og_list, indir, outdir)
    time_tree = PhyloTree(calib_tree)
    time_tree.unroot()
    species_list = []
    time_dic = {}
    time_blengths = get_blengths(time_tree, time_tree)
    temp_blengths = {}
    for k, v in time_blengths.items():
        temp_blengths[str(k)] = time_blengths[k]
    time_blengths = temp_blengths
#    print time_blengths
#    for leaf in time_tree.get_leaves():
#        species_list.append(leaf.name)
#        time_dic[leaf.name] = leaf.dist

        
#    outfile.write("OG\t%s\twilcox_p\tinmed\toutmed\n" % ("\t".join(str(time_blengths.keys()))))
    header_written = False
    for og, phylo in og_phylo_dic.items():
        og_blengths = get_blengths(phylo, time_tree)
        temp_blengths = {}
        for k, v in og_blengths.items():
            temp_blengths[str(k)] = og_blengths[k]
        og_blengths = temp_blengths
#        print og_blengths
        aa_dic = {}
        node_names = []
        header_names = []
        for k, v in time_blengths.items():
            aa_dic[k] = og_blengths[k] / time_blengths[k]
            node_names.append(k)
            header_names.append(str(k))
        if not header_written:
            outfile.write("OG\t%s\twilcox_p\tinmed\toutmed\n" % ("\t".join(header_names)))
            header_written = True
        outline = "OG_%s\t" % og
        # aa_dic = {}
        # for leaf in phylo.get_leaves():
        #     aa_dic[leaf.name] = leaf.dist
        # aa_time_dic = standardize_to_time(aa_dic, time_tree)
        
        # for species in species_list:
        #    outline = outline + "\t" + str(aa_time_dic[species])
        for node_name in node_names:
            outline = outline + "\t" + str(round(aa_dic[node_name], 5))
        
        kendall, inmed, outmed = kendall_test(og_blengths, time_blengths, fore_list, [])
        outline = "%s\t%s\t%s\t%s\n" % (outline, kendall, round(inmed,5), round(outmed, 5))
        outfile.write(outline)
    outfile.close()


def get_free_dics(indir, og_list):

    og_dnds_dic = {}
    og_ds_dic = {}
    og_dn_dic = {}
    species_list = []
    for og_file in glob("%s/og_*.alt" % (indir)):        
        cur_og = int(og_file.split("og_")[1].split(".alt")[0])
        if cur_og not in og_list:
            continue
        reader = open(og_file, 'rU')
        dnds_tree = False
        ds_tree = False
        dn_tree = False
        ds_dic = {}
        dn_dic = {}
        dnds_dic = {}
        first_line = True
        for line in reader:
            if first_line:
                seq_len = int(line.strip().split()[1])
                if seq_len < 900:
                    break
                first_line = False

            # if dnds_tree:
            #     dnds = PhyloTree(line.strip().replace("#", ":"))
            #     print dnds
            #     for leaf in dnds:
            #         if leaf.dist < 4 and leaf.dist > 0.0001:
            #             if ds_dic[leaf.name] > 0.001 and ds_dic[leaf.name] < 10:
            #                 if dn_dic[leaf.name] > 0.001 and dn_dic[leaf.name] < 10:
            #                     dnds_dic[leaf.name] = leaf.dist
            #         if leaf.name not in species_list:
            #             species_list.append(leaf.name)
                dnds_tree = False
            if ds_tree:
                ds = PhyloTree(line.strip())
                og_ds_dic[cur_og] = ds
                ds_tree = False
            if dn_tree:
                dn = PhyloTree(line.strip())
                og_dn_dic[cur_og] = dn
                dn_tree = False
            if line.strip() == "dS tree:":
                ds_tree = True
                continue
            if line.strip() == "dN tree:":
                dn_tree = True
                continue
            if line.strip() == "w ratios as labels for TreeView:":
                dnds_tree = True
                continue
        if seq_len < 900:
            continue
        og_dnds_dic[cur_og] = dnds_dic
    return og_ds_dic, og_dn_dic


#def ds_time_correlations(og_list, indir, outdir, calib_tree, trait_tree, og_ds_dic):
def ds_time_correlations(og_list, indir, outdir, calib_tree, trait_tree, og_ds_dic, categorical):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    time_tree = PhyloTree(calib_tree)
    time_tree.unroot()
    time_blengths = get_blengths(time_tree, time_tree)
    pval_dic = {}
    rho_dic = {}
    rando_rho_dic = {}
    outfile = open("%s/ds_time_correlations.txt" % outdir, 'w')
    random_file = open("%s/ds_time_correlations_randomtrait.txt" % outdir, 'w')
    for og, og_ds_tree in og_ds_dic.items():
        og_file = open("%s/og_%s_traits.txt" % (outdir, og), 'w')
        og_blens = get_blengths(og_ds_tree, time_tree)
        standard_dic = {}
        trait_dic = get_blengths(PhyloTree(trait_tree), time_tree)
        traits_list = []
        blens_list = []
        for k, v in og_blens.items():
            standard_dic[k] = v / time_blengths[k]
            traits_list.append(trait_dic[k])
            blens_list.append(standard_dic[k])
            og_file.write("%s\t%s\t%s\n" % (k, trait_dic[k], standard_dic[k]))
        og_file.close()
        if not categorical:
            rho, pval = scipy.stats.spearmanr(traits_list, blens_list)
        else:
            rho, pval = scipy.stats.kendalltau(blens_list, traits_list)
        outfile.write("%s\t%s\t%s\n" % (og, rho, pval))
        pval_dic[og] = pval
        rho_dic[og] = rho
    outfile.close()
    trait_dic = get_blengths(PhyloTree(trait_tree), time_tree)

    for iternum in range(0,1):
        print iternum
        rando_rhos = []
        traits_list = trait_dic.values()
        random.shuffle(traits_list)
        for og, og_ds_tree in og_ds_dic.items():
            og_blens = get_blengths(og_ds_tree, time_tree)
            standard_dic = {}
            blens_list = []
            for k, v in og_blens.items():
                standard_dic[k] = v / time_blengths[k]
                blens_list.append(standard_dic[k])

            if not categorical:
                rando_rho, rando_pval = scipy.stats.spearmanr(traits_list, blens_list)
            else:
                rando_rho, rando_pval = scipy.stats.kendalltau(blens_list, traits_list)
            rando_rhos.append(rando_rho)
        iter_med = numpy.median(rando_rhos)
        random_file.write("%s\t%s\n" % (iternum, iter_med))
        random_file.flush()
        rando_rho_dic[iternum] = iter_med
    random_file.close()
    print numpy.median(rho_dic.values())
    print numpy.median(rando_rho_dic.values())
    print numpy.max(rando_rho_dic.values())
    print numpy.min(rando_rho_dic.values())


def bootstrapping_ds_time_correlations(og_list, indir, outdir, calib_tree, trait_tree, bootstrap_taxa, categorical):       

    og_ds_dic, og_dn_dic = get_free_dics(indir, og_list)
    ds_time_correlations(og_list, indir, outdir, calib_tree, trait_tree, og_ds_dic, categorical)
    if bootstrap_taxa == False:
        return
    leaf_list = []
    for leaf in PhyloTree(calib_tree).get_leaves():
        leaf_list.append(leaf.name)

    for taxon in leaf_list:
        print taxon
        cur_leaf_list = []
        for leaf in leaf_list:
            if leaf != taxon:
                cur_leaf_list.append(leaf)
        boot_time_tree = PhyloTree(calib_tree)
        boot_time_tree.prune(cur_leaf_list, preserve_branch_length = True)
        boot_time_tree = boot_time_tree.write(format = 5)
        boot_trait_tree = PhyloTree(trait_tree)
        boot_trait_tree.prune(cur_leaf_list)
        boot_trait_tree = boot_trait_tree.write(format = 5)
        taxon_og_ds_dic = {}
        for og, ds_tree in og_ds_dic.items():
            cur_tree = copy.deepcopy(ds_tree)
            cur_tree.prune(cur_leaf_list)
            taxon_og_ds_dic[og] = cur_tree
        ds_time_correlations(og_list, indir, "%s/%s_excluded" % (outdir, taxon), boot_time_tree, boot_trait_tree, taxon_og_ds_dic, categorical)

    

def marine_test(phylo_dic, full_tree, fore_list, exclude_list, outdir):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    full_phylo = PhyloTree(full_tree)
    gw_blengths = {}
    og_blengths = {}
    p_dic = {}
    for cur_og, og_phylo in phylo_dic.items():
        cur_blengths = get_blengths(og_phylo, full_phylo)
        og_blengths[cur_og] = cur_blengths
        for name, dist in cur_blengths.items():
            if name not in gw_blengths.keys():
                gw_blengths[name] = []
            gw_blengths[name].append(dist)
    mean_blengths = {}
    for name, dist_list in gw_blengths.items():
        mean_blengths[name] = numpy.median(dist_list)
    og_list = []
    test_results_p = []
    in_mean_dic = {}
    out_mean_dic = {}

    for cur_og, cur_blengths in og_blengths.items():
        p_val, in_mean, out_mean = wilcox_test(cur_blengths, mean_blengths, fore_list, exclude_list)
        if p_val:
            og_list.append(cur_og)
            test_results_p.append(p_val)
            in_mean_dic[cur_og] = in_mean
            out_mean_dic[cur_og] = out_mean
            p_dic[cur_og] = p_val
    for node, dist_list in gw_blengths.items():
        outfile = open("%s/%s_blens.txt" % (outdir, node), 'w')
        for blen in dist_list:
            outfile.write("%s\n" % (blen))
        outfile.close()
#    pval_corr = smm.multipletests(test_results_p, alpha = 0.05, method = 'fdr_bh')[1]
    pval_corr = test_results_p
    p_dic = {}
    for x in range(len(pval_corr)):
        p_dic[og_list[x]] = pval_corr[x]
    sorted_pvals = sorted(p_dic.items(), key = lambda x: x[1])
    outfile = open("%s/marine_tests.txt" % outdir, 'w')
    outfile.write("OG\tpval\t\n")
    for pval in sorted_pvals:
        outfile.write("%s\t%s\t%s\t%s\t%s\n" % (pval[0], pval[1], p_dic[pval[0]], in_mean_dic[pval[0]], out_mean_dic[pval[0]]))
    outfile.close()

def wilcox_test(blength_dic, standard_dic, fore_list, exclude_list):
    ingroup = []
    outgroup = []
    for name, dist in blength_dic.items():
        if dist < 0.00001:
            continue
#        if standard_dic[name] < 4:
#            continue
        if dist / standard_dic[name] < 0.0001:
            continue
        if name in fore_list:
            ingroup.append(dist / standard_dic[name])
        elif name in exclude_list:
            continue
        else:
            outgroup.append(dist / standard_dic[name])
#       if not dist:
#           print outgroup
#           print blength_dic
#           print standard_dic
    if len(ingroup) <=2 or len(outgroup) <=2:
        return False, -9, -9
#    print ingroup
#    print outgroup
    u_val, p_val = scipy.stats.ranksums(ingroup, outgroup)
    return p_val, numpy.mean(ingroup), numpy.mean(outgroup)

def kendall_test(blength_dic, standard_dic, fore_list, exclude_list):
    ingroup = []
    outgroup = []
    for name, dist in blength_dic.items():
        if dist < 0.00001:
            continue
#        if standard_dic[name] < 4:
#            continue
        if dist / standard_dic[name] < 0.0001:
            continue
        if name in fore_list:
            ingroup.append(dist / standard_dic[name])
        elif name in exclude_list:
            continue
        else:
            outgroup.append(dist / standard_dic[name])
#       if not dist:
#           print outgroup
#           print blength_dic
#           print standard_dic
    if len(ingroup) <=2 or len(outgroup) <=2:
        return False, -9, -9
#    print ingroup
#    print outgroup
#    u_val, p_val = scipy.stats.ranksums(ingroup, outgroup)
    cor_list = ingroup + outgroup
    phen_list = []
    for x in range(len(ingroup)):
        phen_list.append(1)
    for x in range(len(outgroup)):
        phen_list.append(0)
    u_val, p_val = scipy.stats.kendalltau(cor_list, phen_list)
    return p_val, numpy.median(ingroup), numpy.median(outgroup)

        
def get_blengths(og_phylo, full_tree):
    #get the mean branch lengths for every branch encountered
    blength_dic = {}
    node_count = 0
    og_leaves = []
    for leaf in og_phylo.get_leaves():
        og_leaves.append(leaf.name)
    for node in full_tree.traverse("postorder"):
        if node.is_leaf():
            sibling_good = False
            if node.name in og_leaves:
                par = node.up
                for sibling in par.traverse("postorder"):
                    if compare_nodes(sibling, node):
                        continue
                    if sibling.name in og_leaves:
                        sibling_good = True
                        break
                grand_par = par.up
                if not grand_par:
                    grand_par_good = True
                else:
                    for child in grand_par.traverse("postorder"):
                        if compare_nodes(node, child):
                            continue
                        else:
                            for leaf in child.get_leaves():
                                if leaf.name in og_leaves:
                                    grand_par_good = True
                                    break
            if sibling_good and grand_par_good:
                if node.name not in blength_dic.keys():
                    blength_dic[node.name] = og_phylo.get_leaves_by_name(node.name)[0].dist
        else:
            node_count += 1
            node.add_features(label = node_count)
            good_child_count = 0
            good_parent = False
            good_grand_parent = False
            for child in node.children:
                for leaf in child.get_leaves():
                    if leaf.name in og_leaves:
                        good_child_count += 1
                        break
            par = node.up
            if not par:
                good_parent = True
                grand_par = False
            else:
                for sibling in par.children:
                    if compare_nodes(sibling, node):
                        continue
                    for leaf in sibling.get_leaves():
                        if leaf.name in og_leaves:
                            good_parent = True
                grand_par = par.up
            if not grand_par:
                good_grand_parent = True
            else:
                for sibling in grand_par.children:
                    if compare_nodes(sibling, par):
                        continue
                    for leaf in sibling.get_leaves():
                        if leaf.name in og_leaves:
                            good_grand_parent = True
            if good_parent and good_grand_parent and good_child_count == 2:
                leaf_list = []
                for leaf in node.get_leaves():
                    for og_leaf in og_phylo.get_leaves():
                        if leaf.name == og_leaf.name:
                            leaf_list.append(leaf.name)
                if node_count not in blength_dic.keys():
                    blength_dic[node_count] = og_phylo.get_common_ancestor(leaf_list).dist
    return blength_dic
"""
def get_blengths(phylo_dic, full_tree):
    #get the mean branch lengths for every branch encountered
    blength_dic = {}
    node_count = 0
    for cur_og, og_phylo in phylo_dic.items():
        node_count = 0
        og_leaves = []
        for leaf in og_phylo.get_leaves():
            og_leaves.append(leaf.name)
        for node in full_tree.traverse("postorder"):
            if node.is_leaf():
                sibling_good = False
                if node.name in og_leaves:
                    par = node.up
                    for sibling in par.traverse("postorder"):
                        if compare_nodes(sibling, node):
                            continue
                        if sibling.name in og_leaves:
                            sibling_good = True
                            break
                    grand_par = par.up
                    if not grand_par:
                        grand_par_good = True
                    else:
                        for child in grand_par.traverse("postorder"):
                            if compare_nodes(node, child):
                                continue
                            else:
                                for leaf in child.get_leaves():
                                    if leaf.name in og_leaves:
                                        grand_par_good = True
                                        break
                if sibling_good and grand_par_good:
                    if node.name not in blength_dic.keys():
#                        blength_dic[node.name] = []
#                    blength_dic[node.name].append(og_phylo.get_leaves_by_name(node.name)[0].dist)
                        blength_dic[node.name] = og_phylo.get_leaves_by_name(node.name)[0].dist
            else:
                node_count += 1
                node.add_features(label = node_count)
                good_child_count = 0
                good_parent = False
                good_grand_parent = False
                for child in node.children:
                    for leaf in child.get_leaves():
                        if leaf.name in og_leaves:
                            good_child_count += 1
                            break
                par = node.up
                if not par:
                    good_parent = True
                    grand_par = False
                else:
                    for sibling in par.children:
                        if compare_nodes(sibling, node):
                            continue
                        for leaf in sibling.get_leaves():
                            if leaf.name in og_leaves:
                                good_parent = True
                    grand_par = par.up
                if not grand_par:
                    good_grand_parent = True
                else:
                    for sibling in grand_par.children:
                        if compare_nodes(sibling, par):
                            continue
                        for leaf in sibling.get_leaves():
                            if leaf.name in og_leaves:
                                good_grand_parent = True
                if good_parent and good_grand_parent and good_child_count == 2:
                    leaf_list = []
                    for leaf in node.get_leaves():
                        for og_leaf in og_phylo.get_leaves():
                            if leaf.name == og_leaf.name:
                                leaf_list.append(leaf.name)
                    if node_count not in blength_dic.keys():
#                        blength_dic[node_count] = []
#                    blength_dic[node_count].append(og_phylo.get_common_ancestor(leaf_list).dist)
                        blength_dic[node_count] = og_phylo.get_common_ancestor(leaf_list).dist
#    for k, v in blength_dic.items():
#        print k
#        print numpy.mean(v)
#    print blength_dic
    return blength_dic
"""                        
def compare_nodes(node1, node2):
    node1_leaves = node1.get_leaves()
    leaf1_names = []
    node2_leaves = node2.get_leaves()
    leaf2_names = []
    for leaf1 in node1_leaves:
        leaf1_names.append(leaf1.name)
    for leaf2 in node2_leaves:
        leaf2_names.append(leaf2.name)
    for leaf1 in leaf1_names:
        if leaf1 not in leaf2_names:
            return False
    return True

def results_exist(outdir, og_num):
    if os.path.exists("%s/og_%s.alt" % (outdir, og_num)):
        reader = open("%s/og_%s.alt" % (outdir, og_num), 'rU')
        if "Time used:" in reader.readlines()[-1]:
            return True
    if os.path.exists("%s/og_%s.anc" % (outdir, og_num)):
        reader = open("%s/og_%s.anc" % (outdir, og_num), 'rU')
        if "Time used:" in reader.readlines()[-1]:
            return True
    return False

def paml_test(og_list, foreground, test_type, indir, outdir, phylogeny_file, num_threads, use_gblocks, min_taxa, exclude_taxa):
    #performs paml test on all OG's in list
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    og_file_list = []
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    for cur_og in og_list:
        if results_exist(outdir, cur_og):
            continue
        if test_type == "model_d":
            prep_paml_files(cur_og, indir, outdir, "model_d", phylogeny_file, test_type, use_gblocks)
        elif test_type == "free":
            no_stops = prep_paml_files(cur_og, indir, outdir, "free", phylogeny_file, test_type, min_taxa, use_gblocks, exclude_taxa)
            if not no_stops:
                print "OG %s appears to have unexpected stop codons" % cur_og
                continue
        elif test_type == "ancestral":
            if foreground == "ncar":
                cur_out_dir = "%s/ncar_%s" % (outdir, cur_og)
            else:
                cur_out_dir = "%s/OG_%s" % (outdir, cur_og)
            is_codons = prep_paml_files(cur_og, indir, outdir, foreground, phylogeny_file, test_type, min_taxa, use_gblocks, exclude_taxa)
        elif test_type == "aaml_blengths":
            is_codons = prep_paml_files(cur_og, indir, outdir, "aaml_blengths", phylogeny_file, test_type, min_taxa, use_gblocks, exclude_taxa)
        elif test_type == "RELAX":
            if os.path.exists("%s/og_%s_relax_unlabeledback.txt" % (outdir, cur_og)):
                finished = False
                reader = open("%s/og_%s_relax_unlabeledback.txt" % (outdir, cur_og), 'rU')
                for line in reader:
                    if line.startswith("Likelihood ratio test"):
                        finished = True
                if finished:
                    continue
            is_codons = prep_paml_files(cur_og, indir, outdir, foreground, phylogeny_file, test_type, min_taxa, use_gblocks, exclude_taxa)
        elif test_type == "aBSREL":
            if os.path.exists("%s/og_%s_absrel.txt" % (outdir, cur_og)):
                finished = False
                reader = open("%s/og_%s_absrel.txt" % (outdir, cur_og), 'rU')
                for line in reader:
                    if line.startswith("Likelihood ratio test"):
                        finished = True
                if finished:
                    continue
            is_codons = prep_paml_files(cur_og, indir, outdir, foreground, phylogeny_file, test_type, min_taxa, use_gblocks, exclude_taxa)

        else:
            taxon_present = prep_paml_files(cur_og, indir, outdir, foreground, phylogeny_file, test_type, use_gblocks)
            if not taxon_present:
                continue
            ifshort = taxon_present
        if is_codons:
#        if ifshort != "too short":
#        paml_tests.ncar_ancestor_reconstruction([cur_og, outdir])

            work_list.append([cur_og, outdir])
    print work_list
#        paml_tests.relax_worker([cur_og, outdir])
#        paml_tests.absrel_worker([cur_og, outdir])
#        paml_tests.aaml_worker([cur_og, outdir])
    if test_type == "bs":
        pool.map_async(paml_tests.branch_site_worker, work_list).get(9999999)
    elif test_type == "branch":
        pool.map_async(paml_tests.branch_worker, work_list).get(9999999)
    elif test_type == "branchpos":
        pool.map_async(paml_tests.branch_positive_worker, work_list).get(9999999)
    elif test_type == "free":
        pool.map_async(paml_tests.free_ratios_worker, work_list).get(9999999)
    elif test_type == "ancestral":
        if foreground == "ncar":
            pool.map_async(paml_tests.ncar_ancestor_reconstruction, work_list).get(9999999)            
        else:
            pool.map_async(paml_tests.ancestor_reconstruction, work_list).get(9999999)
    elif test_type == "aaml_blengths":
        pool.map_async(paml_tests.aaml_worker, work_list).get(9999999)
    elif test_type == "RELAX":
        pool.map_async(paml_tests.relax_worker, work_list).get(9999999)
    elif test_type == "aBSREL":
        pool.map_async(paml_tests.absrel_worker, work_list).get(9999999)
    

def baseml_blengths(ncar_list, aligndir, outdir, treefile, num_threads, remove_list):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    for ncar in ncar_list:
        work_list.append([ncar, aligndir, outdir, treefile, remove_list])
 #       print ncar
#        baseml_worker([ncar, aligndir, outdir, treefile, remove_list])
    pool.map_async(baseml_worker, work_list).get(9999999)

def baseml_worker(param_list):
    ncar = param_list[0]
    aligndir = param_list[1]
    outdir = param_list[2]
    treefile = param_list[3]
    remove_list = param_list[4]
    file_info = ncar.split("_")
    if not os.path.exists("%s/ncar_%s.afa.trimal" % (aligndir, ncar)):
        if os.path.exists("%s/ncar_%s.afa" % (aligndir, ncar)):
            print "trimal failed on %s. Continuing with next locus." % ncar
            return
    reader = SeqIO.parse("%s/ncar_%s.afa.trimal" % (aligndir, ncar), format = 'fasta')
    seq_dic = {}
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    outfile = open("%s/%s.afa" % (outdir, ncar), 'w')
    taxa_list = []
    for seq_name, seq in seq_dic.items():
        cur_taxa = seq_name.split(".")[0]
        if cur_taxa not in remove_list:
            taxa_list.append(cur_taxa)
            outfile.write(">%s\n%s\n" % (cur_taxa, seq))
    outfile.close()
    tree = PhyloTree(treefile)
    tree.prune(taxa_list)
    tree.unroot()
    tree_str = tree.write(format = 5)
    cons_file = open("%s/%s.constraint" % (outdir, ncar), 'w')
    cons_file.write(tree_str)
    cons_file.close()
    cml = baseml.Baseml(alignment = "%s/%s.afa" % (outdir, ncar), tree = "%s/%s.constraint" % (outdir, ncar), out_file = "%s/%s.alt" % (outdir, ncar), working_dir = "%s/%s_working" % (outdir, ncar))
    cml.set_options(runmode=0,fix_blength=0,model=7, clock = 0, Mgene = 0, fix_kappa = 0, kappa = 2, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6, verbose = True)
    try:
        cml.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/baseml", verbose = True)
    except:
        if os.path.exists("%s/%s.alt" % (outdir, ncar)):
            os.remove("%s/%s.alt" % (outdir, ncar))
        print "%s experienced an error" % ncar

def basic_gff(gff_file):
    reader = open(gff_file, 'rU')
    gff_dic = {}
    for line in reader:
        cur_line = line.split()
        if cur_line[2] in ["mRNA", "ncar", "cnee"]:
            cur_name = cur_line[-1].split("Name=")[1][0:-1]
            gff_dic[cur_name] = (cur_line[0], int(cur_line[3]), int(cur_line[4]), cur_line[6])
    return gff_dic


def neighbor_genes(target, coding_dic):
    cur_chrom = target[0]
    cur_start = target[1]
    cur_end = target[2]
    cur_mid = (cur_end + cur_start) / 2.0
    cur_dir = target[3]
    back_list = []
    for gene, coords in coding_dic.items():
        gene_chr = coords[0]
        gene_start = coords[1]
        gene_end = coords[2]
        gene_dir = coords[3]
        if gene_chr == cur_chrom:
            if gene_dir == "+":
                if abs(gene_start - cur_mid) < 1000000:
                    back_list.append(gene)
            else:
                if abs(gene_end - cur_mid) < 1000000:
                    back_list.append(gene)
    return back_list

def neighbor_fourfs(ancestral_dir, back_list, species_list, outfile):
    SPECIES_LIST = species_list
    full_dic = {}
    for species in SPECIES_LIST:
        full_dic[species] = []
    seq_len = 0
    seq_count = 0
    for og in back_list:
        these_species = []
        og_dic = {}

        alignment = "%s/og_%s_working/4fold.nuc" % (ancestral_dir, og)
        reader = open(alignment, 'rU')
        next(reader)
        next(reader)
        seq_line = False
        for line in reader:
            if "codons included" in line:
                break
            if line == "\n":
                seq_line = False
                continue
            if not seq_line:
                cur_id = line.strip().upper()
                seq_line = True
            else:
                cur_seq = line.strip()
                if cur_id in SPECIES_LIST:
#                    full_dic[cur_id].append(cur_seq)
                    og_dic[cur_id] = cur_seq
#                    these_species.append(cur_id)
                seq_line = False
        if len(og_dic) == len(species_list):
            for species, seq in og_dic.items():
                full_dic[species].append(seq)
            seq_count += 1
            seq_len += len(seq)
            
    writer = open(outfile, 'w')
    for species, seq_list in full_dic.items():
        writer.write(">%s\n%s\n" % (species, "".join(seq_list)))
    writer.close()
    return seq_count

def inspecies_ncar_translate(ncar_ortho_dic, inspecies, inspecies_gff):
    inspecies_dic = {}
    for ncar, loci_dic in ncar_ortho_dic.items():
        cur_ncar = loci_dic.get(inspecies, False)
        if cur_ncar:
            inspecies_dic[ncar] = inspecies_gff[cur_ncar]
    return inspecies_dic

def hyphy_noncoding(ncar_list, aligndir, outdir, treefile, num_threads, remove_list, ancestral_dir, ncar_coords, coding_coords, ortho_dic, outspecies):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    
    for ncar in ncar_list:
        work_list.append([ncar, aligndir, outdir, treefile, remove_list, ncar_coords[ncar], coding_coords, ancestral_dir, ortho_dic, outspecies])
#        hyphy_noncoding_worker([ncar, aligndir, outdir, treefile, remove_list, ncar_coords[ncar], coding_coords, ancestral_dir, ortho_dic, outspecies])
    pool.map_async(hyphy_noncoding_worker, work_list).get(9999999)


def hyphy_noncoding_worker(param_list):
    ncar = param_list[0]
    errorfile = open("/Genomics/kocherlab/berubin/comparative/halictids/halictid_cactus_ncars/errfiles/errors_%s.txt" % ncar, 'w')
    errorfile.write("%s\n" % ncar)
 
    aligndir = param_list[1]
    outdir = param_list[2]
    treefile = param_list[3]
    remove_list = param_list[4]
    target_coord = param_list[5]
    coding_coords = param_list[6]
    ancestral_dir = param_list[7]
    ortho_dic = param_list[8]
    outspecies = param_list[9]
    file_info = ncar.split("_")
    if not os.path.exists("%s/ncar_%s.afa.trimal" % (aligndir, ncar)):
        if os.path.exists("%s/ncar_%s.afa" % (aligndir, ncar)):
            print "trimal failed on %s. Continuing with next locus." % ncar
            return
    print ncar
#    print "getting neighbors"
    neighbors = neighbor_genes(target_coord, coding_coords)
#    print "neighbors gotten"
    if len(neighbors) < 4:
        print "too few neighbors: %s" % ncar
        return
    neighbor_ogs = []
    for neighbor in neighbors:
        og_num = get_og_num(neighbor, ortho_dic)
        if not og_num:
            continue
        else:
            alignment = "%s/og_%s_working/4fold.nuc" % (ancestral_dir, og_num)
            if os.path.exists(alignment):
                neighbor_ogs.append(og_num)

    if len(neighbor_ogs) < 4:
        print "too few neighbors: %s" % ncar
        return
    reader = SeqIO.parse("%s/ncar_%s.afa.trimal" % (aligndir, ncar), format = 'fasta')
    seq_dic = {}
    too_gappy = []
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    if not os.path.exists("%s/hyphy_%s/" % (outdir, ncar)):
        os.mkdir("%s/hyphy_%s/" % (outdir, ncar))
    outfile = open("%s/hyphy_%s/%s.afa" % (outdir, ncar, ncar), 'w')
    taxa_list = []

    for seq_name, seq in seq_dic.items():
        cur_taxa = seq_name.split(".")[0]
        if cur_taxa not in remove_list:
            taxa_list.append(cur_taxa)
            outfile.write(">%s\n%s\n" % (cur_taxa, seq.upper()))
    outfile.close()

    num_neighbors = neighbor_fourfs(ancestral_dir, neighbor_ogs, taxa_list, "%s/hyphy_%s/%s_fourf.afa" % (outdir, ncar, ncar))
    if num_neighbors < 4:
        print "too few neighbors: %s" % ncar
        return

    tree = PhyloTree(treefile)
    tree.prune(taxa_list)
    tree.unroot()
    tree_str = tree.write(format = 5)
    cons_file = open("%s/hyphy_%s/%s.tree" % (outdir, ncar, ncar), 'w')
    cons_file.write(tree_str)
    cons_file.close()
    tree = PhyloTree("%s/hyphy_%s/%s.tree" % (outdir, ncar, ncar))
    p_list = []
    null_zeta_list = []
    alt_zeta_back_list = []
    alt_zeta_fore_list = []
    pos_select_list = []
    branch_list = []
    consolidated_file = False
    for leaf in tree.get_leaves():
        skip_run = False
        consolidated_file = False
        if leaf.name == outspecies:
            continue
        if os.path.exists("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf.name)):
            reader = open("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf.name), 'rU')
            filelines = reader.readlines()
            if len(filelines) > 0:
                last_line = filelines[-1].strip()
                if last_line.startswith("|--------|"):
#                    print "already done: %s %s" % (ncar, leaf.name)
                    skip_run = True
        if not skip_run and os.path.exists("%s/hyphy_%s/tip_outputs.txt" % (outdir, ncar)):
            reader = open("%s/hyphy_%s/tip_outputs.txt" % (outdir, ncar), 'rU')
            start_counting = False
            line_count = 0
            for line in reader:
                if line.startswith(ncar) and not start_counting:
                    cur_tip = line.strip().split("_")[1]
 #                   print cur_tip
#                    print leaf.name
                    if cur_tip == leaf.name:
#                        print "start counting: %s" % leaf.name
                        start_counting = True
                        continue
                if start_counting:
                    line_count += 1
                if start_counting:
                    if line.startswith(ncar) or line.strip() == "-------------------------------":
                        if line_count > 10:
#                            print "already done: %s %s" % (ncar, leaf.name)
                            skip_run = True
                            consolidated_file = True
                            break
        cmd = ["HYPHYMP", "CPU=1", "/Genomics/kocherlab/berubin/local/src/TestForPositiveSelection/nonCodingSelection.bf", "Intronicsites",  "%s/hyphy_%s/%s_fourf.afa" % (outdir, ncar, ncar), "%s/hyphy_%s/%s.tree" % (outdir, ncar, ncar), "Yes", leaf.name, "%s/hyphy_%s/%s.afa" % (outdir, ncar, ncar), "Null2-Alternate2", "NEB"]
        if not skip_run:
#            print "starting hyphy: %s %s" % (ncar, leaf.name)
            
#            print leaf.name
            with open("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf.name), 'w') as outfile:
                subprocess.call(cmd, stdout = outfile)
            outfile.close()
        if not consolidated_file:
#            print "no consolidation %s" % leaf.name
            locus_file = open("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf.name), 'rU')
            cur_p, null_zeta, alt_zeta_back, alt_zeta_fore, pos_list = read_hyphy_noncoding_branch(locus_file.readlines())
        else:
            locus_switch = False
            locus_lines = []
            consol_file = open("%s/hyphy_%s/tip_outputs.txt" % (outdir, ncar), 'rU')
            for line in consol_file:
                if line.startswith("%s_%s" % (ncar, leaf.name)):
                    locus_switch = True
                if locus_switch:
                    if line.strip() == "-------------------------------":
                        locus_switch = False
                        break
                    locus_lines.append(line)
            cur_p, null_zeta, alt_zeta_back, alt_zeta_fore, pos_list = read_hyphy_noncoding_branch(locus_lines)

#        cur_p, null_zeta, alt_zeta_back, alt_zeta_fore, pos_list = read_hyphy_noncoding_branch("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf.name))
        p_list.append(cur_p)
        null_zeta_list.append(null_zeta)
        alt_zeta_back_list.append(alt_zeta_back)
        alt_zeta_fore_list.append(alt_zeta_fore)
        pos_select_list.append(len(pos_list))
        branch_list.append(leaf.name)
    node_num = 0
    test_node_list = []
    node_name_dic = {}
    for node in tree.traverse("preorder"):
        skip_run = False
        if node_num == 0:
            node_num += 1
            continue
        leaves = node.get_leaves()
        if len(leaves) == 1:
            node_num += 1
            continue
        leaf_list = []
        for leaf in leaves:
            leaf_list.append(leaf.name)
        test_node_str = "Node" + str(node_num)
        leaf_str = ",".join(leaf_list)
        if os.path.exists("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf_str)):
            reader = open("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf_str), 'rU')
            filelines = reader.readlines()
            if len(filelines) > 0:
                last_line = filelines[-1].strip()
                if last_line.startswith("|--------|"):
#                    print "already done: %s %s" % (ncar, leaf_str) 
                    skip_run = True
        if not skip_run and os.path.exists("%s/hyphy_%s/branch_outputs.txt" % (outdir, ncar)):
            reader = open("%s/hyphy_%s/branch_outputs.txt" % (outdir, ncar), 'rU')
            start_counting = False
            line_count = 0
            for line in reader:
                if line.startswith(ncar) and not start_counting:
                    cur_tip = line.strip().split("_")[1]
                    if cur_tip == leaf_str:
#                        print "start counting: %s" % leaf_str
                        start_counting = True
                        continue
                if start_counting:
                    line_count += 1
                if start_counting:
                    if line.startswith(ncar) or line.strip() == "-------------------------------":
                        if line_count > 10:
#                            print "already done: %s %s" % (ncar, leaf_str)
                            skip_run = True
                            consolidated_file = True
                            break


        cmd = ["HYPHYMP", "CPU=1", "/Genomics/kocherlab/berubin/local/src/TestForPositiveSelection/nonCodingSelection.bf", 'Intronicsites',  "%s/hyphy_%s/%s_fourf.afa" % (outdir, ncar, ncar), "%s/hyphy_%s/%s.tree" % (outdir, ncar, ncar), "Yes", test_node_str, "%s/hyphy_%s/%s.afa" % (outdir, ncar, ncar), "Null2-Alternate2", "NEB"]
        if not skip_run:
            with open("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf_str), 'w') as outfile:
                subprocess.call(cmd, stdout = outfile)
            outfile.close()



        if not consolidated_file:
            locus_file = open("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf_str), 'rU')
            cur_p, null_zeta, alt_zeta_back, alt_zeta_fore, pos_list = read_hyphy_noncoding_branch(locus_file.readlines())
        else:
            locus_switch = False
            locus_lines = []
            consol_file = open("%s/hyphy_%s/branch_outputs.txt" % (outdir, ncar), 'rU')
            for line in consol_file:
                if line.startswith("%s_%s" % (ncar, leaf_str)):
                    locus_switch = True
                if locus_switch:
                    if line.strip() == "-------------------------------":
                        locus_switch = False
                        break
                    locus_lines.append(line)
            cur_p, null_zeta, alt_zeta_back, alt_zeta_fore, pos_list = read_hyphy_noncoding_branch(locus_lines)



#            node_num += 1
#            locus_file = open("%s/hyphy_%s/%s_%s_noncodingselection.txt" % (outdir, ncar, ncar, leaf_str), 'rU')
#            cur_p, null_zeta, alt_zeta_back, alt_zeta_fore, pos_list = read_hyphy_noncoding_branch(locus_file.readlines())
        p_list.append(cur_p)
        null_zeta_list.append(null_zeta)
        alt_zeta_back_list.append(alt_zeta_back)
        alt_zeta_fore_list.append(alt_zeta_fore)
        pos_select_list.append(len(pos_list))
        branch_list.append(leaf_str)
    pval_corr = smm.multipletests(p_list, alpha = 0.1, method = 'fdr_bh')[1]
    outmiddle = open("%s/%s_noncodingselection_summary.txt" % (outdir, ncar), 'w')
    for x in range(len(branch_list)):
        outmiddle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ncar, p_list[x], pval_corr[x], null_zeta_list[x], alt_zeta_back_list[x], alt_zeta_fore_list[x], branch_list[x], pos_select_list[x]))
    outmiddle.close()
    errorfile.write("done")
    errorfile.close()
    return

def read_hyphy_noncoding_branch(locus_lines):
#    reader = open(filename, 'rU')

    cur_p = -1
    null_zeta = -1
    alt_zeta_back = -1
    alt_zeta_fore = -1
    null_model = False
    alt_model = False
    neb_results = False
    neg_dic = {}
    neut_dic = {}
    pos_dic = {}
    pos_list = []
    for line in locus_lines:
        if "*** Null model ***" in line:
            null_model = True
            continue
        if "*** Alternate model ***" in line:
            alt_model = True
            continue
        if null_model:
            if line.startswith("zeta0:"):
                null_zeta = float(line.split(": ")[1])
                null_model = False
                continue
        if alt_model:
            if line.startswith("zeta0:"):
                alt_zeta_back = float(line.split(": ")[1])
                continue
            elif line.startswith("zeta2:"):
                alt_zeta_fore = float(line.split(": ")[1])
                alt_model = False
                continue
        if line.startswith("LRT p-value"):
            cur_p = float(line.split(" = ")[1])
        if line.startswith("calculate NEB"):
            neb_results = True
            continue
        if neb_results:
            cur_line = line.split()
            if len(cur_line) < 2:
                continue
            if is_int(cur_line[1]):
                cur_pos = cur_line[1]
                neg_prob = float(cur_line[3].strip("|").strip("*"))
                neut_prob = float(cur_line[4].strip("|").strip("*"))
                pos_prob = float(cur_line[5].strip("|").strip("*"))
                neg_dic[cur_pos] = neg_prob
                neut_dic[cur_pos] = neut_prob
                pos_dic[cur_pos] = pos_prob
                if "*" in cur_line[5]:
                    pos_list.append(cur_pos)
#    print "%s %s %s" % (filename, cur_p, cur_zeta)
    return cur_p, null_zeta, alt_zeta_back, alt_zeta_fore, pos_list

def is_int(instring):
    try:
        int(instring)
        return True
    except:
        return False

def read_baseml_phylos(ncar_list, indir, outdir, outfile):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    output = open("%s/%s" % (outdir, outfile), 'w')
    for ncar in ncar_list:
        if not os.path.exists("%s/%s.alt" % (indir, ncar)):
            print "%s/%s.alt is missing" % (indir, ncar)
            continue
        cur_blens = read_baseml_blengths("%s/%s.alt" % (indir, ncar))
        if cur_blens == "NA":
            continue
        output.write("%s\t%s\n" % (ncar, cur_blens))
    output.close()

def read_baseml_blengths(infile):
    og_file =open(infile, 'rU')
    aa_tree = False
    first_line = True
    line = og_file.readline()
    while True:
        if first_line:
            seq_len = int(line.strip().split()[1])
            first_line = False
        if aa_tree:
            aa = PhyloTree(line.strip().replace("#", ":"))
            return aa.write(format = 5)
        if line.strip().startswith("tree length = "):
            if "-nan" in line:
                return "NA"
            og_file.readline()
            og_file.readline()
            og_file.readline()
            line = og_file.readline()
            aa_tree = True
            continue
        line = og_file.readline()
        if not line:
            break



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
    return float(syns), float(nsyns)

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

def count_ingene_fourfold(seq1, seq2, in_gene):
    #counts the number of differences at fourfold degenerate sites
    #in two sequences. Takes into account the fourfold degenerate
    #polymorphisms in both sequences so as not to overcount divergence
    fourfold_list = changes.fourfold_codons()
    empty_chars = ["N", "-", "X"]
    x = 0
    fourfold_count = 0
    in_gene_alt_dic = in_gene.coding_fixed_align(seq1)
    while x < len(seq1):
        if seq1[x:x+3] in fourfold_list:
            if seq1[x:x+2] == seq2[x:x+2]:
                if seq1[x+2] != seq2[x+2]:
                    in_alt_list = [seq1[x+2]]
                    out_alt_list = [seq2[x+2]]
                    if x+2 in in_gene_alt_dic.keys():
                        in_alt_list.append(in_gene_alt_dic[x+2])
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

def trimal_automated(inalignment, outalignment):
    #run trimAl on given file
    cmd = ["/Genomics/kocherlab/berubin/local/src/trimal/source/trimal", "-automated1", "-in", inalignment, "-out", outalignment]
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

def fsa_coding_align(og_list, indir, outdir, num_threads, iscoding): 
    #run fsa alignments
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    og_file_list = []
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    for og in og_list:
        already_done = False
        if os.path.exists("%s/og_cds_%s.afa" % (outdir, og)): #check whether the alignment of this OG has already occurred. If so, it won't redo it. This is handy in case the orthogroup alignment crashes at some point.
            reader = SeqIO.parse("%s/og_cds_%s.afa" % (outdir, og), format = 'fasta')
            counter = 0
            for rec in reader:
                counter += 1
                if counter >= 2:
                    already_done = True
                    break
        if not already_done:
            og_file = "%s/og_cds_%s.fa" % (indir, og)
            work_list.append([og_file, outdir, iscoding])
    pool.map_async(fsa_align_worker, work_list).get(9999999)

def fsa_ncar_align(og_list, indir, outdir, num_threads, iscoding): 
    #run fsa alignments
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    og_file_list = []
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    for og in og_list:
        already_done = False
        if os.path.exists("%s/ncar_%s.afa" % (outdir, og)): #check whether the alignment of this OG has already occurred. If so, it won't redo it. This is handy in case the orthogroup alignment crashes at some point.
            reader = SeqIO.parse("%s/ncar_%s.afa" % (outdir, og), format = 'fasta')
            counter = 0
            for rec in reader:
                counter += 1
                if counter >= 2:
                    already_done = True
                    break
        if not already_done:
            og_file = "%s/ncar_%s.fasta" % (indir, og)
#            fsa_align_worker([og_file, outdir, iscoding])
            work_list.append([og_file, outdir, iscoding])
    pool.map_async(fsa_align_worker, work_list).get(9999999)


def fsa_align_worker(param_list):
    og_file = param_list[0]
    outdir = param_list[1]
    iscoding = param_list[2]
    #the worker method for multiprocessing the fsa alignments
    STOP_CODONS = ["TGA", "TAA", "TAG"]
    cur_og = og_file.split("/")[-1]
#    og_num = cur_og.split("_")[2].split(".fa")[0]
    og_num = cur_og.split("_")[-1].split(".fa")[0]
    reader = SeqIO.parse(og_file, format = 'fasta')
    if iscoding:
        form_og_file = "%s/og_cds_%s.fa" % (outdir, og_num)
    else:
        form_og_file = "%s/ncar_%s.fa" % (outdir, og_num)
    seqs_formatted = open(form_og_file, 'w')
    for rec in reader:
        cur_seq = str(rec.seq)
        if iscoding:
            if cur_seq[-3:] in STOP_CODONS:
                cur_seq = cur_seq[:-3] + "NNN"
        seqs_formatted.write(">%s\n%s\n" % (rec.id, cur_seq))
    seqs_formatted.close()
    
    if iscoding:
        cmd = ["/Genomics/kocherlab/berubin/local/src/fsa-1.15.9/bin/fsa", "--nucprot", "%s" % form_og_file]
    else:
        cmd = ["/Genomics/kocherlab/berubin/local/src/fsa-1.15.9/bin/fsa", "%s" % form_og_file]
    FNULL = open(os.devnull, 'w')
    if iscoding:
        with open("%s/og_cds_%s.afa" % (outdir, og_num), 'w') as outfile:
            subprocess.call(cmd, stdout = outfile, stderr = FNULL)
        gblock("%s/og_cds_%s.afa" % (outdir, og_num))
        trimal_automated("%s/og_cds_%s.afa" % (outdir, og_num), "%s/og_cds_%s.afa.trimal" % (outdir, og_num))
        
    else:
        with open("%s/ncar_%s.afa" % (outdir, og_num), 'w') as outfile:
            subprocess.call(cmd, stdout = outfile, stderr = FNULL)
#        gblock("%s/ncar_%s.afa" % (outdir, og_num))
        trimal_automated("%s/ncar_%s.afa" % (outdir, og_num), "%s/ncar_%s.afa.trimal" % (outdir, og_num))
    print "%s alignment complete" % og_num
    return

def jarvis_filtering(og_list, indir, outdir, len_min, num_threads):
    if not os.path.isdir(indir):
        os.mkdir(indir)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    spotProblematicSeqsModules.jarvis_filter(og_list, indir, outdir, len_min, num_threads)
    
def sequence_gap_filtering_worker(infile, originalfile, noparafile, outfile, cur_og, min_seq_prop_kept, max_seq_prop_gap, min_seq_len):
    original_seqs = {}
    reader = SeqIO.parse(originalfile, format = 'fasta')
    for rec in reader:
        original_seqs[rec.id] = str(rec.seq)
    outhandle = open(outfile, 'w')
    nopara_handle = open(noparafile, 'w')
    reader = SeqIO.parse(infile, format = 'fasta')
    seq_dic = {}
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    taxa_list = [seq[0:4] for seq in seq_dic.keys()]
    new_seq_list = []
    seqfilt_taxa = []
    seqfilt_seqs = []
    for seq_name, seq in seq_dic.items():
        cur_seq = str(seq).replace("-","").replace("N", "")
        if len(cur_seq) >= min_seq_len: #check min_seq_len
            if len(cur_seq) >= len(original_seqs[seq_name])*min_seq_prop_kept: #check that at least some fraction (min_seq_prop_kept) of the original sequence is retained in the filtered alignment
                gap_count = str(seq).count("-") + str(seq).count("N")
                if gap_count <= max_seq_prop_gap*len(seq): #check that the sequences doesn't have more than some maximum fraction of gaps or unknown sequences (max_seq_prop_gap)
                    outhandle.write(">%s\n%s\n" % (seq_name, str(seq)))
                    seqfilt_seqs.append(seq_name)
                    seqfilt_taxa.append(seq_name[0:4])
                    new_seq_list.append(seq_name)
    outhandle.close()
    para_taxa = [taxa for taxa in seqfilt_taxa if seqfilt_taxa.count(taxa) > 1]
    for seq_name in seqfilt_seqs:
        if seq_name[0:4] not in para_taxa:
            nopara_handle.write(">%s\n%s\n" % (seq_name, str(seq_dic[seq_name])))
    nopara_handle.close()
    return new_seq_list
                

def sequence_gap_filtering(aligned_dir, filtered_dir, nopara_dir, original_dir, og_list, min_seq_prop_kept, max_seq_prop_gap, min_seq_len, index_file):
###Filter out individual sequences based on what fraction of the sequence is unknown or gaps. Also removes short sequences.
    if not os.path.isdir(filtered_dir):
        os.mkdir(filtered_dir)
    if not os.path.isdir(nopara_dir):
        os.mkdir(nopara_dir)
    index_handle = open(index_file, 'w')
    index_handle.write("#og\tnum_tax\tnum_seq\n")
    for cur_og in og_list:
        infile = "%s/og_cds_%s.afa" % (aligned_dir, cur_og)
        originalfile = "%s/og_cds_%s.fa" % (original_dir, cur_og)
        outfile = "%s/og_cds_%s.afa" % (filtered_dir, cur_og)
        noparafile = "%s/og_cds_%s.afa" % (nopara_dir, cur_og)
        new_seq_list = sequence_gap_filtering_worker(infile, originalfile, noparafile, outfile, cur_og, min_seq_prop_kept, max_seq_prop_gap, min_seq_len)
        tax_list = []
        for seq in new_seq_list:
            tax_list.append(seq[0:4])
        tax_count = len(list(set(tax_list)))
        index_handle.write("%s\t%s\t%s\t%s\n" % (cur_og, tax_count, len(new_seq_list), ",".join([seq.split("-")[0] for seq in new_seq_list])))
    index_handle.close()
    

def alignment_column_filtering(aligned_dir, filtered_dir, og_list, nogap_min_count, nogap_min_prop, nogap_min_taxa, required_taxon_dic, num_threads):
    if not os.path.isdir(filtered_dir):
        os.mkdir(filtered_dir)
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    for cur_og in og_list:
        infile = "%s/og_cds_%s.afa" % (aligned_dir, cur_og)
        outfile = "%s/og_cds_%s.afa" % (filtered_dir, cur_og)
        work_list.append([infile, outfile, cur_og, nogap_min_count, nogap_min_prop, nogap_min_taxa, required_taxon_dic])
#        alignment_column_filtering_worker([infile, outfile, cur_og, nogap_min_count, nogap_min_prop, nogap_min_taxa, required_taxon_dic])
    pool.map_async(alignment_column_filtering_worker, work_list).get(9999999)

def alignment_column_filtering_worker(param_list):#infile, outfile, cur_og, nogap_min_count, nogap_min_prop, required_taxon_dic):
    #traverses alignment and removes columns with fewer than nogap_min_count
    #sequences that are not gaps and less then nogap_min_prop (fraction) 
    #sequences that are not gaps. Need to add a taxon requirement as well. 
    #I am imagining that required_taxon_dic might have keys of numbers and
    #values of lists of the taxa that need to be present in at least the key
    #quantity.
    infile = param_list[0]
    outfile = param_list[1]
    cur_og = param_list[2]
    print cur_og
    nogap_min_count = param_list[3]
    nogap_min_prop = param_list[4]
    nogap_min_taxa = param_list[5]
    required_raxon_dic = param_list[6]
    seq_dic = {}
    species_list = []
    bad_columns = []
    align = AlignIO.read(infile, format = 'fasta')
    for column in range(len(align[0])):
        gap_count = align[:,column].count("-") + align[:,column].count("N")
        nogap_count = len(align) - gap_count * 1.0
        nogap_species_list = []
        for rec in align:
            if rec[column] != "-" and rec[column] != "N":
                nogap_species_list.append(rec.id)
        if nogap_count < nogap_min_count:
            bad_columns.append(column)
        if nogap_count / len(align) < nogap_min_prop:
            bad_columns.append(column)
        if len(list(set(nogap_species_list))) < nogap_min_taxa:
            bad_columns.append(column)
    new_seq_dic = {}        
    for rec in align:
        new_seq_dic[rec.id] = []
    for column in range(len(align[0])):
        if column in bad_columns:
            continue
        for rec in align:
            new_seq_dic[rec.id].append(rec[column])
    outfile = open(outfile, 'w')
    for k, v in new_seq_dic.items():
        outfile.write(">%s\n%s\n" % (k, "".join(v)))
    outfile.close()
    

def gblock(inalignment):
    #run gblocks on given file
    cmd = ["/Genomics/kocherlab/berubin/local/src/Gblocks_0.91b/Gblocks", inalignment, "-t=c", "-b5=h"]
    subprocess.call(cmd)


def get_bee_cds():
    #get dictionary containing CDS for all genomes
    seq_dic = {}
    for species in ["AAUR", "APUR", "AVIR", "HLIG", "HQUA", "HRUB", "LALB", "LCAL", "LFIG", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "NMEL"]:
        reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release_v2.1/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (species, species), format = 'fasta')
        seq_dic[species] = {}
        for rec in reader:
            seq_dic[species][rec.id] = str(rec.seq)
#    species = "LALB"
#    reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/LALB/LALB_v3/LALB/LALB_OGS_v1.0_longest_isoform.cds.fasta", format = 'fasta')
#    seq_dic["LALB"] = {}
#    for rec in reader:
#        seq_dic[species][rec.id] = str(rec.seq)
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/annotation/reference_genomes/Dufourea_novaeangliae_v1.1.cds.fa", format = 'fasta')
    species = "Dnov"
    seq_dic[species] = {}
    for rec in reader:
        seq_dic[species][rec.id] = str(rec.seq)
#    reader = SeqIO.parse("/Genomics/kocherlab/berubin/data/nmel/Nmel_v1.0.cds.fa", format = 'fasta')
#    species = "Nmel"
#    seq_dic[species] = {}
#    for rec in reader:
#        seq_dic[species][rec.id] = str(rec.seq)
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/data/mgen/Mgen_v1.0.cds.fa", format = 'fasta')
    species = "Mgen"
    seq_dic[species] = {}
    for rec in reader:
        seq_dic[species][rec.id] = str(rec.seq)

    return seq_dic

def get_bumble_cds():
    seq_dic = {}
    for species in ["BIMP", "PMIB"]:
        reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (species, species), format = 'fasta')
        seq_dic[species] = {}
        for rec in reader:
            seq_dic[species][rec.id] = str(rec.seq)
    return seq_dic

def get_acacia_cds():
    seq_dic = {}
    for species in ["PGRA", "PCON", "PPAL", "PDEN", "PELO", "PCUB", "PFLA", "PPSW"]:
        reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (species, species), format = 'fasta')
        seq_dic[species] = {}
        for rec in reader:
            seq_dic[species][rec.id] = str(rec.seq)
    return seq_dic

def get_PGRAants_cds(base_dir):
    seq_dic = get_cds(base_dir)
    for species in ["PGRA"]:
        reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (species, species), format = 'fasta')
        seq_dic[species] = {}
        for rec in reader:
            seq_dic[species][rec.id[0:-3]] = str(rec.seq)
    return seq_dic

def get_cds_files(cds_dic):
    seq_dic = {}
    for k, v in cds_dic.items():
        reader = SeqIO.parse(v, format = 'fasta')
        seq_dic[k] = {}
        for rec in reader:
            seq_dic[k][rec.id] = str(rec.seq)
    return seq_dic

def get_cds(base_dir):
    seq_dic = {}
    for seqfile in glob("%s/*.fna" % (base_dir)):
        reader = SeqIO.parse(seqfile, format = 'fasta')
        cur_species = seqfile.split(".")[0].split("/")[-1]
        seq_dic[cur_species] = {}
        for rec in reader:
            seq_dic[cur_species][rec.id] = str(rec.seq)
    return seq_dic

def get_prot_seqs():
    #get dictionary containing protein sequences for all genomes
    seq_dic = {}
    for species in ["APUR", "HLIG", "LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LFIG", "AAUR", "AVIR", "HRUB", "HQUA"]:
        reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v1.0_longest_isoform.pep.fasta" % (species, species), format = 'fasta')
        seq_dic[species] = {}
        for rec in reader:
            seq_dic[species][rec.id] = str(rec.seq)
    species = "LALB"
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/LALB/LALB_v3/LALB/LALB_OGS_v1.0_longest_isoform.pep.fasta", format = 'fasta')
    seq_dic["LALB"] = {}
    for rec in reader:
        seq_dic[species][rec.id] = str(rec.seq)
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/annotation/reference_genomes/Dufourea_novaeangliae_v1.1.pep.fa", format = 'fasta')
    species = "Dnov"
    seq_dic[species] = {}
    for rec in reader:
        seq_dic[species][rec.id] = str(rec.seq)
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/data/nmel/Nmel_v1.0.pep.fa", format = 'fasta')
    species = "Nmel"
    seq_dic[species] = {}
    for rec in reader:
        seq_dic[species][rec.id] = str(rec.seq)
    return seq_dic


def write_orthos(ortho_file, seq_dic, paras_allowed, outdir, indexfile):
    #read/parse orthology file and write files containing all sequences.
    #also create an index file that lists the number of taxa in each OG.
    #writes OG's with paralogs in them but just doesn't write sequences
    #from the species with the paralogs.
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if not os.path.isdir(outdir + "_prots"):
        os.mkdir(outdir + "_prots")
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
        outfile = open("%s/og_cds_%s.fa" % (outdir, counter), 'w')
        protfile = open("%s_prots/og_cds_%s.faa" % (outdir, counter), 'w')
        genes_list = []
        for seqs in cur_line[3:]:
            if "*" in seqs:
                continue
            if "," in seqs:
                continue
            cur_seqs = seqs.split(",")
            for seq in cur_seqs:
                genes_list.append(seq)
                cur_species = seq[0:4]
                outfile.write(">%s\n%s\n" % (seq, seq_dic[cur_species][seq].upper()))
                protfile.write(">%s\n%s\n" % (seq, str(Seq.Seq(seq_dic[cur_species][seq].upper()).translate())))
        ref_file.write("%s\t%s\t%s\t%s\n" % (counter, cur_line[0], num_paras, "\t".join(genes_list)))
        outfile.close()
        counter += 1
    ref_file.close()

def write_orthogroups(ortho_dic, seq_dic, outdir, indexfile, min_taxa):
    #read/parse orthology file and write files containing all sequences.
    #also create an index file that lists the number of taxa in each OG.
    #writes OG's with paralogs in them but just doesn't write sequences
    #from the species with the paralogs.
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if not os.path.isdir(outdir + "_prots"):
        os.mkdir(outdir + "_prots")
    ref_file = open(indexfile, 'w')
    ref_file.write("#og\tnum_taxa\tnum_seqs\n")
    counter = 0
    noparas_ortho_dic = {}
    for og, gene_dic in ortho_dic.items():
        noparas_ortho_dic[og] = {}
        for cur_species, genes in gene_dic.items():
            if "," not in genes[0]:
                noparas_ortho_dic[og][cur_species] = genes[0]
    for og, gene_dic in noparas_ortho_dic.items():
        if len(gene_dic.keys()) < min_taxa:
            continue
        ref_file.write("%s\t%s\t%s\n" % (og, len(gene_dic.keys()), len(gene_dic.keys())))
        outfile = open("%s/og_cds_%s.fa" % (outdir, og), 'w')
        protfile = open("%s_prots/og_cds_%s.faa" % (outdir, og), 'w')
        
        for cur_species, genes in gene_dic.items():
            outfile.write(">%s\n%s\n" % (genes, seq_dic[cur_species][genes].upper()))
            protfile.write(">%s\n%s\n" % (genes, str(Seq.Seq(seq_dic[cur_species][genes].upper()).translate())))
        outfile.close()
        protfile.close()

def write_orthoparagroups(ortho_dic, seq_dic, outdir, indexfile, min_taxa):
    #read/parse orthology file and write files containing all sequences.
    #also create an index file that lists the number of taxa in each OG.
    #Does write orthogroups with paralogs but does some filtering:
    #1. Removes sequences less than half the median sequence length.
    #2. Removes species with more than three sequences in orthogroup
    #3. Does not print orthogroups with more than 1.5x the number of sequences as the number of species in the orthogroup.
    #These filters are done in that order!
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if not os.path.isdir(outdir + "_prots"):
        os.mkdir(outdir + "_prots")
        
    ref_file = open(indexfile, 'w')
    ref_file.write("#og\tnum_taxa\tnum_seqs\n")
    counter = 0
    paras_ortho_dic = {}
    for og, gene_dic in ortho_dic.items():
        paras_ortho_dic[og] = {}
        for cur_species, genes in gene_dic.items():
#            if "," in genes[0]:
            paras_ortho_dic[og][cur_species] = genes[0].split(",")
    for og, gene_dic in paras_ortho_dic.items():
        gene_dic = remove_shorter_seqs(gene_dic, seq_dic)
        if len(gene_dic.keys()) < min_taxa:
            continue
        outfile = open("%s/og_cds_%s.fa" % (outdir, og), 'w')
        protfile = open("%s_prots/og_cds_%s.faa" % (outdir, og), 'w')
        seq_count = 0        
        for cur_species, genes in gene_dic.items():
            for gene in genes:
                seq_count += 1
                outfile.write(">%s\n%s\n" % (gene, seq_dic[cur_species][gene].upper()))
                protfile.write(">%s\n%s\n" % (gene, str(Seq.Seq(seq_dic[cur_species][gene].upper()).translate())))
        ref_file.write("%s\t%s\t%s\n" % (og, len(gene_dic.keys()), seq_count))
        outfile.close()
        protfile.close()
    ref_file.close()

def remove_shorter_seqs(gene_dic, seq_dic):
    seq_lens = []
    new_gene_dic = []
    for species, gene_list in gene_dic.items():
        for gene in gene_list:
            seq_lens.append(len(seq_dic[species][gene]))
    av_len = numpy.median(seq_lens)
    gene_count = 0
    for species, gene_list in gene_dic.items(): #remove any gene less than half the median length
        new_genes = []
        for gene in gene_list:
            if len(seq_dic[species][gene]) > av_len / 2.0:
                new_genes.append(gene)
                gene_count += 1
        if len(new_genes) > 0 and len(new_genes) <= 3: #remove any species with more than 3 sequences
            gene_dic[species] = new_genes
        else:
            gene_dic.pop(species)
    if gene_count < 1.5*len(gene_dic.keys()): #only keep orthogroups with less than 1.5*(the number of species). otherwise return an empty dic.
        return gene_dic
    else:
        return {}

def concatenate_for_raxml(input_dir, outfile, og_list, species_list):
    #take list of alignments and concatenate them all into a
    #RAxML formatted fasta file of amino acid sequences. 
    #Uses trimAl filtered alignments.
    SPECIES_LIST = species_list
    full_dic = {}
    for species in SPECIES_LIST:
        full_dic[species] = []
    seq_len = 0
    seq_count = 0
    for og in og_list:
        alignment = "%s/og_cds_%s.afa.trimal" % (input_dir, og)
        reader = SeqIO.parse(alignment, format = 'fasta')
        for rec in reader:
            cur_species = rec.id[0:4]
            cur_seq = Seq.Seq(str(rec.seq).replace("-", "N")).translate()
            full_dic[cur_species].append(str(cur_seq))
        seq_count += 1
        seq_len += len(cur_seq)
    print "total genes used: " + str(seq_count)
    writer = open(outfile, 'w')
    writer.write("%s %s\n" % (len(SPECIES_LIST), seq_len))
    for species, seq_list in full_dic.items():
        writer.write("%s\n%s\n" % (species, "".join(seq_list)))
    writer.close()

def concatenate_fourf_for_raxml(ancestral_dir, outfile, og_list, species_list):
    #take list of alignments and concatenate them all into a
    #RAxML formatted fasta file of amino acid sequences. 
    #Uses trimAl filtered alignments.
    SPECIES_LIST = species_list
    full_dic = {}
    for species in SPECIES_LIST:
        full_dic[species] = []
    seq_len = 0
    seq_count = 0
    for og in og_list:
        alignment = "%s/og_%s_working/4fold.nuc" % (ancestral_dir, og)
        reader = open(alignment, 'rU')
        next(reader)
        next(reader)
        seq_line = False
        for line in reader:

            if "codons included" in line:
                break
            if line == "\n":
                seq_line = False
                continue
            if not seq_line:
                cur_id = line.strip()
                seq_line = True
            else:
                cur_seq = line.strip()
                full_dic[cur_id].append(cur_seq)
                seq_line = False
        seq_count += 1
        seq_len += len(cur_seq)
    print "total genes used: " + str(seq_count)
    writer = open(outfile, 'w')
    writer.write("%s %s\n" % (len(SPECIES_LIST), seq_len))
    for species, seq_list in full_dic.items():
        writer.write("%s\n%s\n" % (species, "".join(seq_list)))
    writer.close()



def gene_trees(og_list, aligndir, outdir, constrained, constraint_tree, num_threads, prots_or_nucs):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    work_list = []
    for og in og_list:
        work_list.append([og, aligndir, outdir, constrained, constraint_tree])
    pool = multiprocessing.Pool(processes = num_threads)
    if prots_or_nucs == "nucs":
        pool.map_async(gene_tree_worker, work_list).get(9999999)
    elif prots_or_nucs == "prots":
        pool.map_async(protein_tree_worker, work_list).get(9999999)

def gene_tree_worker(param_list):
    og = param_list[0]
    aligndir = param_list[1]
    outdir = param_list[2]
    constrained = param_list[3]
    constraint_tree = param_list[4]
    seq_dic = {}
    reader = SeqIO.parse("%s/og_cds_%s.afa" % (aligndir, og), format = 'fasta')
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
        seqlen = len(str(rec.seq))
    outfile = open("%s/og_cds_%s.afa" % (outdir, og), 'w')
    outfile.write("%s %s\n" % (len(seq_dic.keys()), seqlen))
    for k, v in seq_dic.items():
        outfile.write("%s\n%s\n" % (k[0:4], v))
    outfile.close()
    switch_dir = os.getcwd()
    os.chdir(outdir)
    if not constrained:
        cmd = ["raxmlHPC-PTHREADS-SSE3", '-f', 'a', '-x', '12345', '-p', '12345', '-m', 'GTRGAMMA', '-#', '10', '-T', '2', '-s', "og_cds_%s.afa" % (og), '-n', "og_%s.tree" % (og)]
    else:
        cmd = ["raxmlHPC-PTHREADS-SSE3", '-f', 'a', '-x', '12345', '-p', '12345', '-m', 'GTRGAMMA', '-#', '10', '-T', '2', '-s', "og_cds_%s.afa" % (og), '-n', "og_%s.tree" % (og), '-g', constraint_tree, '-o', 'LALB,DNOV']
    subprocess.call(cmd)
    os.chdir(switch_dir)

def protein_tree_worker(param_list):
    og = param_list[0]
    aligndir = param_list[1]
    outdir = param_list[2]
    constrained = param_list[3]
    constraint_tree = param_list[4]
    seq_dic = {}
    reader = SeqIO.parse("%s/og_cds_%s.afa-gb" % (aligndir, og), format = 'fasta')
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq).replace("-", "N")
        seqlen = len(str(rec.seq))
    outfile = open("%s/og_cds_%s.afa" % (outdir, og), 'w')
    outfile.write("%s %s\n" % (len(seq_dic.keys()), seqlen / 3))
    for k, v in seq_dic.items():
        outfile.write("%s\n%s\n" % (k[0:4], str(Seq.Seq(v).translate())))
    outfile.close()
    switch_dir = os.getcwd()
    os.chdir(outdir)
    if not constrained:
        cmd = ["raxmlHPC-PTHREADS-SSE3", '-f', 'a', '-x', '12345', '-p', '12345', '-m', 'PROTGAMMAWAG', '-#', '10', '-T', '2', '-s', "og_cds_%s.afa" % (og), '-n', "og_%s.tree" % (og)]
    else:
#        cmd = ["raxmlHPC-PTHREADS-SSE3", '-f', 'a', '-x', '12345', '-p', '12345', '-m', 'GTRGAMMA', '-#', '10', '-T', '2', '-s', "og_cds_%s.afa" % (og), '-n', "og_%s.tree" % (og), '-g', constraint_tree, '-o', 'AMEL']
        cmd = ["raxmlHPC-PTHREADS-SSE3", '-f', 'a', '-x', '12345', '-p', '12345', '-m', 'PROTGAMMAWAG', '-#', '10', '-T', '2', '-s', "og_cds_%s.afa" % (og), '-n', "og_%s.tree" % (og), '-g', constraint_tree, '-o', 'LALB,DNOV']
    subprocess.call(cmd)
    os.chdir(switch_dir)


def discordance(og_list, aligndir, genetreedir, outdir, species_tree, num_threads):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    work_list = []
    for og in og_list:
        work_list.append([og, aligndir, genetreedir, outdir, species_tree])
    pool = multiprocessing.Pool(processes = num_threads)
    pool.map_async(discordance_worker, work_list).get(9999999)

def discordance_worker(param_list):
    #paralogous sequences are removed be default. I don't know
    #what purpose this would serve if we weren't removing them.
    og = param_list[0]
    print og
    aligndir = param_list[1]
    genetreedir = param_list[2]
    outdir = param_list[3]
    species_tree = param_list[4]
    taxa_list = []
    reader = SeqIO.parse("%s/og_cds_%s.afa" % (aligndir, og), format = 'fasta')
    outalign = open("%s/og_cds_%s.afa" % (outdir, og), 'w')
    seq_dic = {}
    for rec in reader:
        outalign.write(">%s\n%s\n" % (rec.id[0:4], str(rec.seq)))
        taxa_list.append(rec.id[0:4])
    outalign.close()
    tree = PhyloTree(species_tree)
    tree.prune(taxa_list)
    tree.unroot()
    outfile = open("%s/og_%s_consensus.tree" % (outdir, og), 'w')
    outfile.write(tree.write(format = 5))
    outfile.close()
    FNULL = open(os.devnull, 'w')
    cmd = ["/Genomics/kocherlab/berubin/local/src/FastTree", "-gamma", "-nt", "-gtr", "-nome", "-mllen", "-intree", "%s/RAxML_bestTree.og_%s.tree" % (genetreedir, og), "-log", "%s/og_%s_novel.log" % (outdir, og), "%s/og_cds_%s.afa" % (outdir, og)]
    with open("%s/og_%s_novel.ftlen" % (outdir, og), 'w') as outfile:
        subprocess.call(cmd, stdout = outfile, stderr = FNULL)
    outfile.close() 

    cmd = ["/Genomics/kocherlab/berubin/local/src/FastTree", "-gamma", "-nt", "-gtr", "-nome", "-mllen", "-intree", "%s/og_%s_consensus.tree" % (outdir, og), "-log", "%s/og_%s_consensus.log" % (outdir, og), "%s/og_cds_%s.afa" % (outdir, og)]
    with open("%s/og_%s_consensus.ftlen" % (outdir, og), 'w') as outfile:
        subprocess.call(cmd, stdout = outfile, stderr = FNULL)
    outfile.close()
    cmd = ["perl", "/Genomics/kocherlab/berubin/local/src/GammaLogToPaup.pl", "%s/og_%s_consensus.log" % (outdir, og), "%s/og_%s_novel.log" % (outdir, og)]
    with open("%s/og_%s_consensus_novel.txt" % (outdir, og), 'w') as outfile:
        subprocess.call(cmd, stdout = outfile, stderr = FNULL)
    outfile.close()
    cmd = ["/Genomics/kocherlab/berubin/local/src/consel/src/makermt", "--paup", "%s/og_%s_consensus_novel.txt" % (outdir, og)]
    subprocess.call(cmd, stderr = FNULL, stdout = FNULL)
    cmd = ["/Genomics/kocherlab/berubin/local/src/consel/src/consel", "%s/og_%s_consensus_novel" % (outdir, og)]
    subprocess.call(cmd, stderr = FNULL, stdout = FNULL)
    cmd = ["/Genomics/kocherlab/berubin/local/src/consel/src/catpv", "%s/og_%s_consensus_novel" % (outdir, og)]
    with open("%s/og_%s_consensus_novel.consel" % (outdir, og), 'w') as outfile:
        subprocess.call(cmd, stdout = outfile, stderr = FNULL)
    outfile.close()

def read_discordance(consel_dir, og_list, outpath):
    new_og_list = []
    outfile = open("%s/consel_consistency.txt" % outpath, 'w')
    outfile.write("OG\tP\tFDR_P\tnumtax\tbp\n")
    p_dic = {}
    p_list = []
    uncorr_dic = {}
    seq_len_dic = {}
    num_taxa_dic = {}
    for cur_og in og_list:
        reader = SeqIO.parse("%s/og_cds_%s.afa" % (consel_dir, cur_og), format = 'fasta')
        num_taxa = 0
        for rec in reader:
            num_taxa += 1
        seq_len_dic[cur_og] = len(rec.seq)
        num_taxa_dic[cur_og] = num_taxa
#                too_short = True
#                break
#        too_short = False
#        if not too_short:

        consist = read_consel("%s/og_%s_consensus_novel.consel" % (consel_dir, cur_og)) 
        p_list.append(consist)
        uncorr_dic[cur_og] = consist
#            outfile.write("%s\t%s\n" % (cur_og, consist))
#            if consist == "consistency":
#                new_og_list.append(cur_og)
    
    pval_corr = smm.multipletests(p_list, alpha = 0.05, method = 'fdr_i')[1]
    for x in range(len(pval_corr)):
        p_dic[og_list[x]] = pval_corr[x]
    sorted_pvals = sorted(p_dic.items(), key = lambda x: x[1])
    for pval in sorted_pvals:
        outfile.write("%s\t%s\t%s\t%s\t%s\n" % (pval[0], uncorr_dic[pval[0]], pval[1], num_taxa_dic[pval[0]], seq_len_dic[pval[0]]))
        outfile.flush()
    outfile.close()
    return

def read_consel(consel_path):
    reader = open(consel_path, 'rU')
    read_aus = False
    for line in reader:
        if "rank" in line:
            read_aus = True
            continue
        if read_aus:
            cur_line = line.split()
            if len(cur_line) == 0:
                continue
            if int(cur_line[2]) == 2:
                alt_p = float(cur_line[4])
            if int(cur_line[2]) == 1:
                cons_p = float(cur_line[4])
    return cons_p
#    if cons_p < 0.01:       
#        return "discrepancy"
#    else:
#        return "consistency"
    

def read_ortho_index(index_file, min_taxa, paras_allowed):
    #get list of all of the OG's with the minumum taxa
    reader = open(index_file, 'rU')
    og_list = []
    for line in reader:
        if line.startswith("#"):
            continue
        cur_line = line.split()
        if not paras_allowed:
            if int(cur_line[2]) > int(cur_line[1]):
                continue
        if int(cur_line[1]) >= min_taxa:
            og_list.append(int(cur_line[0]))
    return og_list

def remove_aligned_paras(og_list, indir, outdir, parafree_index):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    new_index = open(parafree_index, 'w')
    new_index.write("#og\tnum_tax\tnum_seq\tseqs\n")
    for og in og_list:
        seq_dic = {}
        species_list = []
        reader = SeqIO.parse("%s/og_cds_%s.afa" % (indir, og), format = 'fasta')
        for rec in reader:
            cur_species = rec.id[0:4]
            species_list.append(cur_species)
            seq_dic[rec.id] = str(rec.seq)
        outfile = open("%s/og_cds_%s.afa" % (outdir, og), 'w')
        seq_list = []
        for seq_name, seq in seq_dic.items():
            cur_species = seq_name[0:4]
            if species_list.count(cur_species) == 1:
                outfile.write(">%s\n%s\n" % (seq_name, seq))
                seq_list.append(seq_name)
        outfile.close()
        new_index.write("%s\t%s\t%s\t%s\n" % (og, len(seq_list), len(seq_list), ",".join(seq_list)))
    new_index.close()


def limit_list(og_list, lower_bound, higher_bound):
    #limit the number of OG's to be examined
    new_list = []
    for og in og_list:
        if og < lower_bound:
            continue
        if og >= higher_bound:
            continue
        else:
            new_list.append(og)
    return new_list

def limit_list_ncars(og_list, lower_bound, higher_bound):
    #limit the number of OG's to be examined
    new_list = []
    for ncar_id in og_list:
        if type(ncar_id) is str:
            m = re.search("\d", ncar_id)
            ncar_num = int(ncar_id[m.start():])
        if ncar_num < lower_bound:
            continue
        if ncar_num >= higher_bound:
            continue
        else:
            new_list.append(ncar_id)
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
        if float(cur_line[2]) < 80:
            continue
        out_tuple = (int(cur_line[6]), int(cur_line[7]))
        in_tuple = (int(cur_line[8]), int(cur_line[9]))
        coord_pairs[in_tuple] = out_tuple
    if len(coord_pairs.keys()) == 0:
        return (-1, -1), (-1, -1)
    in_tuple_list = coord_pairs.keys()
    for in_coords in coord_pairs.keys():
        in_tuple_list = overlap(in_coords, in_tuple_list)
    out_tuple_list = coord_pairs.values()
 #   for out_coords in coord_pairs.values():
 #       out_tuple_list = overlap(out_coords, out_tuple_list)
    longest_intuple = in_tuple_list[0]
    for in_tuple in in_tuple_list:
        if (in_tuple[1]-in_tuple[0]) > (longest_intuple[1] - longest_intuple[0]):
            longest_intuple = in_tuple
    for intuple in coord_pairs.keys():
        if intuple[0] == longest_intuple[0]:
            starttuple = intuple
        if intuple[1] == longest_intuple[1]:
            endtuple = intuple
    outstart = coord_pairs[starttuple][0]
    outend = coord_pairs[endtuple][1]
#    longest_outtuple = out_tuple_list[0]
#    for out_tuple in out_tuple_list:
#        if (out_tuple[1]-out_tuple[0]) > (longest_outtuple[1] - longest_outtuple[0]):
#            longest_outtuple = out_tuple
    longest_outtuple = (outstart, outend)
    return longest_intuple, longest_outtuple

         
def overlap(mytuple, tuplelist):
    included = False
    i = 0
    merged = []
    while i < len(tuplelist):
        cur_tuple = tuplelist[i]
        if mytuple[0] <= cur_tuple[0] and mytuple[1] >= cur_tuple[0]-50:
#            merged.append(mytuple)
#            merged.append(cur_tuple)
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
   
def lasihali_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file):
    #trims taxa and adds foreground tags for PAML analysis tree
    ingroup_list = ["HLIG", "HRUB", "LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LALB", "LFIG"]
    tree = PhyloTree(phylogeny_file)
    tree.prune(taxa_list)
    tree.unroot()
    leaf_list = []
    for leaf in tree:
        leaf_list.append(leaf.name)
    trimmed_ingroup = []
    for species in ingroup_list:
        if species in leaf_list:
            trimmed_ingroup.append(species)
#    target_branch = tree.get_common_ancestor(trimmed_ingroup).dist            
    tree.get_common_ancestor(trimmed_ingroup).add_features(foreground = "foreground")
#    tree_str = tree.write(format = 5)
    tree_str = tree.write(format = 5, features = ["foreground"])
    tree_list = tree_str.split("[&&NHX:foreground=foreground]")
    beginning = tree_list[0].rsplit(":", 1)
    if "clade" in foreground:
        new_tree = beginning[0] + "$1:" + beginning[1] + tree_list[1]
    else:
        new_tree = beginning[0] + "#1:" + beginning[1] + tree_list[1]
#    tree_str = tree_str.replace("):%s" % target_branch, ")#1:%s" % target_branch)
    outfile = open("%s/og_%s.tree" % (outdir, orthogroup), 'w')
    outfile.write(new_tree)
    outfile.close()

def lasiaugo_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file):
    #trims taxa and adds foreground tags for PAML analysis tree
    ingroup_list = ["LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LALB", "LFIG"]
    augo_list = ["AAUR", "APUR"]
    tree = PhyloTree(phylogeny_file)
    tree.prune(taxa_list)
    tree.unroot()
    leaf_list = []
    for leaf in tree:
        leaf_list.append(leaf.name)
    trimmed_ingroup = []
    for species in ingroup_list:
        if species in leaf_list:
            trimmed_ingroup.append(species)
    trimmed_augos = []
    for species in augo_list:
        if species in leaf_list:
            trimmed_augos.append(species)
#    target_branch = tree.get_common_ancestor(trimmed_ingroup).dist            
    tree.get_common_ancestor(trimmed_ingroup).add_features(foreground = "foreground")
    tree.get_common_ancestor(trimmed_augos).add_features(foreground = "augo_fore")
#    tree_str = tree.write(format = 5)
    tree_str = tree.write(format = 5, features = ["foreground"])

    tree_list = tree_str.split("[&&NHX:foreground=foreground]")
    beginning = tree_list[0].rsplit(":", 1)
    if "clade" in foreground:
        new_tree = beginning[0] + "$1:" + beginning[1] + tree_list[1]
    else:
        new_tree = beginning[0] + "#1:" + beginning[1] + tree_list[1]
    tree_list = new_tree.split("[&&NHX:foreground=augo_fore]")
    beginning = tree_list[0].rsplit(":", 1)
    if "clade" in foreground:
        new_tree = beginning[0] + "$1:" + beginning[1] + tree_list[1]
    else:
        new_tree = beginning[0] + "#1:" + beginning[1] + tree_list[1]

#    tree_str = tree_str.replace("):%s" % target_branch, ")#1:%s" % target_branch)
    outfile = open("%s/og_%s.tree" % (outdir, orthogroup), 'w')
    outfile.write(new_tree)
    outfile.close()


def augo_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file):
    #trims taxa and adds foreground tags for PAML analysis tree
    ingroup_list = ["AAUR", "APUR"]
    tree = PhyloTree(phylogeny_file)
    tree.prune(taxa_list)
    tree.unroot()
    leaf_list = []
    for leaf in tree:
        leaf_list.append(leaf.name)
    trimmed_ingroup = []
    for species in ingroup_list:
        if species in leaf_list:
            trimmed_ingroup.append(species)
    tree.get_common_ancestor(trimmed_ingroup).add_features(foreground = "foreground")
    tree_str = tree.write(format = 5, features = ["foreground"])
    tree_list = tree_str.split("[&&NHX:foreground=foreground]")
    beginning = tree_list[0].rsplit(":", 1)
    if "clade" in foreground:
        new_tree = beginning[0] + "$1:" + beginning[1] + tree_list[1]
    else:    
        new_tree = beginning[0] + "#1:" + beginning[1] + tree_list[1]
    outfile = open("%s/og_%s.tree" % (outdir, orthogroup), 'w')
    outfile.write(new_tree)
    outfile.close()

def hali_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file):
    #trims taxa and adds foreground tags for PAML analysis tree
    ingroup_list = ["HLIG", "HRUB", "HQUA"]
    tree = PhyloTree(phylogeny_file)
    tree.prune(taxa_list)
    tree.unroot()
    leaf_list = []
    for leaf in tree:
        leaf_list.append(leaf.name)
    trimmed_ingroup = []
    for species in ingroup_list:
        if species in leaf_list:
            trimmed_ingroup.append(species)
    tree.get_common_ancestor(trimmed_ingroup).add_features(foreground = "foreground")
    tree_str = tree.write(format = 5, features = ["foreground"])
    tree_list = tree_str.split("[&&NHX:foreground=foreground]")
    beginning = tree_list[0].rsplit(":", 1)
    if "clade" in foreground:
        new_tree = beginning[0] + "$1:" + beginning[1] + tree_list[1]
    else:
        new_tree = beginning[0] + "#1:" + beginning[1] + tree_list[1]
    outfile = open("%s/og_%s.tree" % (outdir, orthogroup), 'w')
    outfile.write(new_tree)
    outfile.close()

def lasi_branch_fore(taxa_list, foreground, orthogroup, outdir, phylogeny_file):
    #trims taxa and adds foreground tags for PAML analysis tree
    ingroup_list = ["LCAL", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LALB", "LFIG"]
#    ingroup_list = ["LCAL", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LALB", "LFIG"]
    tree = PhyloTree(phylogeny_file)
    tree.prune(taxa_list)
    tree.unroot()
    leaf_list = []
    for leaf in tree:
        leaf_list.append(leaf.name)
    trimmed_ingroup = []
    for species in ingroup_list:
        if species in leaf_list:
            trimmed_ingroup.append(species)
    tree.get_common_ancestor(trimmed_ingroup).add_features(foreground = "foreground")
    tree_str = tree.write(format = 5, features = ["foreground"])
    tree_list = tree_str.split("[&&NHX:foreground=foreground]")
    beginning = tree_list[0].rsplit(":", 1)
    if "clade" in foreground:
        new_tree = beginning[0] + "$1:" + beginning[1] + tree_list[1]
    else:
        new_tree = beginning[0] + "#1:" + beginning[1] + tree_list[1]
    outfile = open("%s/og_%s.tree" % (outdir, orthogroup), 'w')
    outfile.write(new_tree)
    outfile.close()

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
    species_syn_vars = {}
    species_nsyn_vars = {}
    species_fourf_vars = {}
    window_metrics = {}
    for og_num in og_list:
        print og_num
        node_index_dic, node_seqs = read_ancestral_seq("%s/og_%s_working/rst" % (indir, og_num))
        inseq_dic = {}
        reader = SeqIO.parse("%s/og_cds_%s.1.fas" % (aligndir, og_num), format = 'fasta')
        for rec in reader:
#            inseq_dic[rec.id[0:-5]] = str(rec.seq)
            inseq_dic[rec.id[0:4]] = str(rec.seq)
        if len(inseq_dic) < 4:
            continue
        nearest_anc = get_nearest_anc(node_index_dic, inseq_dic)

#        print inseq_dic.keys()
#        print node_seqs
        for species, seq in inseq_dic.items():
            cur_metric = window_metric(seq, node_seqs[nearest_anc[species]])
            syn_pos, nsyn_pos, fourf_pos = sub_positions(seq, node_seqs[nearest_anc[species]])
#            syn_pos, nsyn_pos, fourf_pos = sub_positions(seq, inseq_dic["ECOL"])
            if species not in species_syn_dists.keys():
                species_syn_dists[species] = []
                species_nsyn_dists[species] = []
                species_fourf_dists[species] = []
                species_syn_vars[species] = []
                species_nsyn_vars[species] = []
                species_fourf_vars[species] = []
                window_metrics[species] = []
            if cur_metric:
                window_metrics[species].append(cur_metric)
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
            if len(syn_dists) > 5:
                species_syn_vars[species].append(numpy.var(syn_dists))
            if len(nsyn_dists) > 5:
                species_nsyn_vars[species].append(numpy.var(nsyn_dists))
            if len(fourf_dists) > 5:
                species_fourf_vars[species].append(numpy.var(fourf_dists))
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
        outfile = open("%s/%s_nsyn_vars.txt" % (outdir, species), 'w')
        for var in species_nsyn_vars[species]:
            outfile.write(str(var) + "\n")
        outfile.close()
        outfile = open("%s/%s_syn_vars.txt" % (outdir, species), 'w')
        for var in species_syn_vars[species]:
            outfile.write(str(var) + "\n")
        outfile.close()
        outfile = open("%s/%s_fourf_vars.txt" % (outdir, species), 'w')
        for var in species_fourf_vars[species]:
            outfile.write(str(var) + "\n")
        outfile.close()
        outfile = open("%s/%s_window_metric.txt" % (outdir, species), 'w')
        for met in window_metrics[species]:
            outfile.write(str(met) + "\n")
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

def window_metric(inseq, ancestor):
    window_size = 60
    step_size = 30
    min_diff = 3
    num_windows = 0
    index = 0
    metric_sum = 0
    while index + step_size < len(inseq):
        cur_inseq = inseq[index : index+window_size]
        cur_anc = ancestor[index : index+window_size]
        syn_pos, nsyn_pos, fourf_pos = sub_positions(cur_inseq, cur_anc)
        num_syns = len(syn_pos)
        num_nsyns = len(nsyn_pos)
        num_fourfs = len(fourf_pos)
        if num_syns + num_nsyns >= min_diff:
            num_windows += 1
            metric_sum += abs(num_nsyns - num_syns)
        index += step_size
    if num_windows > 3:
        return metric_sum / num_windows
    else:
        return False



def sub_positions(inseq, ancestor):
    fourfold_list = changes.fourfold_codons()
    empty_chars = ["N", "-", "X"]
    insyns = []
    innsyns = []
    fourf = []
    x = 0
#    print len(inseq)
#    print len(ancestor)
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

def target_taxa_in_og(ortho_dic, target_taxa, og_list):
    new_ogs = []
    for og_num in og_list:
        og_complete = True
        for taxa in target_taxa:
            if taxa not in ortho_dic[og_num]:
                og_complete = False
                break
        if og_complete:
            new_ogs.append(og_num)
    return new_ogs

def make_taxa_dic(infile):
    #This reads in the minimum required taxonomy file. It returns two
    #dictionaries. manda_dic has tuples of taxon names as keys and 
    #the number of those taxa required as values. multi_dic has a tuple
    #of tuples as keys. The idea here is that the nested tuples represent
    #sets of species that must appear together. The values are the number
    #of nested tuples that must appear.
    reader = open(infile, 'rU')
    manda_dic = {}
    multi_dic = {}
    remove_list = []
    for line in reader:
        cur_line = line.split()
        if int(cur_line[0]) < 0:
            remove_list.append(cur_line[1])
        if "(" in cur_line[1]:
            taxa_list = cur_line[1].split("),(")
            first_pair = taxa_list[0].replace(")","").replace("(", "").split(",")
            first_tuple = (first_pair[0],)
            for taxon in first_pair[1:]:
                first_tuple = first_tuple + (taxon,)
            taxa_tuple = (first_tuple,)

            for taxon in taxa_list[1:]:
                next_pair = taxon.replace(")","").replace("(", "").split(",")
                next_tuple = (next_pair[0],)
                for next_taxon in next_pair[1:]:
                    next_tuple = next_tuple + (next_taxon,)
                taxa_tuple = taxa_tuple + (next_tuple,)
            multi_dic[taxa_tuple] = int(cur_line[0])
                
        else:
            taxa_list = cur_line[1].split(",")
            taxa_tuple = (taxa_list[0],)
            for taxon in taxa_list[1:]:
                taxa_tuple = taxa_tuple + (taxon,)
            manda_dic[taxa_tuple] = int(cur_line[0])
    return manda_dic, multi_dic, remove_list
            
def min_taxa_membership(manda_dic, multi_dic, remove_list, index_file, min_taxa, exclude_paras):
    reader = open(index_file, 'rU')
    og_list = []
    for line in reader:
        if line.startswith("#"):
            continue
        
        include = True
        cur_line = line.split()
        if int(cur_line[1]) == 0:
            continue
        taxa_to_remove = []
        taxa_list = [seq[0:4] for seq in cur_line[3].split(",")]
        for taxa in taxa_list:
            if taxa in remove_list:
                taxa_to_remove.append(taxa)
        if int(cur_line[1]) - len(taxa_to_remove) < min_taxa:
            continue
        if int(cur_line[2]) > int(cur_line[1]):
            if int(cur_line[1]) - (int(cur_line[2]) - int(cur_line[1])) - len(remove_list) < min_taxa:
                continue #excludes loci where there aren't enough taxa after paralog removal

        if exclude_paras:
            taxa_noparas = [taxa for taxa in taxa_list if taxa_list.count(taxa)<2]
        else:
            taxa_noparas = list(set(taxa_list))
        for these_taxa, min_count in manda_dic.items():
            mycount = 0
            for this_taxon in these_taxa:
                if this_taxon in taxa_noparas:
                    mycount += 1
            if mycount < min_count:
                include = False
                break
        if not include:
            continue
        for these_sets, min_count in multi_dic.items():
            mycount = 0
            for this_set in these_sets:
                set_total = 0
                for this_taxon in this_set:
                    if this_taxon in taxa_noparas:
                        set_total += 1
                if set_total == len(this_set):
                    mycount += 1
            if mycount < min_count:
                include = False
        if include:
            og_list.append(int(cur_line[0]))
    return og_list

def ncar_min_taxa_membership(manda_dic, multi_dic, remove_list, index_file, min_taxa, exclude_paras):
    reader = open(index_file, 'rU')
    og_list = []
    for line in reader:
        if line.startswith("#"):
            continue
        include = True
        cur_line = line.split()
        if int(cur_line[2]) == 0:
            continue
        taxa_to_remove = []
        taxa_list = [seq[0:4] for seq in cur_line[3].split(",")]
        for taxa in taxa_list:
            if taxa in remove_list:
                taxa_to_remove.append(taxa)
#        if int(cur_line[1]) - len(taxa_to_remove) < min_taxa:
        if int(cur_line[2]) - len(taxa_to_remove) < min_taxa:
            continue
        if exclude_paras:
            taxa_noparas = [taxa for taxa in taxa_list if taxa_list.count(taxa)<2]
        else:
            taxa_noparas = list(set(taxa_list))
        for these_taxa, min_count in manda_dic.items():
            mycount = 0
            for this_taxon in these_taxa:
                if this_taxon in taxa_noparas:
                    mycount += 1
            if mycount < min_count:
                include = False
                break

        if not include:
            continue
        for these_sets, min_count in multi_dic.items():
            mycount = 0
            for this_set in these_sets:
                set_total = 0
                for this_taxon in this_set:
                    if this_taxon in taxa_noparas:
                        set_total += 1
                if set_total == len(this_set):
                    mycount += 1
            if mycount < min_count:
                include = False
        if include:
            og_list.append(cur_line[1])
    return og_list


# def min_taxa_membership(ortho_dic, target_taxa, og_list, min_number):
#     new_ogs = []
#     for og_num in og_list:
#         og_complete = True
#         target_count = 0
#         for taxa in target_taxa:
#             if taxa in ortho_dic[og_num]:
#                 target_count += 1
#         if target_count >= min_number:
#             new_ogs.append(og_num)
#     return new_ogs

def og_taxa_list(ortho_dic, og_num):
    new_ogs = []
    for og_num in og_list:
        og_complete = True
        target_count = 0
        for taxa in target_taxa:
            if taxa in ortho_dic[og_num]:
                target_count += 1
        if target_count >= min_number:
            new_ogs.append(og_num)
    return new_ogs


def make_go_database(ortho_dic, data_species, outpath, gaf_file):
    outfile = open("%s.gaf" % outpath, 'w')
    og_file = open("%s.genelist" % outpath, 'w')
    used_ogs = {}
    used_gos = {}
    for species in data_species:
        print species
        species_dic = {}
        for k, v in ortho_dic.items():
            if species in v.keys():
                species_dic[v[species][0]] = k
#        reader = open("/Genomics/kocherlab/berubin/annotation/interproscan/%s.gaf" % species, 'rU')
#        reader = open("/Genomics/kocherlab/berubin/acacias/selection/gilliamella/interproscan/%s.gaf" % species, 'rU')
#        reader = open("/Genomics/kocherlab/berubin/acacias/selection/bartonella/interproscan/%s.gaf" % species, 'rU')
        reader = open(gaf_file, 'rU')
        for line in reader:
            cur_gene = line.split()[1]
            if cur_gene in species_dic.keys():
                cur_og = species_dic[cur_gene]
                cur_go = line.split()[3]                
                if cur_og not in used_gos.keys():
                    used_gos[cur_og] = []
                    used_ogs[cur_og] = []
                if cur_go not in used_gos[cur_og]:
                    used_gos[cur_og].append(cur_go)
#                if cur_og not in used_ogs.keys():
                    used_ogs[cur_og].append(line.replace(cur_gene, str(cur_og)))
#            cur_gene = line.split()[1][:-3]
#            cur_og = get_og_num(cur_gene, ortho_dic)
#            if cur_og:
#                if cur_og not in used_ogs.keys():
#                    used_ogs[cur_og] = line.replace(cur_gene, str(cur_og))
#                    print cur_og
    for og, outline_list in used_ogs.items():
        for outline in outline_list:
            outfile.write(outline)
        og_file.write("%s\n" % og)
    outfile.close()
    og_file.close()

#population_file is the OGs with human representatives
#pop_condition_file is the OGs with human representatives that are autism-associated
#target_file is OGs evolving faster in social species
#target_back_file is OGs tested 
def hypergeom_test(population_file, pop_condition_file, target_file, target_back_file):
    human_reps = []
    reader = open(population_file, 'rU')
    for line in reader:
        if line.startswith("OG_"):
            line = line.replace("OG_", "")
        cur_line = line.split()
        human_reps.append(cur_line[0])
    human_condition = []
    reader = open(pop_condition_file, 'rU')
    for line in reader:
        if line.startswith("OG_"):
            line = line.replace("OG_", "")
        cur_line = line.split()
        human_condition.append(cur_line[0])
    target_genes = []
    target_condition_genes = []
    reader = open(target_file, 'rU')
    for line in reader:
        if line.startswith("OG_"):
            line = line.replace("OG_", "")
        cur_line = line.split()
        cur_gene = cur_line[0]
        if cur_gene in human_reps:
            target_genes.append(cur_gene)
            if cur_gene in human_condition:
                target_condition_genes.append(cur_gene)
    target_backs = []
    target_condition_backs = []
    reader = open(target_back_file, 'rU')
    for line in reader:
        if line.startswith("OG_"):
            line = line.replace("OG_", "")
        cur_line = line.split()
        cur_gene = cur_line[0]
        if cur_gene in human_reps:
            target_backs.append(cur_gene)
            if cur_gene in human_condition:
                target_condition_backs.append(cur_gene)
    M = len(target_backs)
    n = len(target_condition_backs)
    N = len(target_genes)
    x = len(target_condition_genes)
    pval = scipy.stats.hypergeom.sf(x-1, M, n, N)
    folde = (float(x) / N) / (float(n) / M)
#    print "Total population: %s" % M
#    print "Total population with condition: %s" % n
#    print "Subset of population: %s" % N
#    print "Subset with condition: %s" % x
#    print "Probability this many or more: %s" % pval
#    print "Fold-enrichment: %s" % (folde)
    print "M,n,N,x,p,fe: %s,%s,%s,%s,%s,%s" % (M, n, N, x, pval, folde)
    print "Genes overlapping: %s" % (",".join(target_condition_genes))

# def hypergeom(population_file, pop_condition_file, subset_file, ortho_dic, outlier_num):
#     reader = open(population_file, 'rU')
#     pop_list = []
#     gene_og_dic = {}
#     if os.path.exists("%s_og_index" % population_file):
#         reader = open("%s_og_index" % population_file, 'rU')
#         for line in reader:
#             cur_line = line.split()
#             pop_list.append(int(cur_line[1]))
#             gene_og_dic[cur_line[0]] = int(cur_line[1])
#     else:
#         og_index = open("%s_og_index" % population_file, 'w')
#         for line in reader:
#             cur_gene = line.split()[0]
#             cur_og = int(get_og_num(cur_gene[0:-3], ortho_dic))
#             if cur_og:
#                 pop_list.append(cur_og)
#                 gene_og_dic[cur_gene] = cur_og
#                 og_index.write("%s\t%s\n" % (cur_gene, cur_og))
# #        print cur_og
# #            print cur_gene
#         og_index.close()
#     reader = open(pop_condition_file, 'rU')
#     pop_cond_list = []
#     for line in reader:
#         cur_gene = line.split()[0]
#         if cur_gene in gene_og_dic.keys():
#             pop_cond_list.append(gene_og_dic[cur_gene])
#     reader = open(subset_file, 'rU')
#     subset_list = []
#     for line in reader:
#         if outlier_num > 0:
#             if len(subset_list) > outlier_num:
#                 break
#         if line.startswith("#"):
#             continue
#         if line.startswith("OG\t"):
#             continue
#         if line.startswith("baseMean"):
#             continue
# #        if float(line.split()[1]) > 0.05:
# #            continue
#         if line.startswith("OG_"):
#             cur_og = int(line.split()[0].split("OG_")[1])
#         else:
#             cur_og = int(line.split()[0])
#         if cur_og in pop_list:
#             subset_list.append(cur_og)
#     subset_cond_list = []
#     for sub in subset_list:
#         if sub in pop_cond_list:
#             subset_cond_list.append(sub)
#     print "Total population: %s" % len(pop_list)
#     print "Total population with condition: %s" % len(pop_cond_list)
#     print "Subset of population: %s" % len(subset_list)
#     print "Subset with condition: %s" % len(subset_cond_list)
#     print "Probability this many or more: %s" % scipy.stats.hypergeom.sf(len(subset_cond_list)-1, len(pop_list), len(pop_cond_list), len(subset_list))
#     print "Probability this many or fewer: %s" % scipy.stats.hypergeom.cdf(len(subset_cond_list)+1, len(pop_list), len(pop_cond_list), len(subset_list))
#     print subset_cond_list

def rer_hypergeom_test(population_file, pop_condition_file, rer_file, outlier_num, pcut, fast_or_slow):
    reader = open(population_file, 'rU')
    human_reps = []
    for line in reader:
        if line.startswith("OG_"):
            line = line.replace("OG_", "")
        cur_line = line.split()
        human_reps.append(cur_line[0])
    reader = open(pop_condition_file, 'rU')
    human_condition = []
    for line in reader:
        if line.startswith("OG_"):
            line = line.replace("OG_", "")
        cur_line = line.split()
        human_condition.append(cur_line[0])
    reader = open(rer_file, 'rU')
    target_genes = []
    target_condition_genes = []
    target_backs = []
    target_condition_backs = []
    for line in reader:
        if "NA" in line:
            continue
        if line.startswith("#"):
            continue
        if line.startswith("OG\t"):
            continue
        if line.startswith("baseMean"):
            continue
        if line.startswith("Rho"):
            continue
        cur_line = line.split()
        if line.startswith("OG_"):
            cur_og = cur_line[0].split("OG_")[1]
        else:
            cur_og = cur_line[0]
        cur_dir = float(cur_line[1])
        cur_p = float(cur_line[3])
        if cur_og in human_reps:
            target_backs.append(cur_og)
            if cur_og in human_condition:
                target_condition_backs.append(cur_og)
            if cur_p < pcut:
                if fast_or_slow == "fast" and cur_dir > 0:
                    target_genes.append(cur_og)
                    if cur_og in human_condition:
                        target_condition_genes.append(cur_og)
                elif fast_or_slow == "slow" and cur_dir < 0:
                    target_genes.append(cur_og)
                    if cur_og in human_condition:
                        target_condition_genes.append(cur_og)
    M = len(target_backs)
    n = len(target_condition_backs)
    N = len(target_genes)
    x = len(target_condition_genes)
    pval = scipy.stats.hypergeom.sf(x-1, M, n, N)
    folde = (float(x) / N) / (float(n) / M)
    print "Direction,pcut: %s,%s" % (fast_or_slow, pcut)
#    print "Total population: %s" % M
#    print "Total population with condition: %s" % n
#    print "Subset of population: %s" % N
#    print "Subset with condition: %s" % x
#    print "Probability this many or more: %s" % pval
#    print "Fold-enrichment: %s" % (folde)
#    print ",".join(target_condition_genes)
    print "M,n,N,x,p,fe: %s,%s,%s,%s,%s,%s" % (M, n, N, x, pval, folde)
    print "Genes overlapping: %s" % (",".join(target_condition_genes))
    

def rer_hypergeom(population_file, pop_condition_file, subset_file, ortho_dic, outlier_num, pcut, fast_or_slow, outfile):
    reader = open(population_file, 'rU')
    pop_list = []
    gene_og_dic = {}
    if os.path.exists("%s_og_index" % population_file):
        reader = open("%s_og_index" % population_file, 'rU')
        for line in reader:
            cur_line = line.split()
            pop_list.append(int(cur_line[1]))
            gene_og_dic[cur_line[0]] = int(cur_line[1])
    else:
        og_index = open("%s_og_index" % population_file, 'w')
        for line in reader:
            cur_gene = line.split()[0]

            if cur_gene.startswith("LALB"):
                cur_og = int(get_og_num(cur_gene, ortho_dic))
            else:
                cur_og = int(get_og_num(cur_gene.split("-")[0], ortho_dic))
            if cur_og:
                pop_list.append(cur_og)
                gene_og_dic[cur_gene] = cur_og
                og_index.write("%s\t%s\n" % (cur_gene, cur_og))
#        print cur_og
#            print cur_gene
        og_index.close()
    reader = open(pop_condition_file, 'rU')
    pop_cond_list = []
    for line in reader:
        cur_gene = line.split()[0]
        if cur_gene in gene_og_dic.keys():
            pop_cond_list.append(gene_og_dic[cur_gene])
    reader = open(subset_file, 'rU')
    subset_list = []
    tested_ogs = []
    for line in reader:
        if "NA" in line:
            break
        if outlier_num > 0:
            if len(subset_list) > outlier_num:
                break
        if line.startswith("#"):
            continue
        if line.startswith("OG\t"):
            continue
        if line.startswith("baseMean"):
            continue
        if line.startswith("Rho"):
            continue
#        if float(line.split()[1]) > 0.05:
#            continue
        if line.startswith("OG_"):
            cur_og = int(line.split()[0].split("OG_")[1])
        else:
            cur_og = int(line.split()[0])
        tested_ogs.append(cur_og)
        if cur_og in pop_list:
            cur_line = line.split()
            if float(cur_line[3]) < pcut:
                if fast_or_slow == "fast":
                    if float(cur_line[1]) < 0:
                        continue
                    else:
                        subset_list.append(cur_og)
                elif fast_or_slow == "slow":
                    if float(cur_line[1]) > 0:
                        continue
                    else:
                        subset_list.append(cur_og)
                else:
                    subset_list.append(cur_og)
    subset_cond_list = []
    for sub in subset_list:
        if sub in pop_cond_list:
            subset_cond_list.append(sub)
    tested_pop_list = []
    tested_pop_cond_list = []
    for tested_og in tested_ogs:
        if tested_og in pop_list:
            tested_pop_list.append(tested_og)
        if tested_og in pop_cond_list:
            tested_pop_cond_list.append(tested_og)
        
    outtie = open(outfile, 'w')
    # outtie.write("Total population: %s\n" % len(pop_list))
    # outtie.write("Total population with condition: %s\n" % len(pop_cond_list))
    # outtie.write("Subset of population: %s\n" % len(subset_list))
    # outtie.write("Subset with condition: %s\n" % len(subset_cond_list))
    # outtie.write("Probability this many or more: %s\n" % scipy.stats.hypergeom.sf(len(subset_cond_list)-1, len(pop_list), len(pop_cond_list), len(subset_list)))
    # outtie.write("Probability this many or fewer: %s\n" % scipy.stats.hypergeom.cdf(len(subset_cond_list)+1, len(pop_list), len(pop_cond_list), len(subset_list)))
    # outtie.close()
    outtie.write("Total population: %s\n" % len(tested_pop_list))
    outtie.write("Total population with condition: %s\n" % len(tested_pop_cond_list))
    outtie.write("Subset of population: %s\n" % len(subset_list))
    outtie.write("Subset with condition: %s\n" % len(subset_cond_list))
    outtie.write("Probability this many or more: %s\n" % scipy.stats.hypergeom.sf(len(subset_cond_list)-1, len(tested_pop_list), len(tested_pop_cond_list), len(subset_list)))
    outtie.write("Probability this many or fewer: %s\n" % scipy.stats.hypergeom.cdf(len(subset_cond_list)+1, len(tested_pop_list), len(tested_pop_cond_list), len(subset_list)))
    for og in subset_cond_list:
        outtie.write("%s\n" % (og))
    outtie.close()

def goatools_termfinder(forelist, backlist, database_file, go_dir):
    if not os.path.isdir("%s" % (go_dir)):
        os.mkdir("%s" % go_dir)
    fore_file = open("%s/forelist.txt" % (go_dir), 'w')
    for og_num in forelist:
        fore_file.write("OG_%s\n" % (og_num))
    fore_file.close()
    back_file = open("%s/backlist.txt" % (go_dir), 'w')
    for og_num in backlist:
        back_file.write("OG_%s\n" % (og_num))
    back_file.close()
    cmd = ["/Genomics/kocherlab/berubin/local/src/goatools/goatools/scripts/find_enrichment.py", "%s/forelist.txt" % go_dir, "%s/backlist.txt" % go_dir, database_file, "--obo", "/Genomics/kocherlab/berubin/local/src/goatools/goatools/scripts/go-basic.obo", "--pval", "10", "--outfile", "%s_ps.txt" % (go_dir)]
    subprocess.call(cmd)
    reader = open("%s_ps.txt" % (go_dir), 'rU')
    p_dic = {}
    go_deets = {}
    p_list = []
    go_list = []
    uncorr_dic = {}
    for line in reader:
        if line.startswith("#"):
            continue
        cur_line = line.strip().split("\t")
        if cur_line[1] in ["CC", "MF"]:
            continue
        go_count = int(cur_line[8])
        if go_count > 2:
            cur_p = float(cur_line[6])
            go_deets[cur_line[0]] = cur_line
            p_list.append(cur_p)
            go_list.append(cur_line[0])
            uncorr_dic[cur_line[0]] = cur_p
    pval_corr = smm.multipletests(p_list, alpha = 0.05, method = 'fdr_bh')[1]
    for x in range(len(pval_corr)):
        p_dic[go_list[x]] = pval_corr[x]
    sorted_pvals = sorted(p_dic.items(), key = lambda x: x[1])
    outfile = open("%s_filtered_ps.txt" % go_dir, 'w')
    outfile.write("GOterm\tpval\tFDR_p\tNS\tenrichm\tname\tratio_in_study\tratio_in_pop\tdepth\tstudy_count\tstudy_items\n")
    for pval in sorted_pvals:
        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (pval[0], uncorr_dic[pval[0]], pval[1], "\t".join(go_deets[pval[0]][1:6]), "\t".join(go_deets[pval[0]][7:9]), go_deets[pval[0]][-1]))
    outfile.close()


def run_termfinder(forelist, backlist, database_file, ortho_dic, go_dir):
#    reader = open(backfile, 'rU')
#    back_list = []
#    for line in reader:
#        cur_og = int(line.strip())
#        back_list.append(cur_og)
    if not os.path.isdir("%s" % (go_dir)):
        os.mkdir("%s" % go_dir)

    testing_db = open("%s/testing_db.gaf" % go_dir, 'w')
    reader = open(database_file, 'rU')
    good_backs = []
    for line in reader:
        if int(line.split()[1]) in backlist:
#        if line.split()[1] in backlist:
            #go termfinder fails when an OG is named "0"
            testing_db.write(line.replace("\t0\t0\t", "\t9999999\t9999999\t"))
            if int(line.split()[1]) not in good_backs:
#            if line.split()[1] not in good_backs:
                good_backs.append(int(line.split()[1]))
#                good_backs.append(line.split()[1])
    testing_db.close()

    fore_file = open("%s/forelist.txt" % (go_dir), 'w')
    for og_num in forelist:
        if og_num in good_backs:
            if og_num == 0:
                fore_file.write("9999999\n")
            else:
                fore_file.write("%s\n" % (og_num))
    fore_file.close()
#    cmd = ["/Genomics/kocherlab/berubin/local/src/GO-TermFinder-0.86/examples/analyze.pl", "%s/testing_db.gaf" % go_dir, str(len(good_backs)), "/Genomics/kocherlab/berubin/annotation/interproscan/go-basic.obo", "%s/forelist.txt" % (go_dir)]
    cmd = ["/Genomics/kocherlab/berubin/local/src/GO-TermFinder-0.86/examples/analyze.pl", "%s/testing_db.gaf" % go_dir, str(len(good_backs)), "/Genomics/kocherlab/berubin/annotation/trinotate/go-basic_7_24_18.obo", "%s/forelist.txt" % (go_dir)]
    print cmd
    subprocess.call(cmd)
    reader = open("%s/forelist.txt.terms" % (go_dir), 'rU')
    cp_file = open("%s.txt" % go_dir, 'w')
    pval_file = open("%s_ps.txt" % go_dir, 'w')
    inog = False
    readogs = False
    p_dic = {}
    uncorp_dic = {}
    frac_dic = {}
    term_dic = {}
#    seq_dic = get_prot_seqs()
    for line in reader:
        cp_file.write(line)
        if line.startswith("GOID"):
            cur_go = line.split()[1].replace(":", "_").strip()
            inog = True
        if line.startswith("TERM"):
            cur_term = line.strip().split("\t")[1]
            term_dic[cur_go] = cur_term
        if line.startswith("CORRECTED P-VALUE"):
            if "P" not in line.strip().split()[-1]:
                p_dic[cur_go] = float(line.strip().split()[-1])
            else:
                p_dic[cur_go] = 1.0
        if line.startswith("UNCORRECTED P-VALUE"):
            uncorp_dic[cur_go] = float(line.strip().split()[-1])

        if line.startswith("NUM_ANNOTATIONS"):
            cur_line = line.split()
            inhits = float(cur_line[1])
            inback = float(cur_line[3])
            outhits = float(cur_line[8])
            outback = float(cur_line[10])
            frac_line = "%s\t%s\t%s\t%s\t%s\t%s" % (inhits, inback, outhits, outback, round(inhits/inback, 4), round(outhits/outback, 4))
            frac_dic[cur_go] = frac_line
        if "The genes annotated to this node are:" in line:
            readogs = True
            continue
        if inog and readogs:
            # cur_line = line.strip().split(", ")
            # outfile = open("%s/%s.fasta" % (go_dir, cur_go), 'w')
            # for og in cur_line:
            #     for species in ["AAUR", "LCAL", "LLEU", "LMAL", "HLIG", "HRUB", "APUR", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "LFIG", "AVIR"]:
            #         if species in ortho_dic[int(og)] and "," not in ortho_dic[int(og)][species][0]:
            #             og_seq = seq_dic[species][ortho_dic[int(og)][species][0]]
            #             outfile.write(">%s\n%s\n" % (og, og_seq))
            #             break
            # outfile.close()
            readogs = False
            inog = False
    cp_file.close()
    sorted_ps = sorted(uncorp_dic.items(), key = lambda x: x[1])
    for go_ps in sorted_ps:
        pval_file.write("%s\t%s\t%s\t%s\t%s\n" % (go_ps[0].replace("_", ":"), round(p_dic[go_ps[0]], 8), round(go_ps[1],8), frac_dic[go_ps[0]], term_dic[go_ps[0]]))
    pval_file.close()
        
def rer_termfinder(forefile, backfile, database_file, ortho_dic, go_dir, sig_column, sig_cutoff, fast_or_slow):
    fore_list = rer_external_list(forefile, int(sig_column), float(sig_cutoff), fast_or_slow)
    back_list = external_list(backfile, -1, -1)
    run_termfinder(fore_list, back_list, database_file, ortho_dic, go_dir)

def rer_goatools(forefile, backfile, database_file, go_dir, sig_column, sig_cutoff, fast_or_slow):

    fore_list = rer_external_list(forefile, int(sig_column), float(sig_cutoff), fast_or_slow)
    back_list = external_list(backfile, -1, -1)
    goatools_termfinder(fore_list, back_list, database_file, go_dir)

def og_list_goatools(forefile, backfile, database_file, outdir):
    fore_list = external_list(forefile, -1, -1)
    back_list = external_list(backfile, -1, -1)
    goatools_termfinder(fore_list, back_list, database_file, outdir)

def rer_external_list(forefile, sig_column, sig_cutoff, fast_or_slow):
    reader = open(forefile, 'rU')
    fore_list = []
    for line in reader:
        if line.startswith("OG\t"):
            continue
        if line.startswith("baseMean"):
            continue
        if line.startswith("Rho"):
            continue
        cur_line = line.strip().split()
        if sig_column > -1:
            if cur_line[sig_column] == "NA":
                continue
            #this line is for RER output files. <0 means that it is only looking at faster genes. >0 means only slower
            if fast_or_slow == "fast":
                if float(cur_line[1]) < 0:
                    continue
            if fast_or_slow == "slow":
                if float(cur_line[1]) > 0:
                    continue
            if float(cur_line[sig_column]) < sig_cutoff:
                if line.startswith("OG_"):
                    fore_list.append(int(line.strip().split()[0][3:]))
                else:
                    fore_list.append(int(line.strip().split()[0]))
        else:
            if line.startswith("OG_"):
                fore_list.append(int(line.strip().split()[0][3:]))
            else:
                fore_list.append(int(line.strip().split()[0]))

    return fore_list


def external_list(forefile, sig_column, sig_cutoff):
    reader = open(forefile, 'rU')
    fore_list = []
    for line in reader:
        if line.startswith("OG\t"):
            continue
        if line.startswith("baseMean"):
            continue
        if line.startswith("Rho"):
            continue
        if "NA" in line:
            continue #When RERconverge can't include a locus for whatever reason, it writes the results as NAs. These shouldn't be included because they weren't actually tested.
        cur_line = line.strip().split()
        if sig_column > -1:
            if cur_line[sig_column] == "NA":
                continue
            #this line is for RER output files. <0 means that it is only looking at faster genes. >0 means only slower
#            if float(cur_line[1]) > 0:
#                continue
            if float(cur_line[sig_column]) < sig_cutoff:
                if line.startswith("OG_"):
                    fore_list.append(int(line.strip().split()[0][3:]))
                else:
                    fore_list.append(int(line.strip().split()[0]))
        else:
            if line.startswith("OG_"):
                fore_list.append(int(line.strip().split()[0][3:]))
            elif line.startswith("GB"):
                fore_list.append(line.strip().split()[0])
            else:
                fore_list.append(int(line.strip().split()[0]))
    return fore_list
                       
def og_list_termfinder(forefile, backfile, database_file, ortho_dic, go_dir, sig_column, sig_cutoff):
    fore_list = external_list(forefile, int(sig_column), float(sig_cutoff))
    back_list = external_list(backfile, -1, -1)
    run_termfinder(fore_list, back_list, database_file, ortho_dic, go_dir)

def yn_estimates(og_list,indir, outdir, tree_file, min_taxa, use_gblocks, remove_list):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if not os.path.isdir("%s_results" % outdir):
        os.mkdir("%s_results" % outdir)
    dnds_dic = {}
    dn_dic = {}
    ds_dic = {}
    species_pairs = []
    for cur_og in og_list:
        good_align = prep_paml_files(cur_og, indir, outdir, "yn", tree_file, "yn", min_taxa, use_gblocks, remove_list)
        print cur_og
        if not good_align:
            continue
        if good_align == "too short":
            continue
        try: #had to put this in because sometimes yn just results in nans
            pair_dic = paml_tests.pairwise_yn(cur_og, outdir)
        except IndexError:
            print "OG_%s failed for unknown reasons" % cur_og
            continue
  
#        print pair_dic["SPRA"]["SLEU"]
        dnds_dic[cur_og] = {}
        dn_dic[cur_og] = {}
        ds_dic[cur_og] = {}
        for gene in pair_dic.keys():
            for other_gene in pair_dic[gene]:
                species_tuple = (gene[0:4], other_gene[0:4])
                reverse_tuple = (other_gene[0:4], gene[0:4])
                if species_tuple not in species_pairs and reverse_tuple not in species_pairs:
                    species_pairs.append(species_tuple)
                    cur_tuple = species_tuple
                elif species_tuple in species_pairs:
                    cur_tuple = species_tuple
                elif reverse_tuple in species_pairs:
                    cur_tuple = reverse_tuple
                if float(pair_dic[gene][other_gene]["YN00"]["omega"]) == 99 or float(pair_dic[gene][other_gene]["YN00"]["dN"]) == 0.0: # or float(pair_dic[gene][other_gene]["YN00"]["omega"]) == 0.0:
                    dnds_dic[cur_og][cur_tuple] = "NA"
                    dn_dic[cur_og][cur_tuple] = "NA"
                    ds_dic[cur_og][cur_tuple] = "NA"
                    continue
                if cur_tuple not in dnds_dic[cur_og].keys() and reverse_tuple not in dnds_dic[cur_og].keys():
                    dnds_dic[cur_og][cur_tuple] = pair_dic[gene][other_gene]["YN00"]["omega"]
                    dn_dic[cur_og][cur_tuple] = pair_dic[gene][other_gene]["YN00"]["dN"]
                    ds_dic[cur_og][cur_tuple] = pair_dic[gene][other_gene]["YN00"]["dS"]
    for species_pair in species_pairs:
        outfile = open("%s_results/%s_%s.yn" % (outdir, species_pair[0], species_pair[1]), 'w')
        dnfile = open("%s_results/%s_%s_dn.yn" % (outdir, species_pair[0], species_pair[1]), 'w')
        dsfile = open("%s_results/%s_%s_ds.yn" % (outdir, species_pair[0], species_pair[1]), 'w')
        for og_num in og_list:
            if og_num in dnds_dic.keys():
                if species_pair in dnds_dic[og_num].keys():
                    outfile.write("%s\t%s\n" % (og_num, dnds_dic[og_num][species_pair]))
                    dnfile.write("%s\t%s\n" % (og_num, dn_dic[og_num][species_pair]))
                    dsfile.write("%s\t%s\n" % (og_num, ds_dic[og_num][species_pair]))
        outfile.close()
    return dnds_dic



def gene_flanks(gff_file, genome_dic):
    gene_dic = {}
    for line in gff_file:
        cur_line = line.split()
        if cur_line[2] == "gene":
            cur_name = cur_line[-1].split("Name=")[1][0:-1]
            cur_start = int(cur_line[3])
            cur_end = int(cur_line[4])
            cur_scaf = cur_line[0]
            gene_dic[cur_name] = (cur_scaf, cur_start, cur_end)
    flank_dic = {}
    for gene_name, coords in gene_dic.items():
        scaf_seq = genome_dic[coords[0]]
        start_coord = coords[1] - 2000
        end_coord = coords[2] + 2000
        if coords[1] < 2000:
            start_coord = 0
        if end_coord > len(scaf_seq):
            end_coord = len(scaf_seq)
        scafs = scaf_seq[start_coord : coords[1]] + scaf_seq[coords[2] : end_coord]
        scafs = scafs.replace("N", "")
        flank_dic[gene_name] = scafs
    return flank_dic

def codon_bias(outdir, species_list):
    official_dir = "/Genomics/kocherlab/berubin/official_release"
    if not os.path.isdir("%s/codon_bias" % (outdir)):
        os.mkdir("%s/codon_bias/" % (outdir))
    summary_file = open("%s/codon_bias/codon_bias.summ" % outdir, 'w')
    for target_species in species_list:
        print target_species
        if not os.path.isdir("%s/codon_bias/%s" % (outdir, target_species)):
            os.mkdir("%s/codon_bias/%s" % (outdir, target_species))
        if target_species == "LALB":
            gff_file = open("%s/%s/%s_v3/%s/%s_OGS_v1.0_longest_isoform.gff3" % (official_dir, target_species, target_species, target_species, target_species), 'rU')
            genome_seq = SeqIO.parse("%s/%s/%s_v3/%s/%s_genome_v3.0.fasta" % (official_dir, target_species, target_species, target_species, target_species), format = 'fasta')

            cds_seqs = SeqIO.parse("%s/%s/%s_v3/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (official_dir, target_species, target_species, target_species, target_species), format = 'fasta')
        else:
            gff_file = open("%s/%s/%s_OGS_v1.0_longest_isoform.gff3" % (official_dir, target_species, target_species), 'rU')            
            genome_seq = SeqIO.parse("%s/%s/%s_genome_v1.0.fasta" % (official_dir, target_species, target_species), format = 'fasta')
            cds_seqs = SeqIO.parse("%s/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (official_dir, target_species, target_species), format = 'fasta')
        seq_dic = {}
        for rec in genome_seq:
            seq_dic[rec.id] = str(rec.seq)
        flank_dic = gene_flanks(gff_file, seq_dic)
        ncp_dic = {}
        outfile = open("%s/codon_bias/%s_ncp.txt" % (outdir, target_species), 'w')
        outfile.write("OG\tNcp\n")
        for rec in cds_seqs:
            if len(flank_dic[rec.id[:-3]]) < 1000 or len(rec.seq) < 300:
                continue
            ncp = paml_tests.encprime(rec.id[:-3], rec.seq, flank_dic[rec.id[:-3]], "%s/codon_bias/%s" % (outdir, target_species))
            ncp_dic[rec.id[:-3]] = ncp
            outfile.write("%s\t%s\n" % (rec.id[:-3], ncp))
            outfile.flush()
        outfile.close()
        summary_file.write("%s\t%s\n" % (target_species, numpy.mean(ncp_dic.values())))
        summary_file.flush()
    summary_file.close()
        
def codon_bias(outdir, species_list, og_list, ortho_dic):
    official_dir = "/Genomics/kocherlab/berubin/official_release"
    if not os.path.isdir("%s/codon_bias" % (outdir)):
        os.mkdir("%s/codon_bias/" % (outdir))
    summary_file = open("%s/codon_bias/codon_bias.summ" % outdir, 'w')
    for target_species in species_list:
        print target_species
        if not os.path.isdir("%s/codon_bias/%s" % (outdir, target_species)):
            os.mkdir("%s/codon_bias/%s" % (outdir, target_species))
        if target_species == "LALB":
            gff_file = open("%s/%s/%s_v3/%s/%s_OGS_v1.0_longest_isoform.gff3" % (official_dir, target_species, target_species, target_species, target_species), 'rU')
            genome_seq = SeqIO.parse("%s/%s/%s_v3/%s/%s_genome_v3.0.fasta" % (official_dir, target_species, target_species, target_species, target_species), format = 'fasta')

            cds_seqs = SeqIO.parse("%s/%s/%s_v3/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (official_dir, target_species, target_species, target_species, target_species), format = 'fasta')
        else:
            gff_file = open("%s/%s/%s_OGS_v1.0_longest_isoform.gff3" % (official_dir, target_species, target_species), 'rU')            
            genome_seq = SeqIO.parse("%s/%s/%s_genome_v1.0.fasta" % (official_dir, target_species, target_species), format = 'fasta')
            cds_seqs = SeqIO.parse("%s/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (official_dir, target_species, target_species), format = 'fasta')
        seq_dic = {}
        for rec in genome_seq:
            seq_dic[rec.id] = str(rec.seq)
        flank_dic = gene_flanks(gff_file, seq_dic)
        ncp_dic = {}
        outfile = open("%s/codon_bias/%s_ncp.txt" % (outdir, target_species), 'w')
        outfile.write("OG\tNcp\n")
        og_outfile = open("%s/codon_bias/%s_ncp_ogs.txt" % (outdir, target_species), 'w')
        og_outfile.write("OG\tNcp\n")
        for rec in cds_seqs:
            if len(flank_dic[rec.id[:-3]]) < 1000 or len(rec.seq) < 300:
                continue
            ncp = paml_tests.encprime(rec.id[:-3], rec.seq, flank_dic[rec.id[:-3]], "%s/codon_bias/%s" % (outdir, target_species))
            ncp_dic[rec.id[:-3]] = ncp
            outfile.write("%s\t%s\n" % (rec.id[:-3], ncp))
            og_num = get_og_num(rec.id[:-3], ortho_dic)
            if og_num in og_list:
                og_outfile.write("%s\t%s\n" % (rec.id[:-3], ncp))
#            outfile.flush()
#            og_outfile.flush()
        outfile.close()
        og_outfile.close()
        summary_file.write("%s\t%s\n" % (target_species, numpy.mean(ncp_dic.values())))
        summary_file.flush()
    summary_file.close()
        
def read_hyphy_relax(og_list, indir, outdir, forestring):
    outfile = open("%s/compiled_relax_%s.txt" % (outdir, forestring), 'w')
    outfile.write("OG\tK\tP\tFDR_P\tinterpretation\n")
    k_list = []
    p_list = []
    for cur_og in og_list:
#        print cur_og
        if not os.path.exists("%s/og_%s_relax_unlabeledback.txt" % (indir, cur_og)):
            print "%s/og_%s_relax_unlabeledback.txt is missing" % (indir, cur_og)
            continue
        og_file =open("%s/og_%s_relax_unlabeledback.txt" % (indir, cur_og), 'rU')
        interpret_line = False
        for line in og_file:
            if line.startswith("* Relaxation/intensification parameter"):
                cur_k = float(line.strip().split()[-1])
            elif line.startswith("Likelihood ratio test"):
                cur_p = float(line.strip().replace("**.", "").split()[-1])
                interpret_line = True
                p_list.append(cur_p)
                continue
            if interpret_line:
                interpretation = line
                if interpretation.startswith(">Evidence for"):
                    interpretation = interpretation.split("*")[1].split(" ")[0]
                else:
                    interpretation = "nonsignificant"
                k_list.append([cur_og, cur_k, cur_p, interpretation])
                break
        if interpret_line == False:
            print "missing %s" % cur_og
    pval_corr = smm.multipletests(p_list, alpha = 0.1, method = 'fdr_bh')[1]
    for x in range(len(pval_corr)):
        k_list[x].insert(3, pval_corr[x])
    k_list.sort(key = lambda x: x[-2])
    for locus in k_list:
        outfile.write("\t".join([str(i) for i in locus]) + "\n")
    outfile.close()

def read_hyphy_absrel(og_list, indir, outdir):
    outfile = open("%s/compiled_absrel_%s.txt" % (outdir, "all"), 'w')
    outfile.write("OG\tP\tFDR_P\talltests_P\tomega\tprop\tnode\n")
    uncorr_dic = {}
    k_list = []
    p_list = []
    for cur_og in og_list:
        print cur_og
        if not os.path.exists("%s/og_cds_%s.afa.ABSREL.json" % (indir, cur_og)):
            print "%s/og_cds_%s.afa.ABSREL.json is missing" % (indir, cur_og)
            continue
        cur_p_dic = {}
        with open("%s/og_cds_%s.afa.ABSREL.json" % (indir, cur_og), 'rU') as readfile:
            data = json.load(readfile)
        tree = PhyloTree(str(data["input"]["trees"]["0"]) + ";", format = 1)
        for testname, testinfo in data["branch attributes"]["0"].items():
            corr_p = testinfo["Corrected P-value"]
            uncorr_p = testinfo["Uncorrected P-value"]
            omega = testinfo["Baseline MG94xREV omega ratio"]
            adaptive_prop = 0
            p_list.append(float(uncorr_p))
            if testinfo["Rate classes"] == 2:
                adaptive_prop = testinfo["Rate Distributions"][1][1]
            if "Node" in testname:
                daughters = []
                node = tree.search_nodes(name = testname)[0]
                for leaf in node.get_leaves():
                    daughters.append(leaf.name)
                cur_p_dic[",".join(daughters)] = uncorr_p
#                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (cur_og, uncorr_p, corr_p, omega, adaptive_prop, ",".join(daughters)))
                k_list.append([cur_og, uncorr_p, corr_p, omega, adaptive_prop, ",".join(daughters)])
            else:
                cur_p_dic[testname] = uncorr_p
#                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (cur_og, uncorr_p, corr_p, omega, adaptive_prop, testname))
                k_list.append([cur_og, uncorr_p, corr_p, omega, adaptive_prop, testname])

        uncorr_dic[cur_og] = cur_p_dic

    pval_corr = smm.multipletests(p_list, alpha = 0.1, method = 'fdr_bh')[1]
    for x in range(len(pval_corr)):
        k_list[x].insert(3, pval_corr[x])
    k_list.sort(key = lambda x: x[1])
    for locus in k_list:
        outfile.write("\t".join([str(i) for i in locus]) + "\n")
    outfile.close()

def read_hyphy_noncoding(ncar_list, indir, outdir):
    outfile = open("%s/compiled_hyphynoncoding_all.txt" % (outdir), 'w')
    outfile.write("NCAR\tP\tNCARcorr_P\talltests_P\tnull_zeta\talt_zeta_back\talt_zeta_fore\tnode\n")
    line_list = []
    p_list = []
    neg_list = []
    for cur_ncar in ncar_list:
        print cur_ncar
        if not os.path.exists("%s/%s_noncodingselection_summary.txt" % (indir, cur_ncar)):
            print "%s/%s_noncodingselection_summary.txt is missing" % (indir, cur_ncar)
            continue
        ncar_file = open("%s/%s_noncodingselection_summary.txt" % (indir, cur_ncar), 'rU')
        for line in ncar_file:
            cur_line = line.split()
            if len(cur_line) < 2:
                continue
            if "," in cur_line[1]:
                continue
            line_list.append(cur_line)
            p_list.append(float(cur_line[1]))

#            if float(cur_line[1]) < 0:
#                neg_list.append(cur_line)
#            else:

    pval_corr = smm.multipletests(p_list, alpha = 0.1, method = 'fdr_bh')[1]
    for x in range(len(pval_corr)):
        line_list[x].insert(3, pval_corr[x])
    line_list.sort(key = lambda x: x[1])
    for locus in line_list:
        outfile.write("\t".join([str(i) for i in locus]) + "\n")
#    for locus in neg_list:
#        outfile.write("\t".join([str(i) for i in locus]) + "\n")
    outfile.close()


def read_frees_subset(indir, outdir, database_file, go_dir, get_dn_ds, time_tree, og_list, subset_list):
    #reads free ratios files and gets dn/ds ratios
    #can be easily extended to get dn and ds but those are low quality
    soc_sol_pairs = {"LMAR":"LFIG", "LZEP":"LVIE", "LPAU":"LOEN", "AAUR":"APUR"}
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    og_dnds_dic = {}
    og_ds_dic = {}
    og_dn_dic = {}
    species_list = []
    for og_file in glob("%s/og_*.alt" % (indir)):        
        cur_og = int(og_file.split("og_")[1].split(".alt")[0])
        if cur_og not in og_list:
            continue
        reader = open(og_file, 'rU')
        dnds_tree = False
        ds_tree = False
        dn_tree = False
        ds_dic = {}
        dn_dic = {}
        dnds_dic = {}
        first_line = True
        for line in reader:
            if first_line:
                seq_len = int(line.strip().split()[1])
                first_line = False
            if dnds_tree:
                dnds = PhyloTree(line.strip().replace("#", ":"))
                dnds.prune(subset_list, preserve_branch_length = True)
                for leaf in dnds:
                    if leaf.dist < 4 and leaf.dist > 0.0001:
                        if ds_dic[leaf.name] > 0.001 and ds_dic[leaf.name] < 10:
                            if dn_dic[leaf.name] > 0.001 and dn_dic[leaf.name] < 10:
                                dnds_dic[leaf.name] = leaf.dist
                    if leaf.name not in species_list:
                        species_list.append(leaf.name)
                dnds_tree = False
            if ds_tree:
                ds = PhyloTree(line.strip())
                ds.prune(subset_list, preserve_branch_length = True)
                for leaf in ds:
                    ds_dic[leaf.name] = leaf.dist
                ds_tree = False
            if dn_tree:
                dn = PhyloTree(line.strip())
                dn.prune(subset_list, preserve_branch_length = True)
                for leaf in dn:
                    dn_dic[leaf.name] = leaf.dist
                dn_tree = False
            if line.strip() == "dS tree:":
                ds_tree = True
                continue
            if line.strip() == "dN tree:":
                dn_tree = True
                continue
            if line.strip() == "w ratios as labels for TreeView:":
                dnds_tree = True
                continue
        if seq_len < 300:
            continue
        og_dnds_dic[cur_og] = dnds_dic
        if get_dn_ds:
            taxa_pres = dn_dic.keys()
            if len(taxa_pres) == 0:
                continue
            if "Dnov" in taxa_pres:
                taxa_pres.remove("Dnov")
            if "Nmel" in taxa_pres:
                taxa_pres.remove("Nmel")
            tree = PhyloTree(time_tree)
            tree.prune(taxa_pres, preserve_branch_length = True)
            dn_time_dic = standardize_to_time(dn_dic, tree)
            taxa_pres = ds_dic.keys()
            if "Dnov" in taxa_pres:
                taxa_pres.remove("Dnov")
            if "Nmel" in taxa_pres:
                taxa_pres.remove("Nmel")
            tree = PhyloTree(time_tree)
            tree.prune(taxa_pres, preserve_branch_length = True)
            ds_time_dic = standardize_to_time(ds_dic, tree)
            og_ds_dic[cur_og] = ds_time_dic
            og_dn_dic[cur_og] = dn_time_dic
#            og_ds_dic[cur_og] = ds_dic
#            og_dn_dic[cur_og] = dn_dic
    sum_file = open("%s/dnds_summary.txt" % outdir, 'w')
    for species in species_list:
        outfile = open("%s/%s_free_dnds.txt" % (outdir, species), 'w')
        dnds_list = []
        for og in og_dnds_dic.keys():
            if species in og_dnds_dic[og].keys():
                outfile.write("%s\t%s\n" % (og, og_dnds_dic[og][species]))
                dnds_list.append(og_dnds_dic[og][species])
            else:
                outfile.write("%s\tNA\n" % (og))
            
        outfile.close()
        sum_file.write("%s\t%s\t%s\n" % (species, numpy.mean(dnds_list), numpy.var(dnds_list)))
    if get_dn_ds:
        sum_file = open("%s/dn_ds_summary.txt" % outdir, 'w')
        for species in species_list:
            dn_list = []
            ds_list = []
            dnfile = open("%s/%s_free_dn.txt" % (outdir, species), 'w')
            dsfile = open("%s/%s_free_ds.txt" % (outdir, species), 'w')
            for og in og_dn_dic.keys():
                if species in og_dn_dic[og].keys():
                    if og_dn_dic[og][species] < 1:
                        dnfile.write("%s\t%s\n" % (og, og_dn_dic[og][species]))
                        dn_list.append(og_dn_dic[og][species])
                    else:
                        dnfile.write("%s\tNA\n" % (og))
                    if og_ds_dic[og][species] < 1:
                        dsfile.write("%s\t%s\n" % (og, og_ds_dic[og][species]))
                        ds_list.append(og_ds_dic[og][species])
                    else:
                        dsfile.write("%s\tNA\n" % (og))
                else:
                    dnfile.write("%s\tNA\n" % (og))
                    dsfile.write("%s\tNA\n" % (og))
            sum_file.write("%s\t%s\t%s\t%s\t%s\n" % (species, numpy.mean(dn_list), numpy.var(dn_list), numpy.mean(ds_list), numpy.var(ds_list)))
            dnfile.close()
            dsfile.close()
        sum_file.close()

    soc_larger_file = open("%s/soc_larger.txt" % outdir, 'w')
    sol_larger_file = open("%s/sol_larger.txt" % outdir, 'w')
    soc_list = []
    sol_list = []
    background_ogs = []
    for og in og_dnds_dic.keys():
        if og not in og_list:
            continue
        sol_larger = True
        soc_larger = True
        background_og = True
        for soc, sol in soc_sol_pairs.items():
            if soc in og_dnds_dic[og].keys() and sol in og_dnds_dic[og].keys():
                if og_dnds_dic[og][soc] <= og_dnds_dic[og][sol]:
                    soc_larger = False
                if og_dnds_dic[og][sol] <= og_dnds_dic[og][soc]:
                    sol_larger = False
            else:
                soc_larger = False
                sol_larger = False
                background_og = False
        if soc_larger:
            soc_larger_file.write("%s\n" % og)
            soc_list.append(og)
        if sol_larger:
            sol_larger_file.write("%s\n" % og)
            sol_list.append(og)
        if background_og:
            background_ogs.append(og)
    print len(background_ogs)
    soc_larger_file.close()
    sol_larger_file.close()
#    run_termfinder(soc_list, background_ogs, database_file, ortho_dic, "%s_soc" % go_dir)
#    run_termfinder(sol_list, background_ogs, database_file, ortho_dic, "%s_sol" % go_dir)
    sum_file.close()

def gc_content(og_list, indir, outdir):
    gc_count_dic = {}
    len_dic = {}
    for cur_og in og_list:
        reader = SeqIO.parse("%s/og_cds_%s.afa-gb" % (indir, cur_og), format = 'fasta')
        for rec in reader:
            cur_sample = rec.id[0:4]
            if cur_sample not in gc_count_dic.keys():
                gc_count_dic[cur_sample] = 0.0
                len_dic[cur_sample] = 0.0
            gc_count_dic[cur_sample] += str(rec.seq).count("G") + str(rec.seq).count("C")
            len_dic[cur_sample] += len(rec.seq)
    outfile = open("%s/mean_gc.txt" % outdir, 'w')
    for species, gc_count in gc_count_dic.items():
        outfile.write("%s\t%s\n" % (species, gc_count / len_dic[species]))
    outfile.close()
