from Bio.Align.Applications import MuscleCommandline
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
#from Bio.Seq import Seq
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
        
def read_species_pickle(target_species):
    pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species), 'rb')
    gene_objects = pickle.load(pickle_file)
    pickle_file.close()
    return gene_objects

def get_species_data(target_species, num_threads):
    #This is the big method. Harvests gene coordinates from GFF3 files
    #and creates Gene objects with all of their characteristics.
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
    reader = SeqIO.parse("%s/%s/%s_genome_v1.0.fasta" % (official_dir, target_species, target_species), format = 'fasta')
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    cds_dic = {}
    reader = SeqIO.parse("%s/%s/%s_OGS_v1.0_longest_isoform.cds.fasta" % (official_dir, target_species, target_species), format = 'fasta')

    for rec in reader:
        cds_dic[rec.id[:-3]] = str(rec.seq)

#    gff_file = gffParser(open("%s/%s/%s_testset.gff3" % (official_dir, target_species, target_species), 'rU'))
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
#        gene_objects[gene_name] = Gene(gene_name, gene_dic[gene_name][0], gene_dic[gene_name][1])
###
#        print harvest_worker([gene_name, gene_dic, gff_file, seq_dic, target_species])
#        work_queue.put([gene_name, gene_dic, gff_file, seq_dic, target_species])
        work_list.append([gene_name, gene_dic, gff_file, seq_dic, target_species])
#        work_list.append((gene_name, gene_dic, target_species))
#        name_list.append(gene_name)
#        dic_list.append(gene_dic)
#        gff_list.append(gff_file)
#        seq_list.append(seq_dic)
#        species_list.append(target_species)
#        print work_queue.qsize()
#        print work_queue.empty()
#        print gene_name
#    print work_list[0]
#    import itertools
#    params = itertools.izip(name_list, dic_list, gff_list, seq_list, species_list)
#    print len(params)
#    print len(params[0])
    gene_list = pool.map(harvest_worker, work_list)
    print "gene data gotten"
    print str(datetime.datetime.now())
    pool.terminate()
#    print gene_list
    print len(gene_list)
    for gene in gene_list:
        gene_objects[gene.name] = gene
    print "Dumping %s pickle" % target_species
    pickle_file = open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (target_species, target_species), 'wb')
    pickle.dump(gene_objects, pickle_file)
    pickle_file.close()
#    return gene_objects

#    gene_list = pool.map(harvest_worker, itertools.izip(name_list, dic_list, gff_list, seq_list, species_list))
"""
    jobs = []
    print work_queue.empty()
    print work_queue.qsize()
    for i in range(num_threads):
        print i
        worker = Worker(work_queue, result_queue, harvest_worker)
        jobs.append(worker)
        worker.start()
#        print worker
#        print worker.work_queue.qsize()
    print work_queue.qsize()
    print len(jobs)
    for j in jobs:
        print j
        print j.work_queue.qsize()
        j.join()
    print result_queue.qsize()
#    while result_queue.qsize() > 0:
    while not result_queue.empty():
        cur_gene = result_queue.get()
        gene_objects[cur_gene.name] = cur_gene
"""


###
"""
    out_dic = {}
    print "get coords"
    print str(datetime.datetime.now())
    gene_objects = get_gene_coords(gff_file, gene_dic, gene_objects, seq_dic)
    print "get cds dic"
    print str(datetime.datetime.now())
    gene_objects = get_cds_dic(gff_file, gene_dic, gene_objects)
    print "get 3prime"
    print str(datetime.datetime.now())
    gene_objects = get_utr_dic(gff_file, gene_dic, "three", gene_objects)
    print "get 5prime"
    print str(datetime.datetime.now())
    gene_objects = get_utr_dic(gff_file, gene_dic, "five", gene_objects)
    print "get cds sequences"
    print str(datetime.datetime.now())
#    make_cds_sequences(gene_objects, cds_dic)
    gene_objects = make_cds_sequences(gene_objects, num_threads)
    print "done reading gff"
    print str(datetime.datetime.now())
    return gene_objects
"""
def harvest_worker(attribute_list):#gene_name, gene_dic, gff_file, seq_dic, target_species):
    gene_name = attribute_list[0]
    gene_dic = attribute_list[1]
    gff_file = attribute_list[2]
    seq_dic = attribute_list[3]
    target_species = attribute_list[4]
    vcf_reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_filtered_miss.vcf.gz" % (target_species, target_species))
    cur_gene = Gene(gene_name, gene_dic[gene_name][0], gene_dic[gene_name][1])
    gene_deets = gff_file.getGene(gene_dic[gene_name][0], gene_name)[0]
    cur_gene.set_start(gene_deets["start"])
    cur_gene.set_end(gene_deets["end"])
    cur_gene.set_sequence(seq_dic[gene_dic[gene_name][0]], 3000)
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
    cur_gene.get_flank_sequence(3000)
    cur_gene.get_utr_sequence()
    cur_gene.get_intron_sequence()
    cur_gene.get_genotypes(vcf_reader, 3000)
 #   print cur_gene.name
    return cur_gene



"""
def make_cds_sequences(gene_objects, cds_dic):
    #Add sequence data to all of the Gene objects
    for gene in gene_objects.values():
        print gene.name
        gene.get_cds_sequence(cds_dic)
        gene.get_flank_sequence(3000)
        gene.get_utr_sequence()
        gene.get_intron_sequence()
"""
def make_cds_sequences(gene_objects, num_threads):
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()
    for gene in gene_objects.values():
        work_queue.put([gene])
    jobs = []
    for i in range(num_threads):
        worker = Worker(work_queue, result_queue, cds_sequence_worker)
        jobs.append(worker)
        worker.start()
    for j in jobs:
        j.join()
    output_list = []
    print result_queue.qsize()
    while result_queue.qsize() > 0:
        output_list.append(result_queue.get())
        print output_list
    gene_dic = {}
    for gene in output_list:
        print gene.name
        gene_dic[gene.name] = gene
    return gene_dic
        
def cds_sequence_worker(gene):
    gene.get_cds_sequence()
    gene.get_flank_sequence(3000)
    gene.get_utr_sequence()
    gene.get_intron_sequence()
#    print gene.name
    return gene    

def get_gene_coords(gff_file, gene_dic, gene_objects, seq_dic):
    #Add coordinate data for all of the Gene objects
    for gene_name in gene_dic.keys():
        gene_deets = gff_file.getGene(gene_dic[gene_name][0], gene_name)[0]
        gene_objects[gene_name].set_start(gene_deets["start"])
        gene_objects[gene_name].set_end(gene_deets["end"])
        gene_objects[gene_name].set_sequence(seq_dic[gene_dic[gene_name][0]], 3000)
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

def gene_vcf_dic(species, num_threads):
    #This launches the construction of Gene objects and pickles them.
    if os.path.exists("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (species, species)):
        print "Reading %s pickle" % species
        gene_dic = pickle.load(open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (species, species), 'rb'))
    else:
        print "Building %s pickle" % species
        gene_dic = get_species_data(species, num_threads)
        print "make vcf reader"
        print str(datetime.datetime.now())
#        reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_testset.vcf.gz" % (species, species))
        reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_filtered_miss.vcf.gz" % (species, species))
        gene_count = 0
        for gene_name, gene_object in gene_dic.items():
#            if gene_name == "LLEU_00225":
#            if gene_name == "LMAL_09600":
            print gene_name
            print gene_count
            gene_count += 1
            gene_object.get_genotypes(reader, 3000)
        print "Dumping %s pickle" % species
        pickle.dump(gene_dic, open("/scratch/tmp/berubin/resequencing/%s/genotyping/%s_gene_vcf_dic.pickle" % (species, species), 'wb'))
    return gene_dic    

def hka_test(inspecies, outspecies, seq_type, ortho_dic, out_path, num_threads):
    get_species_data(inspecies, num_threads)
    get_species_data(outspecies, num_threads)
#    in_gene_dic = gene_vcf_dic(inspecies, num_threads)
#    out_gene_dic = gene_vcf_dic(outspecies, num_threads)
    og_map = open("%s/%s_%s_og_map.txt" % (out_path, inspecies, outspecies), 'w')
    hka_table = open("%s/%s_%s_%s_hka_table.txt" % (out_path, inspecies, outspecies, seq_type), 'w')
    hka_line_list = []
    numloci = 0
    og_list = []
#    work_queue = multiprocessing.Queue()
#    result_queue = multiprocessing.Queue()
    for og_num in ortho_dic.keys():
        if inspecies not in ortho_dic[og_num] or outspecies not in ortho_dic[og_num]:
            continue
        if ortho_dic[og_num][inspecies][0].count(inspecies) > 1 or ortho_dic[og_num][outspecies][0].count(outspecies) > 1:
            continue
        if len(ortho_dic[og_num][inspecies]) == 1 and len(ortho_dic[og_num][outspecies]):
            og_list.append(og_num)
#    og_list = [2110]
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    in_gene_dic = read_species_pickle(inspecies)
    out_gene_dic = read_species_pickle(outspecies)
    for og in og_list:
        inspecies_gene = ortho_dic[og][inspecies][0][:-3]
        outspecies_gene = ortho_dic[og][outspecies][0][:-3]
        in_gene = in_gene_dic[inspecies_gene]
        out_gene = out_gene_dic[outspecies_gene]
        og_map.write("OG_%s\t%s\t%s\n" % (og, inspecies_gene, outspecies_gene))
#        work_queue.put([inspecies, outspecies, seq_type, og, in_gene, out_gene])
        work_list.append([inspecies, outspecies, seq_type, og, in_gene, out_gene])
    og_map.close()
    hka_lines = pool.map(hka_worker, work_list)
    pool.terminate()
    for line in hka_lines:
        if "too_short" not in line:
            hka_table.write(line)
    hka_table.close()
#    jobs = []
#    for i in range(num_threads):
#        worker = Worker(work_queue, result_queue, hka_worker)
#        jobs.append(worker)
#        worker.start()
#    for j in jobs:
#        j.join()
#        """
#    try:
#        for j in jobs:
#            j.join()
#    except KeyboardInterrupt:
#        for j in jobs:
#            j.terminate()
#            j.join()
#        """
#    output_list = []
#    while result_queue.qsize() > 0:
        
#        output_list.append(result_queue.get())
#    for line in output_list:
#        hka_table.write(line)
#    hka_table.close()
    
#    for result_str in result_queue:
#        hka_table.write(result_str)

#def hka_worker(inspecies, outspecies, seq_type, og_num, in_gene, out_gene):

def hka_worker(attribute_list):
    inspecies = attribute_list[0]
    outspecies = attribute_list[1]
    seq_type = attribute_list[2]
    og_num = attribute_list[3]
    in_gene = attribute_list[4]
    out_gene = attribute_list[5]

#    outfile = open("/Genomics/kocherlab/berubin/local/developing/selection_pipeline/tempfile.tmp", 'w')
#    print in_gene.name
#    print out_gene.name
#    print str(datetime.datetime.now())
    inpoly, outpoly, inseq, outseq, insample, outsample = gather_hka_data(in_gene, out_gene, seq_type)
#    print "data gathered"
#    print str(datetime.datetime.now())
    if len(inseq) < 200 or len(outseq) < 200:
        return "OG_%s\ttoo_short\n" % (og_num)
    if in_gene.strand == -1:
        inseq = str(Seq.Seq(inseq).reverse_complement())
    if out_gene.strand == -1:
        outseq = str(Seq.Seq(outseq).reverse_complement())
    if len(outseq) < len(inseq):
        inseq = inseq[0:len(outseq)]
    if len(inseq) < len(outseq):
        outseq = outseq[0:len(inseq)]
#    print "start aligning"
#    print str(datetime.datetime.now())
    average_diff, align_len = muscle_pairwise_diff_count(inseq, outseq, inspecies, outspecies, og_num)
    if 1.0 * average_diff / align_len >0.15:
        return "OG_%s\ttoo_different\t%s\n" % (og_num, 1.0 * average_diff / align_len)
#    print "done aligning"
#    print str(datetime.datetime.now())

#    print >>outfile, "OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, 1.0, len(inseq), len(outseq), align_len, insample*2, outsample*2, inpoly, outpoly, average_diff)
#    outfile.write("OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, 1.0, len(inseq), len(outseq), align_len, insample*2, outsample*2, inpoly, outpoly, average_diff))
#    outfile.close()
#    return og_num
    return "OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, 1.0, len(inseq), len(outseq), align_len, insample*2, outsample*2, inpoly, outpoly, average_diff)

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
        inseq = "".join(in_gene_dic[inspecies_gene].intron_dic.values())
        outseq = "".join(out_gene_dic[outspecies_gene].intro_dic.values())
        insample = in_gene_dic[inspecies_gene].average_sample_size(in_gene_dic[inspecies_gene].intron_dic)
        outsample = out_gene_dic[inspecies_gene].average_sample_size(out_gene_dic[outspecies_gene].intron_dic)
    elif seq_type == "three_prime_utr":
        inpoly = in_gene_dic[inspecies_gene].utr3_subs
        outpoly = out_gene_dic[outspecies_gene].utr3_subs
        inseq = "".join(in_gene_dic[inspecies_gene].utr3_dic.values())
        outseq = "".join(out_gene_dic[outspecies_gene].utr3_dic.values())
        insample = in_gene_dic[inspecies_gene].average_sample_size(in_gene_dic[inspecies_gene].utr3_dic)
        outsample = out_gene_dic[inspecies_gene].average_sample_size(out_gene_dic[outspecies_gene].utr3_dic)
    elif seq_type == "five_prime_utr":
        inpoly = in_gene_dic[inspecies_gene].utr5_subs
        outpoly = out_gene_dic[outspecies_gene].utr5_subs
        inseq = "".join(in_gene_dic[inspecies_gene].utr5_dic.values())
        outseq = "".join(out_gene_dic[outspecies_gene].utr5_dic.values())
        insample = in_gene_dic[inspecies_gene].average_sample_size(in_gene_dic[inspecies_gene].utr5_dic)
        outsample = out_gene_dic[inspecies_gene].average_sample_size(out_gene_dic[outspecies_gene].utr5_dic)

    return inpoly, outpoly, inseq, outseq, insample, outsample    
"""
def hka_test(inspecies, outspecies, seq_type, ortho_dic, out_path):
    in_gene_dic = gene_vcf_dic(inspecies)
    out_gene_dic = gene_vcf_dic(outspecies)
    hka_table = open("%s/%s_%s_%s_hka_table.txt" % (out_path, inspecies, outspecies, seq_type), 'w')
    hka_line_list = []
    numloci = 0
    test_table = open("%s/%s_%s_%s_hka_table_test.txt" % (out_path, inspecies, outspecies, seq_type), 'w')
    for og_num in ortho_dic.keys():
        if inspecies not in ortho_dic[og_num] or outspecies not in ortho_dic[og_num]:
            continue
        if ortho_dic[og_num][inspecies][0].count(inspecies) > 1 or ortho_dic[og_num][outspecies][0].count(outspecies) > 1:
            continue
        if len(ortho_dic[og_num][inspecies]) == 1 and len(ortho_dic[og_num][outspecies]):
#            if og_num != 2110:
#                continue
            print og_num
            inspecies_gene = ortho_dic[og_num][inspecies][0][:-3]
            outspecies_gene = ortho_dic[og_num][outspecies][0][:-3]
            print inspecies_gene
            print outspecies_gene
            print str(datetime.datetime.now())
            inpoly, outpoly, inseq, outseq, insample, outsample = gather_hka_data(in_gene_dic, out_gene_dic, inspecies_gene, outspecies_gene, seq_type)
            print "data gathered"
            print str(datetime.datetime.now())
            if len(inseq) < 200 or len(outseq) < 200:
                continue
            if in_gene_dic[inspecies_gene].strand == -1:
                inseq = str(Seq.Seq(inseq).reverse_complement())
            if out_gene_dic[outspecies_gene].strand == -1:
                outseq = str(Seq.Seq(outseq).reverse_complement())
            print "start aligning"
            print str(datetime.datetime.now())
            average_diff, align_len = muscle_pairwise_diff_count(inseq, outseq)
            print "done aligning"
            print str(datetime.datetime.now())
        numloci += 1
        hka_line_list.append("OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, 1.0, len(inseq), len(outseq), align_len, insample*2, outsample*2, inpoly, outpoly, average_diff))
        test_table.write("OG_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (og_num, 1.0, len(inseq), len(outseq), align_len, insample*2, outsample*2, inpoly, outpoly, average_diff))
    hka_table.write("HKA table\n")
    hka_table.write("#%s against %s\n" % (inspecies, outspecies))
    hka_table.write("%s\n" % numloci)
    hka_table.write("%s\t%s\n" % (inspecies, outspecies))
    for line in hka_line_list:
        hka_table.write(line)
        hka_table.flush()
    hka_table.close()
    #run this: https://bio.cst.temple.edu/~hey/program_files/HKA/HKA_Documentation.htm

def gather_hka_data(in_gene_dic, out_gene_dic, inspecies_gene, outspecies_gene, seq_type):
    #Gather together all of the necessary HKA test data from different
    #types of sequence.
    if seq_type == "flank":
        inpoly = in_gene_dic[inspecies_gene].flank_subs
        outpoly = out_gene_dic[outspecies_gene].flank_subs
        #again, trying to make it faster and less memory
#        inseq = "".join(in_gene_dic[inspecies_gene].flank_dic.values())
#        outseq = "".join(out_gene_dic[outspecies_gene].flank_dic.values())
        inseq = in_gene_dic[inspecies_gene].flank_seq
        outseq = out_gene_dic[outspecies_gene].flank_seq
        insample = in_gene_dic[inspecies_gene].average_sample_size(in_gene_dic[inspecies_gene].flank_dic)
        outsample = out_gene_dic[outspecies_gene].average_sample_size(out_gene_dic[outspecies_gene].flank_dic)
    elif seq_type == "intron":
        inpoly = in_gene_dic[inspecies_gene].intron_subs
        outpoly = out_gene_dic[outspecies_gene].intron_subs
        inseq = "".join(in_gene_dic[inspecies_gene].intron_dic.values())
        outseq = "".join(out_gene_dic[outspecies_gene].intro_dic.values())
        insample = in_gene_dic[inspecies_gene].average_sample_size(in_gene_dic[inspecies_gene].intron_dic)
        outsample = out_gene_dic[inspecies_gene].average_sample_size(out_gene_dic[outspecies_gene].intron_dic)
    elif seq_type == "three_prime_utr":
        inpoly = in_gene_dic[inspecies_gene].utr3_subs
        outpoly = out_gene_dic[outspecies_gene].utr3_subs
        inseq = "".join(in_gene_dic[inspecies_gene].utr3_dic.values())
        outseq = "".join(out_gene_dic[outspecies_gene].utr3_dic.values())
        insample = in_gene_dic[inspecies_gene].average_sample_size(in_gene_dic[inspecies_gene].utr3_dic)
        outsample = out_gene_dic[inspecies_gene].average_sample_size(out_gene_dic[outspecies_gene].utr3_dic)
    elif seq_type == "five_prime_utr":
        inpoly = in_gene_dic[inspecies_gene].utr5_subs
        outpoly = out_gene_dic[outspecies_gene].utr5_subs
        inseq = "".join(in_gene_dic[inspecies_gene].utr5_dic.values())
        outseq = "".join(out_gene_dic[outspecies_gene].utr5_dic.values())
        insample = in_gene_dic[inspecies_gene].average_sample_size(in_gene_dic[inspecies_gene].utr5_dic)
        outsample = out_gene_dic[inspecies_gene].average_sample_size(out_gene_dic[outspecies_gene].utr5_dic)
    return inpoly, outpoly, inseq, outseq, insample, outsample
"""
def muscle_pairwise_diff_count(seq1, seq2, inspecies, outspecies, og_num):
    #Align two sequences using muscle and return the number of 
    #differences between them.
#    seq1 = seq1[0:len(seq2)]
#    seq2 = seq2[0:len(seq1)]
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

def mk_test(inspecies, outspecies, ortho_dic, align_dir, out_path, num_threads):
    in_gene_dic = gene_vcf_dic(inspecies, num_threads)
    out_gene_dic = gene_vcf_dic(outspecies, num_threads)
    mk_table = open("%s/%s_%s_mk_table.txt" % (out_path, inspecies, outspecies), 'w')
    mk_table.write("geneID\tPR\tFR\tPS\tFS\tTsil\tTrepl\tnout\tnpop\n")
    for og_num in ortho_dic.keys():
        if inspecies not in ortho_dic[og_num] or outspecies not in ortho_dic[og_num]:
            continue
        if len(ortho_dic[og_num][inspecies]) == 1 and len(ortho_dic[og_num][outspecies]):
#            if og_num != 2110:
#                continue
            reader = SeqIO.parse("%s/og_cds_%s.1.fas" % (align_dir, og_num), format = 'fasta')
            for rec in reader:
                if rec.id[0:4] == inspecies:
                    inseq = str(rec.seq)
                    inseq_name = rec.id[:-3]
                elif rec.id[0:4] == outspecies:
                    outseq = str(rec.seq)
                    outseq_name = rec.id[:-3]
            print og_num
            print inseq_name
            print outseq_name
            fix_syn, fix_nsyn = count_sub_types(inseq, outseq)
            in_poly_syn = in_gene_dic[inseq_name].syn_count
            in_poly_nsyn = in_gene_dic[inseq_name].nsyn_count
            in_potent_syn = in_gene_dic[inseq_name].potent_syn
            in_potent_nsyn = in_gene_dic[inseq_name].potent_nsyn
            in_n_size = in_gene_dic[inseq_name].average_sample_size(in_gene_dic[inseq_name].cds) * 2
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

            
