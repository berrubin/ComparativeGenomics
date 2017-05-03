import copy
from Bio import SeqIO
from Region import Gene
from gff import gffParser

def ortho_reader(orthofile):
    #returns dictionary of dictionaries of orthologos genes. 
    #Lower level dictionaries have keys of species codes and 
    #values of gene names. Upper level dictionaries have the
    #OG index (line number in orthofile) as keys.
    reader = open(orthofile, 'rU')
    ortho_dic = {}
    for line in reader:
        counter = 0
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
        print gene_name
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
        print gene_name
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

def mk_test(inspecies, outspecies, align_dir):
    ortho_dic = ortho_reader("/Genomics/kocherlab/berubin/annotation/orthology/proteinortho3.proteinortho", 'rU')
    in_gene_dic = get_species_data(inspecies)
    reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_testset.vcf.gz" % (inspecies, inspecies))
    for gene_name, gene_object in in_gene_dic.items():
        gene_object.get_genotypes(reader, 3000)
    
    out_gene_dic = get_species_data(outspecies)
    reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/%s_testset.vcf.gz" % (outspecies, outspecies))
    for gene_name, gene_object in out_gene_dic.items():
        gene_object.get_genotypes(reader, 3000)

    mk_table = open("%s_%s_mk_table.txt", 'w')
    mk_table.write("geneID\tPR\tFR\tPS\tFS\tTsil\tTrepl\tnout\tnpop\n")
    for og_num in ortho_dic.keys():
        if ortho_dic[og_num].keys().count(inspecies) == 1 and ortho_dic[og_num].keys().count(outspecies) == 1:
            reader = SeqIO.parse("%s/og_cds_%s.1.fas" % (align_dir, og_num))
            for rec in reader:
                if rec.id[0:4] == inspecies:
                    inseq = str(rec.seq)
                    inseq_name = rec.id
                elif rec.id[0:4] == outspecies:
                    outseq = str(rec.seq)
                    outseq_name = rec.id
            fix_syn, fix_nsyn = count_sub_types(inseq, outseq)
            in_poly_syn = in_gene_dic[inseq_name].syn_count
            in_poly_nsyn = in_gene_dic[inseq_name].nsyn_count
            in_potent_syn = in_gene_dic[inseq_name].potent_syn
            in_potent_nsyn = in_gene_dic[inseq_name].potent_nsyn
            in_n_size = in_gene_dic[inseq_name].average_n
            mk_table.write("\t".join([og_num, in_poly_nsyn, fix_nsyn, in_poly_syn, fix_syn, in_potent_syn, in_potent_nsyn, "1", in_n_size]) + "\n")
    mk_table.close()
            

def potential_changes_dict():
    #just ran once to create dictionary of the potential syn and nsyn 
    #sites for each codon. That dictionary is returned by potent_dic() below
    #https://github.com/a1ultima/hpcleap_dnds/blob/master/py/scripts/changes.py

    """ Generate a dictionary, with S and N pre-calculated for all 
    possible pairs of codons (key: pair of codons, value: (S,N).
    ARGS:
        nt_to_aa, a dict mapping codons (keys), e.g. 'TTA', to 
            amino-acid letter (values), e.g. 'L'
            
            e.g. geneticCode("standard")
    Notes:
        Sources of formulae:
        http://www.megasoftware.net/mega4/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_nei_gojobori_method.htm
    """
    nt_to_aa = {  'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', 'AGT':'S','ATA':'I','ATC':'I','ATG':'M','ATT':'I','CAA':'Q','CAC':'H','CAG':'Q','CAT':'H','CCA':'P','CCC':'P','CCG':'P', 'CCT':'P','CGA':'R','CGC':'R','CGG':'R','CGT':'R','CTA':'L','CTC':'L','CTG':'L','CTT':'L','GAA':'E','GAC':'D','GAG':'E', 'GAT':'D','GCA':'A','GCC':'A','GCG':'A','GCT':'A','GGA':'G','GGC':'G','GGG':'G','GGT':'G','GTA':'V','GTC':'V','GTG':'V', 'GTT':'V','TAA':'*','TAC':'Y','TAG':'*','TAT':'Y','TCA':'S','TCC':'S','TCG':'S','TCT':'S','TGA':'*','TGC':'C','TGG':'W', 'TGT':'C','TTA':'L','TTC':'F','TTG':'L','TTT':'F'  }

    potential_changes = {   'S': {  'AAA':0.0,'AAC':0.0,'AAG':0.0,'AAT':0.0, 'ACA':0.0, 'ACC':0.0, 'ACG':0.0, 'ACT':0.0, 'AGA':0.0, 'AGC':0.0, 'AGG':0.0, 'AGT':0.0, 'ATA':0.0, 'ATC':0.0, 'ATG':0.0, 'ATT':0.0, 'CAA':0.0, 'CAC':0.0, 'CAG':0.0, 'CAT':0.0, 'CCA':0.0,'CCC':0.0,'CCG':0.0,'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0, 'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0,'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0, 'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0,'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0, 'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0,'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0},

                            'N': {  'AAA':0.0, 'AAC':0.0, 'AAG':0.0, 'AAT':0.0, 'ACA':0.0, 'ACC':0.0, 'ACG':0.0, 'ACT':0.0, 'AGA':0.0, 'AGC':0.0, 'AGG':0.0, 'AGT':0.0,'ATA':0.0,'ATC':0.0,'ATG':0.0,'ATT':0.0,'CAA':0.0,'CAC':0.0,'CAG':0.0,'CAT':0.0,'CCA':0.0,'CCC':0.0,'CCG':0.0, 'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0,'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0, 'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0,'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0, 'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0,'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0, 'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0}}


    # Mutate (substitutions) all possible codons in the given genetic code, and count proportions of mutations that are synonymous and non-synonmyous
    for codon in nt_to_aa.keys():

        # assert (codon not in codon_record)  @DONE: no duplicate entries
        # codon_record.append(codon)  

        # Calculate S and N (i.e. potential synonymous and potential
        # non-synonymous sites) ()

        # i.e. Of possible mutations that can occur on the codon, 
        # what proportion are synonymous (no aa change)?

        # for each bp position in the codon...
        for codon_p in range(0,2+1):

            nts = ['A','G','T','C']  # @DONE: refactor away, A: we can't, since the next line

            nts.remove(codon[codon_p]) # we do not consider self substitutions, e.g. A->A

            # ...and for each nucleotide that the bp can change 
            # into... 
            for nt in nts:

                codon_mutated = list(copy.deepcopy(codon))
                #codon_mutated = codon
                codon_mutated[codon_p] = nt  # mutate the basepair
                codon_mutated = ''.join(codon_mutated)
                
                # ...count how many of them are synonymous.
                if nt_to_aa[codon]==nt_to_aa[codon_mutated]:
                    potential_changes['S'][codon]+=1.0/3.0 #@DONE: Q: but why 1/3? to me it should be 1/(3*4), A: since it represents 1/3 of a "site"
                else:
                    potential_changes['N'][codon]+=1.0/3.0 #@DONE: Q: but why 1/3? to me it should be 1/(3*4), A: since it represents 1/3 of a "site"
    outfile = open("potential_changes_dic.txt", 'w')
    outfile.write(str(potential_changes))
    outfile.close()

def potent_dic():
    return {'S': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.3333333333333333, 'AAA': 0.3333333333333333, 'ATC': 0.6666666666666666, 'AAC': 0.3333333333333333, 'ATA': 0.6666666666666666, 'AGG': 0.6666666666666666, 'CCT': 1.0, 'CTC': 1.0, 'AGC': 0.3333333333333333, 'ACA': 1.0, 'AGA': 0.6666666666666666, 'CAT': 0.3333333333333333, 'AAT': 0.3333333333333333, 'ATT': 0.6666666666666666, 'CTG': 1.3333333333333333, 'CTA': 1.3333333333333333, 'ACT': 1.0, 'CAC': 0.3333333333333333, 'ACG': 1.0, 'CAA': 0.3333333333333333, 'AGT': 0.3333333333333333, 'CAG': 0.3333333333333333, 'CCG': 1.0, 'CCC': 1.0, 'TAT': 0.3333333333333333, 'GGT': 1.0, 'TGT': 0.3333333333333333, 'CGA': 1.3333333333333333, 'CCA': 1.0, 'CGC': 1.0, 'GAT': 0.3333333333333333, 'CGG': 1.3333333333333333, 'CTT': 1.0, 'TGC': 0.3333333333333333, 'GGG': 1.0, 'TAG': 0.3333333333333333, 'GGA': 1.0, 'TAA': 0.6666666666666666, 'GGC': 1.0, 'TAC': 0.3333333333333333, 'GAG': 0.3333333333333333, 'TCG': 1.0, 'TTA': 0.6666666666666666, 'TTT': 0.3333333333333333, 'GAC': 0.3333333333333333, 'CGT': 1.0, 'GAA': 0.3333333333333333, 'TCA': 1.0, 'GCA': 1.0, 'GTA': 1.0, 'GCC': 1.0, 'GTC': 1.0, 'GCG': 1.0, 'GTG': 1.0, 'TTC': 0.3333333333333333, 'GTT': 1.0, 'GCT': 1.0, 'TGA': 0.3333333333333333, 'TTG': 0.6666666666666666, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 1.0}, 'N': {'ACC': 1.9999999999999998, 'ATG': 3.0, 'AAG': 2.6666666666666665, 'AAA': 2.6666666666666665, 'ATC': 2.333333333333333, 'AAC': 2.6666666666666665, 'ATA': 2.333333333333333, 'AGG': 2.333333333333333, 'CCT': 1.9999999999999998, 'CTC': 1.9999999999999998, 'AGC': 2.6666666666666665, 'ACA': 1.9999999999999998, 'AGA': 2.333333333333333, 'CAT': 2.6666666666666665, 'AAT': 2.6666666666666665, 'ATT': 2.333333333333333, 'CTG': 1.6666666666666665, 'CTA': 1.6666666666666665, 'ACT': 1.9999999999999998, 'CAC': 2.6666666666666665, 'ACG': 1.9999999999999998, 'CAA': 2.6666666666666665, 'AGT': 2.6666666666666665, 'CAG': 2.6666666666666665, 'CCG': 1.9999999999999998, 'CCC': 1.9999999999999998, 'TAT': 2.6666666666666665, 'GGT': 1.9999999999999998, 'TGT': 2.6666666666666665, 'CGA': 1.6666666666666665, 'CCA': 1.9999999999999998, 'CGC': 1.9999999999999998, 'GAT': 2.6666666666666665, 'CGG': 1.6666666666666665, 'CTT': 1.9999999999999998, 'TGC': 2.6666666666666665, 'GGG': 1.9999999999999998, 'TAG': 2.6666666666666665, 'GGA': 1.9999999999999998, 'TAA': 2.333333333333333, 'GGC': 1.9999999999999998, 'TAC': 2.6666666666666665, 'GAG': 2.6666666666666665, 'TCG': 1.9999999999999998, 'TTA': 2.333333333333333, 'TTT': 2.6666666666666665, 'GAC': 2.6666666666666665, 'CGT': 1.9999999999999998, 'GAA': 2.6666666666666665, 'TCA': 1.9999999999999998, 'GCA': 1.9999999999999998, 'GTA': 1.9999999999999998, 'GCC': 1.9999999999999998, 'GTC': 1.9999999999999998, 'GCG': 1.9999999999999998, 'GTG': 1.9999999999999998, 'TTC': 2.6666666666666665, 'GTT': 1.9999999999999998, 'GCT': 1.9999999999999998, 'TGA': 2.6666666666666665, 'TTG': 2.333333333333333, 'TCC': 1.9999999999999998, 'TGG': 3.0, 'TCT': 1.9999999999999998}}
