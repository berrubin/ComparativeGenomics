import utils
import collections
from Bio import Seq

class Gene:
    def __init__(self, gene_name, scaf_name, strand):
        self.name = gene_name
        self.scaf = scaf_name
        self.cds = []
        self.introns = []
        self.five_utrs = []
        self.three_utrs = []
#        self.three_flank = []
#        self.five_flank = []
        self.start = -1
        self.end = -1
        self.sequence = {}
#        self.samples = []
        self.strand = strand
        self.species = gene_name[0:4]
        self.alts = {}
        self.syn_count = -1
        self.nsyn_count = -1
        self.potent_syn
        self.potent_nsyn
        
    def add_cds(self, cds_list):
        self.cds = cds_list
        self.get_introns()

    def add_utrs(self, utr_list, prime_end):
        if prime_end == "three":
            self.three_utrs = utr_list
        elif prime_end == "five":
            self.five_utrs = utr_list

    def set_start(self, start_coord):
        self.start = start_coord
    
    def set_end(self, end_coord):
        self.end = end_coord

    def set_sequence(self, scaf_seq, flank_size):
        seq_dic = collections.OrderedDict()
        if self.start < 0:
            self.start = 0
        if self.end > len(scaf_seq):
            self.end = len(scaf_seq)
        cur_start = self.start - flank_size
        cur_end = self.end + flank_size
        index = cur_start
        while index < cur_end:
            seq_dic[index] = scaf_seq[index-1]
            index += 1
        self.sequence = seq_dic

    def get_introns(self):
        intron_list = []
        for x in range(len(self.cds)-1):
            cur_intron = (self.cds[x][1], self.cds[x+1][0])
            intron_list.append(cur_intron)
        self.introns = intron_list
    
    def get_cds_sequence(self, seq_dic):
        cds_dic = collections.OrderedDict()
        for site, nuc in self.sequence.items():
            for exon in self.cds:
                if site >= exon[0] and site <= exon[1]:
                    cds_dic[site] = nuc
        recon_seq = "".join(cds_dic.values())
        if self.strand == -1:
            recon_seq = str(Seq.Seq(recon_seq).reverse_complement())
        if recon_seq != seq_dic[self.name]:
            raise ImportError("Failed to reconstruct CDS for %s" % self.name)
        self.cds = cds_dic

    def syn_and_nsyn(self):
        syn_count = 0
        nsyn_count = 0
        for index in self.alts.keys():
            old_codon = ""
            new_codon = []
            if index in self.cds.keys():

                if self.strand == 1:
                    cds_index = self.cds.keys().index(index)# - self.cds.keys()[0]
                    codon_start = index - (cds_index % 3)
                    for x in range(3):
                        old_codon += self.cds[codon_start + x]
                        new_codon.append(self.cds[codon_start + x])
                    new_codon[cds_index % 3] = str(self.alts[index])
                    old_codon = Seq.Seq(old_codon)
                    new_codon = Seq.Seq("".join(new_codon))
                elif self.strand == -1:
                    cds_index = len(self.cds.keys()) - 1 - self.cds.keys().index(index)
                    codon_start = index + (cds_index % 3)
                    for x in range(3):
                        old_codon += self.cds[codon_start - x]
                        new_codon.append(self.cds[codon_start - x])
                    new_codon[cds_index % 3] = str(self.alts[index])
                    old_codon = Seq.Seq(old_codon).complement()
                    new_codon = Seq.Seq("".join(new_codon)).complement()
                    print new_codon
                    print old_codon
                if str(new_codon.translate()) == str(old_codon.translate()):
                    syn_count += 1
                elif len(new_codon) % 3 == 0:
                    nsyn_count += 1
                else:
                    #This must be a frameshift mutation which is a disaster and is, therefore, probably not a real variant so don't count it
                    continue
        self.syn_count = syn_count
        self.nsyn_count = nsyn_count


    def get_genotypes(self, vcf_reader, flank_size):
        print self.name
        print self.start
#        print self.sequence
        scaf_len = vcf_reader.contigs[self.scaf][1]
        right_flank_index = self.end
        left_flank_index = self.start
#        sample_dic = {}
#        for sample in vcf_reader.samples:
#            sample_dic[sample] = Sample(sample)
        called_site_count = 0
        alt_dic = {}
        called_count = {}
        for rec in vcf_reader.fetch(self.scaf, left_flank_index, right_flank_index):
            if rec.num_called < 2:
                continue
            if rec.num_called == rec.num_hom_ref:
                called_count[rec.POS] = rec.num_called
                continue
            if rec.num_called == rec.num_het or rec.num_called == rec.num_hom_alt:
                self.sequence[rec.POS] = "N"
                continue
            else:
                alt_dic[rec.POS] = rec.ALT[0]
                called_count[rec.POS] = rec.num_called
        self.alts = alt_dic
        self.syn_and_nsyn()
        self.potential_sites()


                
#                    self.add_site_to_samples(rec, sample_dic, rec.POS)
#                else:
#                    self.add_site_to_samples(rec, sample_dic, rec.POS)
#            print rec.REF
#            print rec.ALT
#            print rec.num_called
#            for s in rec.samples:
#                print s.data.GT

    def potential_sites(self):
        potent_dic = utils.potent_dic()
        cur_cds = "".join(self.cds.values())
        if self.strand == -1:
            cur_cds = str(Seq.Seq(cur_cds).reverse_complement())
        x = 0
        syn_sites = 0
        nsyn_sites = 0
        while x < len(cur_cds):
            codon = cur_cds[x:x+3]
            nsyn_sites += potent_dic["N"][codon]
            syn_sites += potent_dic["S"][codon]
            x = x + 3
        self.potent_syn = syn_sites
        self.potent_nsyn = nsyn_sites
            
    def add_site_to_samples(self, vcf_rec, sample_dic, index):
        if vcf_rec.num_called >= 3:
            if vcf_rec.num_called == vcf_rec.num_het or vcf_rec.num_called == vcf_rec.num_hom_alt:
                for s in vcf_rec.samples:
                    sample_dic[s.sample].add_site(index, "N", "N", s.data.GT)
            else:
                for s in vcf_rec.samples:
                    sample_dic[s.sample].add_site(index, vcf_rec.REF, vcf_rec.ALT, s.data.GT)
            

class Sample:
    def __init__(self, mysample):
        self.sample = mysample
        self.genotypes = {}

    def add_site(self, index, ref_allele, alt_allele, genotype):
        self.genotypes[index] = Genotype(ref_allele, alt_allele, genotype)

class Genotype:
    def __init__(self, ref_allele, alt_allele, genotype):
        self.ref = ref_allele
        self.alt = alt_allele
        self.gt = genotype


