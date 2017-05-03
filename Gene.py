import utils
import changes
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
        self.start = -1
        self.end = -1
        self.sequence = {}
        self.strand = strand
        self.species = gene_name[0:4]
        self.alts = {}
        self.refs = {}
        self.syn_count = -1 #PS
        self.nsyn_count = -1 #PR
        self.potent_syn = -1 #Tsil
        self.potent_nsyn = -1 #Trepl
        self.average_n = -1 #npop
        
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
        print self.name
        print self.start
        
        for index in self.alts.keys():
            old_codon = ""
            new_codon = []
            if index in self.cds.keys():
                if len(self.alts[index]) != len(self.refs[index]):
                    nsyn_count += 1
                    continue
                if len(self.alts[index]) < 2:
                    continue
                variant_len = len(self.alts[index])
                if self.strand == 1:
                    cds_index = self.cds.keys().index(index)
                    codon_start = index - (cds_index % 3)
                    
                    for x in range(3 + ((variant_len / 3) * 3)):
                        old_codon += self.cds[codon_start + x]
                        new_codon.append(self.cds[codon_start + x])
                    new_codon[cds_index % 3] = str(self.alts[index])
                    for x in range(variant_len):
                        new_codon[cds_index % 3 + x] = str(self.alts[index])[x]
                    old_codon = Seq.Seq(old_codon)
                    new_codon = Seq.Seq("".join(new_codon))
                elif self.strand == -1:
                    cds_index = len(self.cds.keys()) - 1 - self.cds.keys().index(index)
                    codon_start = index + (cds_index % 3) + ((variant_len / 3) * 3)
                    for x in range(3 + ((variant_len / 3) * 3)):
                        old_codon += self.cds[codon_start - x]
                        new_codon.append(self.cds[codon_start - x])
                    new_codon[cds_index % 3] = str(self.alts[index])
                    for x in range(variant_len):
                        new_codon[cds_index % 3 + x] = str(self.alts[index])[variant_len - x - 1]
                    old_codon = Seq.Seq(old_codon).complement()
                    new_codon = Seq.Seq("".join(new_codon)).complement()
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
#        print self.sequence
        scaf_len = vcf_reader.contigs[self.scaf][1]
        right_flank_index = self.end
        left_flank_index = self.start
        called_site_count = 0
        alt_dic = {}
        called_count = {}
        ref_dic = {}
        for rec in vcf_reader.fetch(self.scaf, left_flank_index, right_flank_index):
            if rec.num_called < 4:
                continue
            elif rec.num_called == rec.num_hom_ref:
                called_count[rec.POS] = rec.num_called
                continue
            elif rec.num_called == rec.num_het or rec.num_called == rec.num_hom_alt:
                self.sequence[rec.POS] = "N"
                continue
            half_missing = 0
            for s in rec.samples:
                if str(s.data.GT) in ["./1", "./0", "1/.", "0/."]:
                    half_missing += 1
            if (rec.num_called - half_missing) == rec.num_hom_alt or (rec.num_called - half_missing) == rec.num_het:
                continue
            else:
                alt_dic[rec.POS] = rec.ALT[0]
                ref_dic[rec.POS] = rec.REF
                called_count[rec.POS] = rec.num_called
        self.average_n = sum(called_count.values()) / len(called_count.values())
        self.alts = alt_dic
        self.refs = ref_dic
        self.syn_and_nsyn()
        self.potential_sites()

    def potential_sites(self):
        potent_dic = changes.potent_dic()
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
            

