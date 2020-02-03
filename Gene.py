import numpy
import utils
import changes
import collections
from Bio import Seq

class Gene:
    def __init__(self, gene_name, scaf_name, strand):
        self.name = gene_name
        self.scaf = scaf_name
        self.cds_coords = []
        self.cds = {}
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
        self.sample_gts = {}
        self.fourfold = -1
        self.syn_count = -1 #PS
        self.nsyn_count = -1 #PR
        self.potent_syn = -1 #Tsil
        self.potent_nsyn = -1 #Trepl
        self.potent_fourfold = -1
        self.average_n = -1 #npop
        self.cds_average_n = 0
        self.intron_dic = {}
        self.first_intron_dic = {}
        self.utr3_dic = {}
        self.utr5_dic = {}
        self.flank_dic = {}
        self.flank_subs = -1 
        self.intron_subs = -1
        self.first_intron_subs = -1
        self.utr3_subs = -1
        self.utr5_subs =-1
        self.called_counts = {}
        self.flank_start = -1
        self.flank_end = -1
        self.average_n = 0
        self.flank_seq = ""
        self.intron_seq = ""
        self.first_intron_seq = ""
        self.neighbors = []
        self.syn_coords = []
        self.nsyn_coords = []
        self.fourf_coords = []
        self.flank_sub_coords = []
        self.sample_nsyn_counts = {}
        self.sample_syn_counts = {}
        self.sample_nsyn_genod = {}
        self.sample_syn_genod = {}
        self.low_impact_count = -1
        self.hi_impact_count = -1
        self.low_impact_coords = []
        self.hi_impact_coords = []
        self.sample_low_counts = {}
        self.sample_hi_counts = {}
        self.sample_low_genod = {}
        self.sample_hi_genod = {}
 
    def get_new_coord(self, aligned_seq, old_coord):
        new_coord = old_coord
        old_index = 0
        new_index = old_coord
        total_gap = 0
        if aligned_seq[:new_index].count("-") == 0:
            return new_index
        while True:
            gap_count = aligned_seq[old_index:new_index].count("-")
            total_gap = total_gap + gap_count
            old_index = new_index
            if gap_count == 0:
                new_index = new_index + 1
            else:
                new_index = new_index + gap_count
            
            if len(aligned_seq[:new_index]) - aligned_seq[:new_index].count("-") == old_coord:
                return new_index

    def aligned_alts(self, aligned_seq, start_coord):
        region_alt_dic = {}
        for k in self.alts.keys():
            if k >= start_coord and k < start_coord + len(aligned_seq) - aligned_seq.count("-"):
                region_alt_dic[k-start_coord] = self.alts[k]
        aligned_vars = {}
        for k in region_alt_dic.keys():
            new_coord = self.get_new_coord(aligned_seq, k)
            if self.strand == 1:
                aligned_vars[new_coord] = region_alt_dic[k]
            else:
                aligned_vars[len(aligned_seq) - new_coord] = str(Seq.Seq(str(region_alt_dic[k])).reverse_complement())
        return aligned_vars
                

    def check_polymorphisms_fixed(self, align_dic, start_coord, in_or_out):
        #this method readjusts the coordinates for polymorphisms into
        #an aligned sequence space
        #"start_coord" is just the beginning of whatever feature is 
        #currently being examined. Use this for noncoding sequences
        cur_seq = align_dic[in_or_out]
        if self.strand == -1:
            cur_seq = str(Seq.Seq(cur_seq).reverse_complement())
        aligned_vars = self.aligned_alts(cur_seq, start_coord)
        return aligned_vars

    def conserved_region_polys(self, seq_type, new_coords):
        #this method is for counting the number of polymorphisms in
        #flanking sequence after it is limited down to the regions
        #that can actually be aligned
        target_start = 0
        if seq_type == "flank":
            if self.strand == 1:
                target_start = new_coords[0] + self.flank_start
                target_end = self.flank_start + (new_coords[1] - new_coords[0])
            else:
                target_end = self.flank_end - new_coords[0]
                target_start = self.flank_end - (new_coords[1] - new_coords[0])
                
        cur_seq = align_dic[in_or_out]
        if self.strand == -1:
            cur_seq = str(Seq.Seq(cur_seq).reverse_complement())
        aligned_vars = self.aligned_alts(cur_seq, start_coord)
        return aligned_vars

            
    def zeroed_tuples(self, tuple_list):
        #take a list of tuples (CDS starts and ends or intron starts
        #and ends) and map them to new tuples that represent their
        #positions in a sequence space without intervening sequence
        #(e.g. aligned CDS)
        new_tuple_dic = {}
        base_index = tuple_list[0]
        new_tuple = (0, base_index[1] - base_index[0])
        new_tuple_dic[base_index] = new_tuple
        intron_length = 0
        prev_coord = tuple_list[0]
        for coord_tuple in tuple_list:
            if coord_tuple == base_index:
                continue
            intron_length += coord_tuple[0] - prev_coord[1] -1
            adjustment = intron_length + base_index[0]
            new_tuple = (coord_tuple[0] - adjustment, coord_tuple[1] - adjustment)
            new_tuple_dic[coord_tuple] = new_tuple
            base_index = (tuple_list[0][0], coord_tuple[1])
            prev_coord = coord_tuple
        return new_tuple_dic

    def new_coding_coord(self, aligned_seq, coord_list, zeroed_tuples, old_coord):
        #return an individual coordinate zeroed to the new coordinate
        #space
        if self.strand == -1:
            cur_seq = str(Seq.Seq(aligned_seq).reverse_complement())
        else:
            cur_seq = aligned_seq
        for coord_tuple in coord_list:
            if old_coord >= coord_tuple[0] and old_coord < coord_tuple[1]:
                target_coords = coord_tuple
        new_tuple_dic = zeroed_tuples
        shift = target_coords[0] - new_tuple_dic[target_coords][0] 
        new_coord = old_coord - shift
        if self.strand == 1:
            return new_coord
        else:
            return len(aligned_seq) - new_coord - 1

    def coding_fixed_align(self, aligned_seq):
        #get all of the polymorphisms in the aligned CDS and return
        #a dictionary where the indices have been adjusted to represent
        #the positions in the alignment
        region_alt_dic = {}
        for k in self.alts.keys():
            for cds_tuple in self.cds_coords:
                if k >= cds_tuple[0] and k < cds_tuple[1]:
                    region_alt_dic[k] = self.alts[k]
        aligned_vars = {}
        zeroed_tuples = self.zeroed_tuples(self.cds_coords)
#        print region_alt_dic
        for k in region_alt_dic.keys():
            old_coord = k
            new_coord = self.new_coding_coord(aligned_seq, self.cds_coords, zeroed_tuples, old_coord)
            if self.strand == 1:
                aligned_vars[new_coord] = region_alt_dic[k]
            else:
                aligned_vars[new_coord] = str(Seq.Seq(str(region_alt_dic[k])).reverse_complement())
#            print "%s\t%s" % (old_coord, new_coord)
        return aligned_vars
                    

    def set_neighbors(self, gene_list):
        #identify the closest neighbors of the gene.
        #it is currently hardcoded to identify 20 neighbors
        #but this can be adjusted to whatever
        print self.name
        neighbor_dic = {}
        closest_neighbors = [self.name]
        for gene in gene_list:
            if gene.scaf != self.scaf:
                continue
            if gene.start > self.start and gene.end > self.end:
                neighbor_dic[gene.start - self.end] = gene.name
            elif gene.start < self.start and gene.end < self.end:
                neighbor_dic[self.start - gene.end] = gene.name
        for distance in sorted(neighbor_dic):
            if len(closest_neighbors) > 20:
                break
            if len(closest_neighbors) <= 20:
                closest_neighbors.append(neighbor_dic[distance])
        self.neighbors = closest_neighbors
        neighbor_dic = {}

    def add_cds(self, cds_list):
        self.cds_coords = cds_list
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
        cur_start = self.start - flank_size
        cur_end = self.end + flank_size
        if cur_start < 0:
            cur_start = 0
        if cur_end > len(scaf_seq):
            cur_end = len(scaf_seq)
        index = cur_start
        while index <= cur_end:
            seq_dic[index] = scaf_seq[index-1]
            index += 1
        self.sequence = seq_dic
        self.flank_start = cur_start
        self.flank_end = cur_end

    def get_introns(self):
        intron_list = []
        for x in range(len(self.cds_coords)-1):
            cur_intron = (self.cds_coords[x][1], self.cds_coords[x+1][0])
            intron_list.append(cur_intron)
        self.introns = intron_list

    def get_intron_sequence(self):
        intron_dic = collections.OrderedDict()
        for site, nuc in self.sequence.items():
            for intron in self.introns:
                if site >= intron[0] and site <= intron[1]:
                    intron_dic[site] = nuc
        self.intron_dic = intron_dic
        self.intron_seq = "".join(intron_dic.values())

    def get_first_intron_sequence(self):
        intron_dic = collections.OrderedDict()
        for site, nuc in self.sequence.items():
            if len(self.introns) > 0:
                if self.strand == 1:
                    intron = self.introns[0]
                else:
                    intron = self.introns[-1]
                if site >= intron[0] and site <= intron[1]:
                    intron_dic[site] = nuc
        self.first_intron_dic = intron_dic
        self.first_intron_seq = "".join(intron_dic.values())


    def get_utr_sequence(self):
        utr_dic = collections.OrderedDict()
        for site, nuc in self.sequence.items():
            for utr in self.three_utrs:
                if site >= utr[0] and site <= utr[1]:
                    utr_dic[site] = nuc
        self.utr3_dic = utr_dic
        utr_dic = collections.OrderedDict()
        for site, nuc in self.sequence.items():
            for utr in self.five_utrs:
                if site >= utr[0] and site <= utr[1]:
                    utr_dic[site] = nuc
        self.utr5_dic = utr_dic

    def get_flank_sequence(self, flank_size):
        flank_dic = collections.OrderedDict()
        for site, nuc in self.sequence.items():
            if self.strand == 1:
                if site < self.start:
                    flank_dic[site] = nuc
            elif self.strand == -1:
                if site > self.end:
                    flank_dic[site] = nuc
        self.flank_dic = flank_dic
        self.flank_seq = "".join(flank_dic.values())
    
    def get_cds_sequence(self):#, seq_dic):
        cds_dic = collections.OrderedDict()
        for site, nuc in self.sequence.items():
            for exon in self.cds_coords:
                if site >= exon[0] and site <= exon[1]:
                    cds_dic[site] = nuc
        #if you want to check the reconstructed sequences against the 
        #sequences in "seq_dic" (i.e. those in a CDS file somwhere), 
        #that can be done with the code below.
#        recon_seq = "".join(cds_dic.values())
#        if self.strand == -1:
#            recon_seq = str(Seq.Seq(recon_seq).reverse_complement())
#        if recon_seq != seq_dic[self.name]:
#            raise ImportError("Failed to reconstruct CDS for %s" % self.name)
        self.cds = cds_dic 
        cds_dic = {}

    def noncoding_subs(self):
        flank_subs = 0
        intron_subs = 0
        first_intron_subs = 0
        utr3_subs = 0
        utr5_subs = 0
        for index in self.alts.keys():
            if index in self.flank_dic.keys():
                if len(self.alts[index]) != len(self.refs[index]):
                    continue
                for x in range(len(self.alts[index])):
                    if x > len(self.refs[index]) - 1:
                        break
                    if str(self.alts[index])[x] != str(self.refs[index])[x]:
                        self.flank_sub_coords.append(index)
                        flank_subs += 1
            if index in self.intron_dic.keys():
                for x in range(len(self.alts[index])):
                    if x > len(self.refs[index]) - 1:
                        break
                    if str(self.alts[index])[x] != str(self.refs[index])[x]:
                        intron_subs += 1
            if index in self.first_intron_dic.keys():
                for x in range(len(self.alts[index])):
                    if x > len(self.refs[index]) - 1:
                        break
                    if str(self.alts[index])[x] != str(self.refs[index])[x]:
                        first_intron_subs += 1            
            if index in self.utr3_dic.keys():                
                for x in range(len(self.alts[index])):
                    if x > len(self.refs[index]) - 1:
                        break
                    if str(self.alts[index])[x] != str(self.refs[index])[x]:
                        utr3_subs += 1
            if index in self.utr5_dic.keys():                
                for x in range(len(self.alts[index])):
                    if x > len(self.refs[index]) - 1:
                        break
                    if str(self.alts[index])[x] != str(self.refs[index])[x]:
                        utr5_subs += 1
        self.flank_subs = flank_subs
        self.intron_subs = intron_subs
        self.first_intron_subs = first_intron_subs
        self.utr3_subs = utr3_subs
        self.utr5_subs = utr5_subs

    def syn_and_nsyn(self):
        #count the number of synonymous and nonsynonymous polymorphisms
        #in this gene
        syn_count = 0
        nsyn_count = 0
        fourfold_count = 0
        low_impact_count = 0
        hi_impact_count = 0
        nsyn_coords = []
        syn_coords = []
        low_impact_coords = []
        hi_impact_coords = []
        fourf_coords = []
        aa_types_dic = changes.aa_types()
        for index in self.alts.keys():
            old_codon = ""
            new_codon = []
            if index in self.cds.keys():
                if len(self.alts[index]) != len(self.refs[index]):
                    nsyn_count += 1
                    continue
                variant_len = len(self.alts[index])
                #########TODO: Do unit tests of all lengths of alts genotypes at all locations within codons (1-3). Many of these have been done but should do systematically.
#                print self.alts[index]
                if self.strand == 1:
                    cds_index = self.cds.keys().index(index)

                    if cds_index + variant_len >= len(self.cds):
                        variant_len = len(self.cds) - cds_index - 1 
                        self.alts[index] = str(self.alts[index])[0:variant_len]

                    codon_start = cds_index - (cds_index % 3)
                    leftover_len = variant_len - (3 - cds_index % 3)
                    additional_len = 3 * (leftover_len / 3)
                    if leftover_len % 3 > 0:
                        additional_len += 3
                    for x in range(3+additional_len):
                        new_index = codon_start + x
                        new_index = self.cds.keys()[new_index]
                        old_codon += self.cds[new_index]
                        new_codon.append(self.cds[new_index])
                    for x in range(variant_len):
                        new_codon[cds_index % 3 + x] = str(self.alts[index])[x]
                    old_codon = Seq.Seq(old_codon)
                    new_codon = Seq.Seq("".join(new_codon))

                elif self.strand == -1:
#                    print self.cds

                    cds_index = len(self.cds.keys()) - 1 - self.cds.keys().index(index) + 1 - (variant_len % 3)

#                    print index
#                    print len(self.cds.keys())
#                    print self.refs[index]
#                    print variant_len
#                    print index
#                    print len(self.cds.keys())
#                    print self.cds.keys().index(index)
                    if cds_index - variant_len < 0:
                        variant_len = cds_index
                        self.alts[index] = str(self.alts[index])[0:variant_len]
                    codon_start = len(self.cds.keys()) - 1 - cds_index + (cds_index % 3) + ((variant_len / 3) * 3)
#                    codon_start = 
                    leftover_len = variant_len - (3 - cds_index % 3)
#                    print "leftover: %s" % leftover_len
                    additional_len = 3 * (leftover_len / 3)
                    if leftover_len % 3 > 0:
                        additional_len += 3

                    for x in range(3 + additional_len):
                        new_index = codon_start - x
                            
                        new_index = self.cds.keys()[new_index]
                        old_codon += self.cds[new_index]
                        new_codon.append(self.cds[new_index])
                    for x in range(variant_len):
#                        print new_codon
#                        print cds_index
#                        print codon_start
#                        print cds_index % 3
#                        print x
#                        print str(self.alts[index])[variant_len + x - 1]
                        new_codon[cds_index % 3 + x] = str(self.alts[index])[variant_len - x - 1]
####the following line is trash. Was an attempted fix for an indexing problem but it did not work across situations
#                        new_codon[(cds_index - variant_len + 1) % 3 + x] = str(self.alts[index])[variant_len - x - 1]
#                        print new_codon
                    old_codon = Seq.Seq(old_codon).complement()
                    new_codon = Seq.Seq("".join(new_codon)).complement()
                old_aa = str(old_codon.translate())
                new_aa = str(new_codon.translate())
                if self.check_fourfold(old_codon, new_codon):
                    fourfold_count += 1
                    fourf_coords.append(index)
                temp_nsyn = 0
#                print codon_start
#                print cds_index
#                print old_aa
#                print str(old_codon)
#                print new_aa
#                print str(new_codon)
#                print variant_len
                for aa in range(len(new_aa)):
                    if old_aa[aa] != new_aa[aa]:
                        nsyn_count += 1
                        nsyn_coords.append(index)
                        temp_nsyn += 1
                        if self.is_low_impact(old_aa, new_aa):
                            low_impact_count += 1
                            low_impact_coords.append(index)
                        else:
                            hi_impact_count += 1
                            hi_impact_coords.append(index)
                            
#                        print index
#                        print old_aa[aa]
#                        print new_aa[aa]

                temp_syn = 0
                for nuc in range(len(str(new_codon))):
                    if str(old_codon)[nuc] != str(new_codon)[nuc]:
                        syn_count += 1
                        if index not in nsyn_coords:
                            syn_coords.append(index)
                        temp_syn 
                syn_count = syn_count - temp_nsyn
        self.syn_count = syn_count
        self.nsyn_count = nsyn_count
        self.fourfold = fourfold_count
#        print nsyn_coords
#        print syn_coords
#        print len(nsyn_coords)
#        print len(syn_coords)
#        print nsyn_count
#        print syn_count
        self.syn_coords = syn_coords
        self.nsyn_coords = nsyn_coords
        self.low_impact_count = low_impact_count
        self.hi_impact_count = hi_impact_count
        self.low_impact_coords = low_impact_coords
        self.hi_impact_coords = hi_impact_coords
#        self.fourf_coords = fourf_coords

    def is_low_impact(self, old_aa, new_aa):
        aa_types_dic = changes.aa_types()
        low_impact = True
        for position in range(len(old_aa)):
            if aa_types_dic[old_aa[position]] != aa_types_dic[new_aa[position]]:
                low_impact = False
        return low_impact

    def check_fourfold(self, old_codon, new_codon):
        #check if two codons are fourfold degenerate sites of the same
        #amino acid
        if str(old_codon) in changes.fourfold_codons():
            if str(old_codon)[0] == str(new_codon)[0] and str(old_codon)[1] == str(old_codon)[1] and str(old_codon)[2] != str(new_codon)[2]:
                return True
        return False

    def average_sample_size(self, target_dic):
        #returns average number of samples genotyped at this gene
        #the complexity of this was reduced to decrease memory use
        #and runtime but if exact numbers are required, that is 
        #possible to obtain by use of self.called_counts
        return self.average_n
        size_sum = 0
        num_sites = 0
        for site in target_dic.keys():
            num_sites += 1
            size_sum += self.called_counts[site]
        if num_sites == 0:
            return 0
        return size_sum / num_sites

    def get_genotypes(self, vcf_reader, flank_size):
        #read in polymorphism data from VCF file
        called_site_count = 0
        alt_dic = {}
        called_count = {}
        ref_dic = {}
        try:
            reader = vcf_reader.fetch(self.scaf, self.flank_start, self.flank_end)
        except ValueError:
            print "can't get that bit of vcf"
            self.alts = {}
            self.refs = {}
            self.called_counts = {}
            self.potential_sites()
            return
        scaf_len = vcf_reader.contigs[self.scaf][1]
        gene_sample_dic = {}
        for rec in reader: 
            if rec.num_called < 4:
                continue
            elif rec.num_called == rec.num_hom_ref:
                called_count[rec.POS] = rec.num_called
                continue
            elif rec.num_called == rec.num_het or rec.num_called == rec.num_hom_alt:
                self.sequence[rec.POS] = "N"
                continue
            elif sum(rec.aaf) <= 0.1:
                called_count[rec.POS] = rec.num_called
                continue
            elif sum(rec.aaf)*rec.num_called*2 == 1:
                called_count[rec.POS] = rec.num_called
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
                gene_sample_dic[rec.POS] = {}
                for s in rec.samples:
                    gene_sample_dic[rec.POS][s.sample] = s['GT']

        #if you need exact counts of sample numbers, use this
#        self.called_counts = called_count
        if len(called_count.values()) == 0:
            self.average_n = 0
#        elif "LALB" in self.name:
#            self.average_n = round(1.0 * sum(called_count.values()) / len(called_count.values()))
#        else:
        self.average_n = round(1.0 * sum(called_count.values()) / len(self.sequence))
        
        self.alts = alt_dic
        self.refs = ref_dic
        self.sample_gts = gene_sample_dic
        self.syn_and_nsyn()
        cds_sample = []
        for site in self.cds.keys():
            if called_count.get(site, None) != None:
                cds_sample.append(called_count[site])
            else:
                cds_sample.append(0)
        self.cds_average_n = 1.0 * sum(cds_sample) / len(cds_sample)
        self.syn_nsyn_by_sample()
        self.potential_sites()
        self.noncoding_subs()
        #reduce memory usage
        self.refs = {}
#        self.alts = {}
        self.cds = {}
        self.flank_dic = {}
        self.intron_dic = {}
        self.sequence = {}

    def syn_nsyn_by_sample(self):
        #These are the numbers of polymorphic sites where individuals are heterozygous
        sample_nsyn_counts = {}
        sample_syn_counts = {}
        sample_low_counts = {}
        sample_hi_counts = {}
        #Below are the numbers of polymorphic sites where an individual has been genotyped
        sample_nsyn_genod = {}
        sample_syn_genod = {}
        sample_low_genod = {}
        sample_hi_genod = {}
        for position, gt_dic in self.sample_gts.items():
            if position in self.nsyn_coords:
                for sample, gt in gt_dic.items():
                    if sample_nsyn_counts.get(sample, None) == None:
                        sample_nsyn_counts[sample] = 0
                        sample_nsyn_genod[sample] = 0
                        sample_hi_counts[sample] = 0
                        sample_hi_genod[sample] = 0
                        sample_low_counts[sample] = 0
                        sample_low_genod[sample] = 0
                    if gt in ["0/1", "1/0"]:
                        sample_nsyn_counts[sample] += 1
                        if position in self.low_impact_coords:
                            sample_low_counts[sample] += 1
                        elif position in self.hi_impact_coords:
                            sample_hi_counts[sample] += 1
                    if gt in ["0/1", "1/0", "1/1", "0/0"]:
                        sample_nsyn_genod[sample] += 1
                        if position in self.low_impact_coords:
                            sample_low_genod[sample] += 1
                        elif position in self.hi_impact_coords:
                            sample_hi_genod[sample] += 1
            elif position in self.syn_coords:
                for sample, gt in gt_dic.items():
                    if sample_syn_counts.get(sample, None) == None:
                        sample_syn_counts[sample] = 0
                        sample_syn_genod[sample] = 0
                    if gt in ["0/1", "1/0"]:
                        sample_syn_counts[sample] += 1
                    if gt in ["0/1", "1/0", "1/1", "0/0"]:
                        sample_syn_genod[sample] += 1
        self.sample_nsyn_counts = sample_nsyn_counts
        self.sample_syn_counts = sample_syn_counts
        self.sample_nsyn_genod = sample_nsyn_genod
        self.sample_syn_genod = sample_syn_genod
        self.sample_low_counts = sample_low_counts
        self.sample_hi_counts = sample_hi_counts
        self.sample_low_genod = sample_low_genod
        self.sample_hi_genod = sample_hi_genod
#        print sample_nsyn_counts
#        print sample_syn_counts

    def potential_sites(self):
        potent_dic = changes.potent_dic()
        fourfold_list = changes.fourfold_codons()
        cur_cds = "".join(self.cds.values())
        if self.strand == -1:
            cur_cds = str(Seq.Seq(cur_cds).reverse_complement())
        x = 0
        syn_sites = 0
        nsyn_sites = 0
        potent_fourfold = 0
        while x < len(cur_cds):
            ambig_codon = False
            codon = cur_cds[x:x+3].upper()
            for ambig in ["N", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B"]:
                if ambig in codon:
                    ambig_codon = True
            if ambig_codon:
                x = x + 3
                ambig_codon = False
                continue
            nsyn_sites += potent_dic["N"][codon]
            syn_sites += potent_dic["S"][codon]
            if codon in fourfold_list:
                potent_fourfold += 1
            x = x + 3
        self.potent_syn = syn_sites
        self.potent_nsyn = nsyn_sites
        self.potent_fourfold = potent_fourfold
            
    def add_site_to_samples(self, vcf_rec, sample_dic, index):
        #I don't use this. this increases complexity, RAM, and runtime.
        #however, identifying individual samples can allow for the
        #construction of individual genotypes which might be useful
        #for something someday.
        if vcf_rec.num_called >= 3:
            if vcf_rec.num_called == vcf_rec.num_het or vcf_rec.num_called == vcf_rec.num_hom_alt:
                for s in vcf_rec.samples:
                    sample_dic[s.sample].add_site(index, "N", "N", s.data.GT)
            else:
                for s in vcf_rec.samples:
                    sample_dic[s.sample].add_site(index, vcf_rec.REF, vcf_rec.ALT, s.data.GT)
            

