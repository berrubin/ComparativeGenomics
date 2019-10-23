import utils
import changes
import collections
from Bio import Seq

class NCAR:
    def __init__(self, gene_name, scaf_name, start, end, strand):
        self.name = gene_name
        self.scaf = scaf_name
        self.cds_coords = []
        self.cds = {}
        self.start = start
        self.end = end
        self.sequence = {}
        self.strand = strand
        self.species = gene_name[0:4]
        self.alts = {}
        self.refs = {}
        self.sample_gts = {}
        self.poly_count = -1
        self.average_n = -1 #npop
        self.called_counts = {}
        self.average_n = 0
        self.neighbors = []
        self.poly_coords = []

    def __str__(self):
        neighbor_names = []
        for neigh in self.neighbors:
            neighbor_names.append(neigh)
        return "Name: %s\nScaffold: %s\nStart coord: %s\nEnd coord: %s\nPoly coords: %s\nPoly count: %s\nAverage #genotypes: %s\nNeighboring genes: %s" % (self.name, self.scaf, self.start, self.end, self.alts.keys(), self.poly_count, self.average_n, ",".join(neighbor_names))

    def potential_poly_sites(self):
        potent_count = 0
        for count in self.called_counts.values():
            if count >= 4:
                potent_count += 1
        return potent_count
 
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
        for k in region_alt_dic.keys():
            old_coord = k
            new_coord = self.new_coding_coord(aligned_seq, self.cds_coords, zeroed_tuples, old_coord)
            if self.strand == 1:
                aligned_vars[new_coord] = region_alt_dic[k]
            else:
                aligned_vars[new_coord] = str(Seq.Seq(str(region_alt_dic[k])).reverse_complement())
        return aligned_vars


    def noncoding_fixed_align(self, aligned_seq):
        #get all of the polymorphisms in the aligned CDS and return
        #a dictionary where the indices have been adjusted to represent
        #the positions in the alignment
        region_alt_dic = {}
        ncar_tuple = (self.start, self.end)
        for k in self.alts.keys():
            if k >= ncar_tuple[0] and k < ncar_tuple[1]:
                region_alt_dic[k] = self.alts[k]
        aligned_vars = {}
        zeroed_tuples = self.zeroed_tuples([ncar_tuple])
        for k in region_alt_dic.keys():
            old_coord = k
            new_coord = self.new_coding_coord(aligned_seq, [ncar_tuple], zeroed_tuples, old_coord)
            if self.strand == 1:
                aligned_vars[new_coord] = region_alt_dic[k]
            else:
                aligned_vars[new_coord] = str(Seq.Seq(str(region_alt_dic[k])).reverse_complement())
        return aligned_vars
                    

    def set_neighbors(self, gene_list):
        #identify the closest neighbors of the gene.
        #it is currently hardcoded to identify 20 neighbors
        #but this can be adjusted to whatever
#        print self.name
        neighbor_dic = {}
        closest_neighbors = []
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
            reader = vcf_reader.fetch(self.scaf, self.start, self.end)
        except ValueError:
            print "can't get that bit of vcf: %s" % self.name
            self.alts = {}
            self.refs = {}
            self.called_counts = {}
#            self.potential_sites()
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
            elif sum(rec.aaf) <= 0.1: #MAF cutoff: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3666677/, http://www.genetics.org/content/158/3/1227.long
                called_count[rec.POS] = rec.num_called
                continue        
            elif sum(rec.aaf)*rec.num_called*2 == 1: #also require allele appears at least twice
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
        #reduce memory usage
        self.refs = {}
        self.poly_count = len(self.alts)
        self.poly_coords = self.alts.keys()
        self.called_counts = called_count
