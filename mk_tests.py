import os
import vcf
from gff import gffParser
import utils

def main():
#    ortho_dic = utils.ortho_reader("/Genomics/kocherlab/berubin/annotation/orthology/proteinortho3.proteinortho")
#    utils.potential_changes_dict()
    for species in ["LMAL"]:
        gene_dic = utils.get_species_data(species)
        reader = vcf.Reader(filename = "/scratch/tmp/berubin/resequencing/%s/genotyping/scaf5_precalls_normal.eff.vcf.gz" % species)

        for gene_name, gene_object in gene_dic.items():
            gene_object.get_genotypes(reader, 3000)

    #    
#        gff_table = gffParser(open("/Genomics/kocherlab/berubin/official_release/LMAL/LMAL_OGS_v1.0_longest_isoform.gff3", 'rU'))
#    print gff_table.geneDict()
    



#reader = vcf.Reader(filename = "/Genomics/kocherlab/berubin/resequencing/selection/LLEU.eff.scaf_0.1Mto1.1M.vcf.recode.vcf.gz")
"""
for rec in reader.fetch("LMAL_scaf_5", 35261, 39382):

#    print rec.ID
#    print rec.REF
#    if rec.ALT[0] == None:
#        continue
    if rec.num_called < 3:
        continue
#    if rec.num_het == 0:
#        continue
#    print rec.QUAL
#    print rec.FILTER
    print rec.samples

#    print rec.genotype
#    print rec.FORMAT
    print rec.POS
    print rec.num_called, rec.num_hom_ref, rec.num_hom_alt, rec.num_het
    print rec.is_monomorphic
    if len(rec.INFO) > 0:
        print rec.INFO
"""
if __name__ == '__main__':
    main()

