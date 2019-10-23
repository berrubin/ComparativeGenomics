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

def fourfold_codons():
    return ["GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "CGG", "GGT", "GGC", "GGA", "GGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "TCT", "TCC", "TCA", "TCG", "ACT", "ACC", "ACA", "ACG", "GTT", "GTC", "GTA", "GTG"]

def aa_types(): #https://en.wikipedia.org/wiki/Conservative_replacement
    return {'G': "aliphatic", 'A': "aliphatic", 'V': "aliphatic", 'L': "aliphatic", 'I': "aliphatic", 'S': "oxyl", 'C': "oxyl", 'U': "oxyl", 'T': "oxyl", 'M': "oxyl", 'P': "cyclic", 'F': "aromatic", 'Y': "aromatic", 'W': "aromatic", 'H': "basic", 'K': "basic", 'R': "basic", 'D': "acidic", 'E': "acidic", 'N': "acidic", 'Q': "acidic", "*": "stop"}
