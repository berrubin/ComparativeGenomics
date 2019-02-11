####This file was downloaded from ftp://parrot.genomics.cn/gigadb/pub/10.5524/101001_102000/101041/Scripts.tar.gz on 1/31/2019
####used the Jarvis et al. (Science 346: 1320-1331) Avian Phylogenomics Project scripts for amino acid alignments, which masks over poorly aligning regions of individual sequences (rather than omitting entire alignment columns)

from numpy import *
import os
import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
import re
import multiprocessing
from multiprocessing import Pool

def jarvis_filter(og_list, indir, outdir, len_min, num_threads):
    window=12
    step=1
    maxGapContent=0.6
    aaMatrix="/Genomics/kocherlab/berubin/local/developing/selection_pipeline/blosum62.txt"
    pool = multiprocessing.Pool(processes = num_threads)

    matrixAA=readMatrix(aaMatrix)
    work_list = []
    for cur_og in og_list:
        outfile = "%s/og_pep_%s_WrongA.txt" % (outdir, cur_og)
        if os.path.exists(outfile):
            reader = open(outfile, 'rU')
            if "END_OF_MASKS" in reader.readlines():
                continue #don't rerun an already finished mask
        infile = "%s/og_cds_%s.afa" % (indir, cur_og)
        translate_file = open("%s/og_pep_%s.afa" % (outdir, cur_og), 'w')
        reader = SeqIO.parse(infile, format = 'fasta')
        for rec in reader:
            cur_seq = Seq.Seq(str(rec.seq).replace("-", "N")).translate()
            translate_file.write(">%s\n%s\n" % (rec.id, str(cur_seq).replace("*","X").replace("J","X").replace("Z","X").replace("B","X")))
        translate_file.close()
        infile = "%s/og_pep_%s.afa" % (outdir, cur_og)
        work_list.append([infile, outfile, window, step, maxGapContent, matrixAA])
    pool.map_async(jarvis_filter_worker, work_list).get(9999999)
#        jarvis_filter_worker(infile, outfile, window, step, maxGapContent, matrixAA)
    for cur_og in og_list:
        jarvis_mask("%s/og_cds_%s.afa" % (indir, cur_og), "%s/og_pep_%s_WrongA.txt" % (outdir, cur_og), "%s/og_cds_%s.afa" % (outdir, cur_og), len_min)

def jarvis_mask(infile, jarvis_out, outfile, len_min):
    reader = open(jarvis_out, 'rU')
    filter_coords = {}
    for line in reader:
        if "END_OF_MASKS" in line:
            break
        cur_line = line.split()
        filter_coords[cur_line[0]] = ((int(cur_line[1])+1)*3, (int(cur_line[2])+1)*3)
    reader = SeqIO.parse(infile, format = 'fasta')
    seq_dic = {}
    for rec in reader:
        seq_dic[rec.id] = list(str(rec.seq))
    for seq_name, mask in filter_coords.items():
        seq_dic[seq_name][mask[0] : mask[1]] = "N"*(mask[1]-mask[0])
    outdic = {}
    for k, v in seq_dic.items():
#        if len(v) - (v.count("N") + v.count("-")) >= len_min:
        outdic[k] = "".join(v)
    outfile = open(outfile, 'w')
    for k, v in outdic.items():
        outfile.write(">%s\n%s\n" % (k, v))
    outfile.close()

def jarvis_filter_worker(param_list): #infile, outfile, window, step, maxGapContent, matrixAA):
    infile = param_list[0]
    print infile.split("_")[-1][0:-4]
    outfile = param_list[1]
    window = param_list[2]
    step = param_list[3]
    maxGapContent = param_list[4]
    matrixAA = param_list[5]
    wrong={} #store all windows for this alignment
    align = AlignIO.read(infile, "fasta")  
    out1 = open(outfile, 'w')
    results={}
    j=0
    while j<len(align[1].seq)-window:
        start=j+1
        end=j+window
        windowAlign=align[:, j:j+window] 
        scoresProb=getAveProbPerWindowPerSeq(windowAlign,maxGapContent)
        scoresMatrixPerCol=getColumnDistScoreToClosest(windowAlign, matrixAA) 
        for record in align:
            name=record.id 
            species=name.split("_")[-1]
            
            seq=str(record.seq)
            startSeq=len(re.sub("-", "", seq[:j]))+1
            endSeq=len(re.sub("-", "", seq[:j+window]))
        
            id=name.split("_")[-1]
            
            P=scoresProb[name][0]
            Z=scoresProb[name][1]
            
            D=scoresMatrixPerCol[name][0]
            Z62=scoresMatrixPerCol[name][1]
            
            if P!="NA" and float(P)<0.7: #ignore windows with too many gaps or those that are very conserved
            
                #model for ALL in W15S3 but actually works better also with W12
                result=1+2.21*P+0.36*Z62+1.32*Z+0.06*D+(-1.5)*D*P+(-0.70)*Z*P+(-0.10)*Z*Z62                
                if result > 0: # it's WRONG or at least very divergent according to the alignment
                    
                    if not wrong.has_key(name):
                        wrong[name]=[]
                
                    if len(wrong[name])>0:

                        if wrong[name][-1][1]>=startSeq: #check for coordinates overlap

                            wrong[name][-1][1]=endSeq
                            wrong[name][-1][3]=end
                            wrong[name][-1][-1]+=1 #n.o window overlapped
                            wrong[name][-1][-2]+=result
                        else:
                            wrong[name].append([startSeq, endSeq, start,end,result,1])
                    else:
                        wrong[name].append([startSeq, endSeq, start,end,result,1])
        
        j+=step
    for name in wrong:
        info=wrong[name]
        
        i=0
        while i <len(info):
            infoA=info[i]
            ave=infoA[4]
            nWindows=infoA[5]
            aveResults=round(ave/nWindows, 2)
            infoA[-2]=aveResults
            
            if aveResults > 3:
                
                infoA[-1]="C"
                a="\t".join(["%s" % el for el in infoA])      
                out1.write("%s\t%s\n"%(name, a))
                        
            i+=1
    out1.write("END_OF_MASKS")
    out1.close()

def getAveProbPerWindowPerSeq(align, pGaps):
    
    scorePerSpeciesPerWindow={} # [species]=[Z-score w1, Z-score w2...]
    ####################################################################################
    # amino acid frequencies, gaps and X will be given a zero prob
    aminoacids=["X","A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-"]

    #start up the aa freqs dicts
    aaFreqsDict={}          
    for element in aminoacids:

        aaFreqsDict[element]=[] 
    

    currentWindowSpeciesProbs={} # [species]=[probsite1, probsite2...]
    aveWindowSpeciesProbs={} # [species]=windowprob
    window=len(align[1].seq)
    
    i=0
    while i<window:
        
        col=align[:, i]
        cleancol=col.replace("-", "").replace("X", "")
        
        # populate the amino acids frequencies per site per column 
        for key in aaFreqsDict:

            if key!="-" and key!="X" and len(cleancol)>0 :
                aaFreqsDict[key]=round(cleancol.count(key)*1.00/len(cleancol), 10)

            else:
                aaFreqsDict[key]="NA"
                
        # save the prob of a species site aa
        for name in align:
            species=name.id
            currentAA=name.seq[i]
            siteProb=aaFreqsDict[currentAA]
            if not currentWindowSpeciesProbs.has_key(species):
                currentWindowSpeciesProbs[species]=[]        
                
#            if siteProb!="NA":#
            currentWindowSpeciesProbs[species].append(siteProb)
#            else:
#                currentWindowSpeciesProbs[species].append(siteProb)#####
                
        i+=1
        
    #calculate ave prob per window  per seq
    allAveWindow=[]
    
    for species in currentWindowSpeciesProbs:
        probs=currentWindowSpeciesProbs[species]
        probsNumeric=[x for x in probs if x!="NA" ]
        if len(probsNumeric)>pGaps*window: #if not mostly-gap seq
        
            aveProb=array(probsNumeric).mean()
            aveWindowSpeciesProbs[species]=round(aveProb, 2)
            allAveWindow.append(aveProb)
            
        else:
            aveWindowSpeciesProbs[species]="NA"
            
    # calculate Z-score
    if allAveWindow>0: #not a gap-only window
        aveWindow=array(allAveWindow).mean()
        stdevWindow=array(allAveWindow).std()
        
    for species in aveWindowSpeciesProbs:
        probSpecies=aveWindowSpeciesProbs[species]
        Z=0
        if not probSpecies=="NA": #not a gap-only window
            if not stdevWindow<0.001: #very conserved window
                dev=abs(aveWindow-probSpecies)
                Z=round(dev/stdevWindow, 2)
            else:
                Z=0
        else:
            Z="NA"
            
        if not scorePerSpeciesPerWindow.has_key(species):
            scorePerSpeciesPerWindow[species]=[]
            
        score=[probSpecies,Z]
        scorePerSpeciesPerWindow[species]=score
      
    return scorePerSpeciesPerWindow

def readMatrix(nameMatrix):
    #takes an half similarity matrix and turns into a dict

    input=open(nameMatrix)
    lines=input.readlines()

    names=lines[0].rstrip().split()

    aaDict={} #aaDict[(aa1,aa2)]=value

    i=1
    while i<len(lines):
        info=lines[i].rstrip().split()
        aa1=info[0]

        j=1
        while j<len(info):
            aa2=names[j-1]
            value=int(info[j])
            pair=(aa1,aa2)
            aaDict[pair]=value
            pairI=(aa2,aa1)
            aaDict[pairI]=value
            j+=1
            
        i+=1
    return aaDict

def scorePairAlignedSeqsWithMatrix(seq1,seq2,matrixAA):
    #takes 2 aa seqs as strings with the same size, 
    #and calculates the score of aligned aas ignoring gaps

    totalScore=[]
    gaps=0
    i=0
    while i<len(seq1):
        if seq1[i]!="-" and seq2[i]!="-" and seq1[i]!="X" and seq2[i]!="X":
            totalScore.append(matrixAA[(seq1[i],seq2[i])])
        else:
            gaps+=1
        i+=1
        
    if gaps==len(seq1):
        
        total="NA"
    
    else:
        total=array(totalScore).mean()
    
    return total

def getPairwiseScore(windowAlign, matrixAA):
    
    scoreSlice={} # [species]=lowestScore
    allSeqs={}
    forAverage=[]
    
    for entry in windowAlign:
        name=entry.id
        allSeqs[name]=str(entry.seq)

    for name in allSeqs:
        allScores=[]
        others=allSeqs.keys()
        others.remove(name)
        
        for alt in others:
            seq1=allSeqs[name]
            seq2=allSeqs[alt]
            score=scorePairAlignedSeqsWithMatrix(seq1, seq2, matrixAA)
            
            if score!="NA":
                allScores.append(score)
            
        if len(allScores)>0:
            score=round(allScores[-1], 2)
            forAverage.append(score)
            scoreSlice[name]=[score]
        
    aveWindow=array(forAverage).mean()
    stdevWindow=array(forAverage).std()
        
    for species in allSeqs:
        if scoreSlice.has_key(species):
            score=round(scoreSlice[species][0], 2)
            
            if not stdevWindow<0.001: #very conserved window
                dev=abs(aveWindow-score)
                Z=round(dev/stdevWindow, 2)
                
            else:
                Z=0
                
            scoreSlice[species].append(Z)
        else:
            scoreSlice[species]=["NA", "NA"]
        
    return scoreSlice

def extractSeqs(fileLocation, alnName, listSpecies):
    
    import re
    seqs={}
    ids={}

    align = open("%s"%(fileLocation))
    info=align.readline()
    while info:
        if re.match(">", info):
            key=info.rstrip()[1:]
            seqs[key]=""
            
            species=key.split("_")[-1]
            ids[species]=key
            
        else:
            seqs[key]+=info.rstrip()
            
        info=align.readline()
        
    for name in listSpecies:
        if ids.has_key(name):
            seqs.pop(ids[name])
            
    out1=open("%s_clean.temp"%(alnName), "w")
    for key in seqs:
        out1.write(">%s\n%s\n"%(key, seqs[key]))
    out1.close()

def getColumnDistScoreToClosest(windowAlign, matrixAA):
    
    scoreSlice={} # [species]=lowestScore
    allSeqs={}
    forAverage=[]
    allScores=[]
    
    for entry in windowAlign:
        name=entry.id
        allSeqs[name]=str(entry.seq)

    for name in allSeqs:
        seq1=allSeqs[name]
        others=allSeqs.keys()
        others.remove(name)
        windowScore=[]
        
        i=0
        while i<len(seq1):
            scoreCol=[]
            #score each position of the alignment 
            for alt in others:
            
                seq2=allSeqs[alt][i]
                score=scorePairAlignedSeqsWithMatrix(seq1[i], seq2, matrixAA)
                if score!="NA":
                    scoreCol.append(score)
            
            if len(scoreCol)>0: 
                scoreCol.sort()
                windowScore.append(scoreCol[-1]) #get the highest score for this column  
                
            i+=1
        
        if len(windowScore)>0:
            scoreAVEwindow=array(windowScore).mean()
            scoreAVEwindow=round(scoreAVEwindow, 2)
            
            allScores.append(scoreAVEwindow)
            scoreSlice[name]=[scoreAVEwindow]   
    
    aveWindow=array(allScores).mean()
    stdevWindow=array(allScores).std()
    
    for species in allSeqs:
        if scoreSlice.has_key(species):
            score=scoreSlice[species][0]
            
            if not stdevWindow<0.001: #very conserved window
                dev=abs(aveWindow-score)
                Z=round(dev/stdevWindow, 2)
            else:
                Z=0
                
            scoreSlice[species].append(Z)
        else:
            scoreSlice[species]=["NA", "NA"]
        
    return scoreSlice
