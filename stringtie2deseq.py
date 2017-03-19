#!/usr/bin/python

"""
Script to convert the transcript level output files made by stringtie for ballgown into files that will work with DEseq2.
"""

import os
import sys
import re
import pandas as pd
from collections import defaultdict


#----------------------------------------------------------------------------------------------------------------------
class Transcript:
    def __init__(self, idx, chr, start, end, strand):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand

        s = re.sub(r"[\"\s]+", "", idx)

        self.id = s
        self.geneId = re.sub(r"\.\d$", "", s)
        self.ucounts = 0
        self.exons = dict() # k= exon ID from 'e_data', v = unique counts for this exon

    def addExon(self, e_id, ucounts):
        self.exons[e_id] = ucounts

    def calcUcounts(self):
        for k in self.exons:
            self.ucounts += self.exons[k]


#----------------------------------------------------------------------------------------------------------------------
def process(fp, tsDict):
    """
    Each folder contains:
        e2t.ctab: exon-to-transcript mapping file
        e_data.ctab: exon level counts
        t_data.ctab: transcript level counts
    """
    ret = dict() # k = transcript ID, v = transcript read count

    ## read in the exon-to-transcript mappping
    e2t = defaultdict(list) # k = transcript number v = list of exon Ids
    e2tF = fp + "/e2t.ctab"
    E2T = open(e2tF, 'r')
    for line in E2T:
        line = line.strip()
        if line.startswith("e_id"): continue;
        (e, t) = line.split("\t")
        e2t[t].append(e)
    E2T.close()

    ## read in the exon read counts counts
    exons = dict()  # k = exon Id, v = counts
    eF = fp + "/e_data.ctab"
    E = open(eF, 'r')
    for line in E:
        line = line.strip()
        if line.startswith("e_id"): continue;
        e_id = line.split("\t")[0]
        e_count = int(line.split("\t")[6])
        exons[e_id] = e_count
    E.close()

    ## read in the transcript summary file
    tsF = fp + "/t_data.ctab"
    T = open(tsF, 'r')
    for line in T:
        line = line.strip()
        if line.startswith("t_id"): continue;
        t_id = line.split("\t")[0]
        tsK = line.split("\t")[5]

        ## get all the exons that map to this transcript
        if tsK in tsDict:
            curTS = tsDict[ tsK ]
            curExons = e2t[ t_id ]
            for e in curExons: curTS.addExon(e, exons[e]);
            curTS.calcUcounts()
            ret[ tsK ] = curTS.ucounts
        else:
            ret[ tsK ] = 0
    T.close()

    return(ret)



def parseGTF(inF):
    ret = dict()
    F = open(inF, "r")
    for line in F:
        line = line.strip()
        if line.startswith("#"): continue;
        data = line.split("\t")

        if data[2] == 'transcript':
            chr = data[0]
            s = data[3]
            e = data[4]
            st = data[6]
            id = data[8].split(";")[1].replace("transcript_id ", "")
            ts = Transcript(id, chr, s, e, st)
            ret[ts.id] = ts
    F.close()
    return(ret)

#----------------------------------------------------------------------------------------------------------------------


if len(sys.argv) < 4:
    sys.stderr.write("\nUSAGE: python " + sys.argv[0] + " GTF_ancher_file ballgown_data_folder output_file_name\n\n")
    sys.exit(0)

gtf=os.path.abspath(sys.argv[1])
ballgownDir=os.path.abspath(sys.argv[2])
outName=sys.argv[3]

print "GTF = " + os.path.abspath(sys.argv[1])
print "Target Dir = " + ballgownDir

#####################
## Main part starts!
#####################
tsDict = parseGTF(gtf)

outDF = pd.DataFrame(None) ## holds final output

for fn in os.listdir(ballgownDir):

    if fn.endswith(".gtf"): continue;

    d = process("/".join([ballgownDir,fn]), tsDict)

    df = pd.DataFrame(d.items(), columns=['transcriptId', 'ucounts'])
    df['sample'] = fn
    outDF = outDF.append(df)


out2 = outDF.pivot(index="transcriptId", columns="sample", values="ucounts")
out2.to_csv(outName, sep="\t", header=True, index=True)

sys.stderr.write("\nDESeq2 file written to " + os.path.abspath(outName) + "\n\n")