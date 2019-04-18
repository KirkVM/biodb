import pickle
from Bio import SeqIO,pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import numpy as np

class PWAInfo:
    def __init__(self,dbacc,curateacc,bpalns,gap_penalty,ext_penalty,penalize_end_gaps):
        ''''sets up aln info

        Arguments:
            dbacc- accession code of sequence to be investigated
            curateacc- acc code of curated sequence (length truncd to known domain(s) size)
            bpalns: a biopython alignment object (list of alignments-
                    -alignment ordering should be performed according to pw(dbseq,curateseq) (can I verify this??)
        '''
        self.dbacc=dbacc
        self.curateacc=curateacc
        self.alns=bpalns #a biopython pw aligns result
        self.gap_penalty=gap_penalty
        self.ext_penalty=ext_penalty
        self.penalize_end_gaps=penalize_end_gaps
        self.dbseqlen=None
        self.curateseqlen=None
        self.seqaln_density=None
        self.vseq_curatedensity=None
        self.score=None
        self.starts=[]
        self.stops=[]
        self.avg_start=None
        self.avg_stop=None
        self.median_start=None
        self.median_stop=None
        self.initialize_info()
        if len(self.alns)>0:
            self.buildinfo()
    def initialize_info(self):
        if len(self.alns)==0:
            self.score=np.nan
        else: #use aln[0] to set up
            aln=self.alns[0]
            self.score=aln[2]
            self.dbseqlen=len(aln[0])-aln[0].count('-')
            self.curateseqlen=len(aln[1])-aln[1].count('-')
            self.seqaln_density=np.zeros((self.dbseqlen))
            self.vseqaln_curatedensity=np.zeros((self.dbseqlen))
            self.avg_start=0.0
            self.avg_stop=0.0
    def buildinfo(self):
        for aln in self.alns:
            dbalnseq=aln[0]
            curatealnseq=aln[1]
            db_pos=-1
            curate_pos=-1
            aln_pos=-1
            aln_length=len(aln[0])
            for aln_pos in range(aln_length):
                dbres=dbalnseq[aln_pos]
                curateres=curatealnseq[aln_pos]
                if dbres!='-':
                    db_pos+=1
                    aln_nbrs=dbalnseq[max(0,aln_pos-3):min(aln_length,aln_pos+4)]
                    aln_nbrsize=len(aln_nbrs)-1 #don't include current position
                    resdensity=(aln_nbrsize-aln_nbrs.count('-'))/aln_nbrsize
                    self.seqaln_density[db_pos]+=resdensity
                    
                    valn_nbrs=curatealnseq[max(0,aln_pos-3):min(aln_length,aln_pos+4)]
                    valn_nbrsize=len(valn_nbrs) #do include current position
                    vresdensity=(valn_nbrsize-valn_nbrs.count('-'))/valn_nbrsize
                    self.vseqaln_curatedensity[db_pos]+=vresdensity
                if curateres!='-':
                    curate_pos+=1
                if curate_pos==0:
                    self.starts.append(db_pos)
                if curate_pos==self.curateseqlen-1:
                    self.stops.append(db_pos)
        self.avg_start=np.mean(self.starts)
        self.avg_stop=np.mean(self.stops)
        self.median_start=np.median(self.starts)
        self.median_stop=np.median(self.stops)
        self.seqaln_density/=len(self.alns)
        self.vseqaln_curatedensity/=len(self.alns)

def biopython_pwaln(dbsr,curatesr,gap_penalty=-10,ext_penalty=-0.5,penalize_end_gaps=(True,False)):
    alns=pairwise2.align.globalds(dbsr.seq,curatesr.seq,blosum62,gap_penalty,ext_penalty,penalize_end_gaps=penalize_end_gaps)
    pwainfo=PWAInfo(dbsr.id,curatesr.id,alns,gap_penalty=gap_penalty,ext_penalty=ext_penalty,penalize_end_gaps=penalize_end_gaps)
    return pwainfo

def get_curate_alstats(dbgbsr,cgbsrs,gap_penalty=-10,ext_penalty=-0.5,penalize_end_gaps=(True,False)):
    pwainfos=[]
    for cgbsr in cgbsrs:
        curpwainfo=biopython_pwaln(dbgbsr,cgbsr,gap_penalty=gap_penalty,ext_penalty=ext_penalty,penalize_end_gaps=penalize_end_gaps)
        pwainfos.append(curpwainfo)
    return pwainfos

def multi_curate_alstats(dbgbsrs,cgbsrs,gap_penalty=-10,ext_penalty=-0.5,penalize_end_gaps=(True,False)):
    allpwainfos=[]
    for dbgbsr in dbgbsrs:
        pwainfos=get_curate_alstats(dbgbsr,cgbsrs,gap_penalty=gap_penalty,ext_penalty=ext_penalty,penalize_end_gaps=penalize_end_gaps)
        allpwainfos.append(pwainfos)
    return allpwainfos
    

#def dostuff(dbfile,dbnames,curatefile,curatenames,numsplits,split):


