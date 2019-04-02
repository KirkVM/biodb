import re,pickle,os
import numpy as np
from kmslib.hmmerkools import hmmsearchparser
from . import seqdbutils

class EndRangeProb:
    def __init__(self,dpm,seqlen,hmmsize):
        self.hmmpositions=None
        self.seqpositions=None
        self.envfractions=None
        self.initial_logistic=None
        self.updated_prob=None
        self.extrange=int(0.26*hmmsize)
    def update_prob(self,curatedEA):
        pass

##HMMER numbering starts at 1!!!
class NTRangeProb(EndRangeProb):
    def __init__(self,dpm,seqlen,hmmsize):
        super().__init__(dpm,seqlen,hmmsize)
        #HMMER is 1-indexed so if aln starts at first position ntalign_min,ntalign_max=0,1
        ntalign_min=max(0, (dpm.align_start-dpm.hmm_start)-self.extrange)
        ntalign_max=dpm.align_start 
        ntalign_span=ntalign_max-ntalign_min
        #seqpositions will now be 0-indexed
        self.seqpositions=np.array(range(ntalign_min,ntalign_max))
        #keeping HMMER positions 1-indexed
        self.hmmpositions=np.array(range(dpm.hmm_start+1-ntalign_span,dpm.hmm_start+1))
        self.envfractions= -(self.hmmpositions-1)/hmmsize

        gamma=1
        env_frac_remaining=(dpm.hmm_start-1)/hmmsize
        fx=-7.5*self.envfractions*np.exp(1-env_frac_remaining)
        for x in range(len(self.envfractions)):
            if self.envfractions[x]<=-0.05*env_frac_remaining:
                fx[x]+=-50*(self.envfractions[x]+0.05*env_frac_remaining)
        self.initial_logistic=1/(1+np.exp(-gamma*(fx)))

##HMMER numbering starts at 1!!!
class CTRangeProb(EndRangeProb):
    def __init__(self,dpm,seqlen,hmmsize):
        super().__init__(dpm,seqlen,hmmsize)
        #HMMER is 1-indexed so if aln stops at last position ntalign_min,ntalign_max=seqlen-1,seqlen
        ctalign_min=dpm.align_stop-1
        ctalign_max=min(seqlen, dpm.align_stop + (hmmsize-dpm.hmm_stop) + self.extrange)
        ctalign_span=ctalign_max-ctalign_min
        #seqpositions will now be 0-indexed
        self.seqpositions=np.array(range(ctalign_min,ctalign_max))
        #keeping HMMER positions 1-indexed
        self.hmmpositions=np.array(range(dpm.hmm_stop,dpm.hmm_stop+ctalign_span))
        self.envfractions=-((hmmsize-self.hmmpositions)/hmmsize)
        
        gamma=1
        env_frac_remaining=(hmmsize-dpm.hmm_stop)/hmmsize
        fx=-7.5*self.envfractions*np.exp(1-env_frac_remaining)
        for x in range(len(self.envfractions)):
            if self.envfractions[x]<=-0.05*env_frac_remaining:
                fx[x]+=-50*(self.envfractions[x]+0.05*env_frac_remaining)
        self.initial_logistic=1/(1+np.exp(-gamma*(fx)))
    def get_boundary(self):
        pass

def build_cctable(dbpath,hmmsearchfpath,pfamcode):#,seqfpath,seqformat='fasta'):
    hsrp=hmmsearchparser.HMMERSearchResParser()
    hsrp.file_read(hmmsearchfpath)
    #motifHT=hsrp.get_motifHT()
    protHT=hsrp.get_protHT() #optional...can add this to then get picture of that sequence's matches
    filteredHT={}
    accRE=re.compile("(\S+)\.(\d+)")
    for k in protHT.keys():
        nrmlaccobj=accRE.match(k)
        if nrmlaccobj:
            pacc,pvrsn=nrmlaccobj.groups()
        else: #can't split into <acc>.<vrsn>
            pacc=k
        for spm in protHT[k]:
            if spm.pacc==pfamcode:
                filteredHT[pacc]=spm
    print(f'{len(filteredHT.keys())} sequences in {os.path.basename(hmmsearchfpath)} match {pfamcode}')

    #my_motif=motifHT[pfamcode] #my motif is a list of spms (SeqProfileMatch objects)
    hmmsfile_mtime=os.stat(hmmsearchfpath).st_mtime
    conn=seqdbutils.gracefuldbopen(dbpath) #open the db after HMMER file read-in goes ok
    c=conn.cursor()
    c.execute('''SELECT acc FROM PROTEINGBS WHERE pklgbsr''')
    pgbsraccs=[x['acc'] for x in c.fetchall()]
    absent_entries=list(set(pgbsraccs).difference(filteredHT.keys()) )
    print(f'{len(absent_entries)} entries in PROTEINGBS but not {os.path.basename(hmmsearchfpath)}')
    c.execute('''CREATE TABLE IF NOT EXISTS CCDATA (acc text, version text, pfamcode text,pfamvrsn int,\
                seq_checksum, hmmsfname text, hmmsf_mtime int, spm glob, ntrp glob, ctrp glob)''')
    #seq_checksum included so could compare from PROTEINGBS table perspective whether need to re-run hmmsearch 
    c.execute('''SELECT acc FROM CCDATA''')
    #find accs that need to be added-
    acc_indb=[x['acc'] for x in c.fetchall()]

    acc2add=list(set(filteredHT.keys()).difference(acc_indb))
    print(f'{len(acc_indb)} entries already in CCDATA, {len(acc2add)} new entries to add')

    for acc in acc2add:
        spm=filteredHT[acc]
        proteinacc=spm.prot_uval
        nrmlaccobj=accRE.match(proteinacc)
        if nrmlaccobj:
            pacc,pvrsn=nrmlaccobj.groups()
        else: #can't split into <acc>.<vrsn>
            pacc=proteinacc
            pvrsn=None
        c.execute('''SELECT * FROM PROTEINGBS WHERE acc = (?)''',(pacc,) )
        gbdbentry=c.fetchone()
        seqlen=len(pickle.loads(gbdbentry['pklgbsr']).seq)
        ntrp=NTRangeProb(spm.dpms_[0],seqlen,spm.psize)
        ctrp=CTRangeProb(spm.dpms_[0],seqlen,spm.psize)
        new_tuple=(pacc,pvrsn,pfamcode,spm.paccvrsn,gbdbentry['seq_checksum'],os.path.basename(hmmsearchfpath),\
                   hmmsfile_mtime,pickle.dumps(spm),pickle.dumps(ntrp),pickle.dumps(ctrp))
        c.execute('''INSERT INTO CCDATA VALUES (?,?,?,?,?,?,?,?,?,?)''',new_tuple)
    conn.commit()
    conn.close()
    #ccvals=find_motifalign(my_motif)