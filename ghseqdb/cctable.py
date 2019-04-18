import re,pickle,os,sys
import numpy as np
from kmslib.hmmerkools import hmmsearchparser
from . import seqdbutils

class CCBoundaries:
    def __init__(self):
        self.curated_pwbounds0


#TODO- addin a per-GH yaml file with appropriate settings for conf_envfrac etc.
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
    #now read in curate_alndf for each sequence?
    conn=seqdbutils.gracefuldbopen(dbpath) #open the db after HMMER file read-in goes ok
    seqdbutils.check_tables_exist(['CCDATA'])
    c=conn.cursor()
    c.execute('''SELECT acc FROM PROTEINGBS WHERE pklgbsr''')
    pgbsraccs=[x['acc'] for x in c.fetchall()]
    absent_entries=list(set(pgbsraccs).difference(filteredHT.keys()) )
    print(f'{len(absent_entries)} entries in PROTEINGBS but not {os.path.basename(hmmsearchfpath)}')
    c.execute('''CREATE TABLE IF NOT EXISTS CCDATA (acc text, version text, pfamcode text,pfamvrsn int,\
                seq_checksum, hmmsfname text, hmmsf_mtime int, spm glob, ccb glob, ccstart int, ccstop int,
                ntenvfrac REAL, ctenvfrac REAL, verified_status text)''')
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
        ccbounder=CCBounder(spm.dpms_[0],seqlen,spm.psize)

        new_tuple=(pacc,pvrsn,pfamcode,spm.paccvrsn,gbdbentry['seq_checksum'],os.path.basename(hmmsearchfpath),\
                   hmmsfile_mtime,pickle.dumps(spm),pickle.dumps(ccbounder),int(ccbounder.ccseqstart),int(ccbounder.ccseqstop),\
                   ccbounder.nthmm_envfrac,ccbounder.cthmm_envfrac,"unknown")
#        print(acc,ccbounder.ccseqstart,ccbounder.ccseqstop,ccbounder.nthmm_envfrac,ccbounder.cthmm_envfrac)
        c.execute('''INSERT INTO CCDATA VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',new_tuple)
        ##ADD something that indicates coverage at start???
    conn.commit()
    conn.close()
    #ccvals=find_motifalign(my_motif)


def extractsrs_cc(dbpath,format='fasta'):
    """returns a file containing all the sequence files in PROTEINGBS table"""
    conn=seqdbutils.gracefuldbopen(dbpath)
    seqdbutils.check_tables_exist(conn,['CCDATA','PROTEINGBS'])
    
    c=conn.cursor()
    c.execute('''SELECT PROTEINGBS.acc,pklgbsr,ccstart,ccstop,ntenvfrac,ctenvfrac FROM PROTEINGBS INNER JOIN CCDATA ON PROTEINGBS.acc=CCDATA.acc''')
    allrows=c.fetchall()
    print(f'found {len(allrows)} non-null entries')
    outsrs=[]
    for gbrow in allrows:
        if gbrow['ntenvfrac']>0.8 and gbrow['ctenvfrac']>0.8:
            gbsr=pickle.loads(gbrow['pklgbsr'])
            newsr=gbsr[gbrow['ccstart']:gbrow['ccstop']]
            newsr.id=gbrow['acc']
            outsrs.append(newsr)
    
    conn.close()
    return outsrs
