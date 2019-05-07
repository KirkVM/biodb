import sqlite3,re,datetime,pickle
from . import seqdbutils
from scrapers import cazydbscrapers
from kmslib.hmmerkools import hmmsearchparser

def build_hmmertable(dbpathstr,hmmsrchfpathstr):
    """adds HMMERSEQDATA table with hmmersearch. hmmsearch results file typically built using
    
    hmmsearch --noali -E 1e-1 --incE 1e-1 [pfamfile] [seqfile] > [outfile]

    Argumenst:
        dbpathstr: path to db
        hmmsrchfpathstr: path to hmmersearch results
    """    

    conn=seqdbutils.gracefuldbopen(dbpathstr) 
    seqdbutils.check_tables_exist(conn,['PROTEINGBS'])
    c=conn.cursor()
    c.execute('''DROP TABLE IF EXISTS HMMERSEQDATA''')
    c.execute('''CREATE TABLE HMMERSEQDATA (acc text, version text, seq_checksum text, \
                motifacc text, motifvrsn txt, motifsize int, align_start int, align_stop int, \
                env_start int, env_stop int, hmm_start int, hmm_stop int, i_eval real, dpm glob)''')

    hsrp=hmmsearchparser.HMMERSearchResParser()
    hsrp.file_read(hmmsrchfpathstr)
    protHT=hsrp.get_protHT() #optional...can add this to then get picture of that sequence's matches
    for proteinacc in protHT.keys():
        c.execute('''SELECT * FROM PROTEINGBS WHERE acc=(?)''',(proteinacc,))
        pgbrow=c.fetchone()
        pvrsn=pgbrow['version']
        pcksum=pgbrow['seq_checksum']
        for spm in protHT[proteinacc]:
            motifacc=spm.pacc
            motifvrsn=spm.paccvrsn
            motifsize=spm.psize
            for dpm in spm.dpms_:
                #align_start_state,align_stop_state,env_start_state,env_stop_state
                pkldpm=pickle.dumps(dpm)
                newtuple=(proteinacc,pvrsn,pcksum,motifacc,motifvrsn,motifsize,dpm.align_start,dpm.align_stop,\
                        dpm.env_start,dpm.env_stop,dpm.hmm_start,dpm.hmm_stop,dpm.i_evalue,pkldpm)
                c.execute('''INSERT INTO HMMERSEQDATA VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',newtuple)
    conn.commit()
    conn.close()