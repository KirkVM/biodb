import os,sys,pickle
import pandas as pd #need xlrd as well
from . import seqdbutils


def build_curatedxtaltable(dbpath,xtalxlfpath,sheet_name="Sheet1"):
    xtaldf=pd.read_excel(xtalxlfpath,sheet_name=sheet_name,dtype={'flagged':str})
    xtalfname=os.path.basename(xtalxlfpath)
    xtalfile_mtime=os.stat(xtalxlfpath).st_mtime
    conn=seqdbutils.gracefuldbopen(dbpath)
    c=conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS CURATEDXTALS (acc text, xtalfname text, xtalf_mtime int, \
                cazypdblist text, pdbid text, pdbchain text, ntccpos int, ctccpos int, ntccseq text, ctccseq text,\
                estimatedcc_hmmstart int,estimatedcc_hmmstop int,estimatedcc_seqstart int,estimatedcc_seqstop)''')
    for _,xrow in xtaldf.dropna(subset=['curated_pdb']).iterrows():
        assert(type(xrow['cc_startnum']) in (float,int)),"cc_startnum type not numeric"
#        print(type(xrow.flagged))
        if str(xrow['flagged']).lower()=='yes': #not a confident assessment, skip for now
            print('flagged',xrow['gb_accession'],'- passing')
            continue
        new_tuple=(xrow['gb_accession'],xtalfname,xtalfile_mtime,xrow['PDBS'],xrow['curated_pdb'],\
                 xrow['chain'],int(xrow['cc_startnum']),int(xrow['cc_stopnum']),xrow['cc_startseq'],xrow['cc_stopseq'],\
                 None,None,None,None)
        c.execute('''INSERT INTO CURATEDXTALS VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',new_tuple)
    conn.commit()
    conn.close()


def mergein_predictedccdata(dbpath):
    conn=seqdbutils.gracefuldbopen(dbpath)
    seqdbutils.check_tables_exist(conn,['CURATEDXTALS','CCDATA'])
    c=conn.cursor()

    #c.execute('''SELECT CCDATA.acc,CCDATA.ccb,CCDATA.ccstart,CCDATA.ccstop,CURATEDXTALS.ntccpos,CURATEDXTALS.ctccpos FROM CCDATA INNER JOIN CURATEDXTALS ON CCDATA.acc = CURATEDXTALS.acc''')
    c.execute('''SELECT * FROM CCDATA INNER JOIN CURATEDXTALS,PROTEINGBS ON CCDATA.acc = CURATEDXTALS.acc''')
    cxccs=c.fetchall()
    for cxcc in cxccs:
        curacc=cxcc['acc']
        c.execute('''SELECT pklgbsr FROM PROTEINGBS WHERE acc=(?)''',(curacc,))
        gbsr=pickle.loads(c.fetchone())
        myseq=gbsr.seq
        #NOW MATCH the seq to see if this works... 
#        curccb=pickle.loads(cxcc['ccb'])
#        print(curccb.)
        try:
            #print(cxcc['acc'],cxcc['ccstart'],cxcc['ccstop'],cxcc['ntccpos'],cxcc['ctccpos'])
            print(cxcc['acc'],cxcc['ccstart'],cxcc['ccstop'],cxcc['ntccpos'],cxcc['ctccpos'])
            print(cxcc['ntccseq'],cxcc['ctccseq'])
        except:
            conn.close()
            sys.exit()
    conn.close()

#tests-numtypes coming out of db
#unique list of accs in PROTEINGBS