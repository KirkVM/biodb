import os,sys,pickle,re,difflib
import pandas as pd #need xlrd as well
from pathlib import Path
from Bio.PDB import MMCIFParser,PPBuilder
from . import seqdbutils

def check_cxstatus(c,xrow,xtalfile_mtime):
    c.execute('''SELECT * FROM CURATEDXTALS WHERE pdbid=(?)''',(xrow['curated_pdb'],))
    prevrow=c.fetchone()
    #print(prevrow['alt_startnum'].isnull(),xrow['alt_startnum'].isnull())
    if xrow['flagged']: #not a confident assessment, skip for now
        print('flagged',xrow['gb_accession'],'- passing')
        if prevrow is not None:
            print(f'Due to flagged status, dropping existing entry {xrow["curated_pdb"]}')
            return 'delete'
        else:
            return 'skip'

    if prevrow is None:
        return 'absent'
    else:
        if (prevrow['ntccpos']==xrow['cc_startnum'] and prevrow['ctccpos']==xrow['cc_stopnum'] and
            prevrow['ntccseq']==xrow['cc_startseq'] and prevrow['ctccseq']==xrow['cc_stopseq'] and
            prevrow['acc']==xrow['gb_accession'] and prevrow['cazypdblist']==xrow['PDBS'] and
            prevrow['enable_fuzzy']==xrow['enable_fuzzy']):
            return 'exists'
        else:
            assert (xtalfile_mtime>prevrow['xtalf_mtime']),"replacing with val from older excel file?"
            return 'update_required'
    return 'prev_status_unclear'

def align_proteingbs_seq(conn,xtalsdbpath):
    #check PROTEINGBS exists then 
    seqdbutils.check_tables_exist(conn,['PROTEINGBS','CURATEDXTALS'])

    assert (xtalsdbpath.exists()),f'not finding xtalsdbpath {xtalsdbpath}, check keyword args'
    xdbconn=seqdbutils.gracefuldbopen(xtalsdbpath)
    seqdbutils.check_tables_exist(xdbconn,['XTALS'])
    xc=xdbconn.cursor()
    c=conn.cursor()
    
    mmcparser=MMCIFParser()
    ppb=PPBuilder()
    
    c.execute('''SELECT * FROM CURATEDXTALS''')
    curatedrows=c.fetchall()
    for curatedrow in curatedrows:
        core_start=None
        core_stop=None
        c.execute('''SELECT * FROM PROTEINGBS WHERE acc=(?)''',(curatedrow['acc'],))
        pgbrow=c.fetchone()
        if (pgbrow['seq_checksum']==curatedrow['pgbsr_seqchecksum'] and \
            curatedrow['pgbsr_ccstart'] is not None and \
            curatedrow['pgbsr_ccstop'] is not None):
            continue #skip this one!
        ntcc=curatedrow['ntccseq']
        ctcc=curatedrow['ctccseq']
        inter_seqlen=curatedrow['ctccpos']-len(curatedrow['ctccseq']) - (curatedrow['ntccpos']+len(curatedrow['ntccseq'])) + 1
        inter_max=inter_seqlen+2
        inter_min=inter_seqlen-2
        seqrestring=f'({ntcc}.{{{inter_min},{inter_max}}}{ctcc})'
#        print(seqrestring)
        seqRE=re.compile(seqrestring) 
        pgbseq=pickle.loads(pgbrow['pklgbsr']).seq
        pgbmatch=seqRE.search(str(pgbseq))
        if pgbmatch:
            #coreseqstr=pgbmatch.group(1)
            start,stop=pgbmatch.span()
            core_start=start
            core_stop=stop-1
            print(f'found {curatedrow["pdbid"]} the easy way by matching genbank sequence {curatedrow["acc"]}:\
                 ({core_start},{core_stop})')
            #print(curatedrow['pdbid'],start,stop-1)
        else:
            print(f'-->match to genbank seq failed for {curatedrow["pdbid"]}...trying xtalfile')
            xc.execute('''SELECT * FROM XTALS WHERE pdbid=(?)''',(curatedrow['pdbid'],))
            xtalrow=xc.fetchone()
            xtalfpath= xtalsdbpath.parent / xtalrow['relpath']
            xtalstruct=mmcparser.get_structure('aname',xtalfpath)
            xtalchain=xtalstruct#[0]#[curatedrow['pdbchain']]
            xtalppchain=ppb.build_peptides(xtalchain)
            xtalseqstr=str(xtalppchain[0].get_sequence())
            xstrmatch=seqRE.search(xtalseqstr)
            if xstrmatch:
                max_ratio=0
                start=None
                stop=None
                xtalcoreseqstr=xstrmatch.group(1)
                for strtcand in range(0,len(pgbseq)-inter_min ):
                    for stopcand in range(strtcand+inter_min,strtcand+inter_max):
                        #stopcand=min(stopcand,len(pgbseq))
                        sm=difflib.SequenceMatcher(a=xtalcoreseqstr,b=str(pgbseq),autojunk=False)
                        if sm.ratio()>max_ratio:
                            core_start=strtcand
                            core_stop=stopcand-1
                print(f'-->found {curatedrow["pdbid"]} by matching xtalpeptide to genbank sequence {curatedrow["acc"]}:\
                        ({core_start},{core_stop})')
            elif curatedrow['enable_fuzzy']:
                print(f'---->Now resorting to fuzzy matching for {curatedrow["pdbid"]}')
                cstarts=[]
                cstops=[]
                for strtcand in range(0,len(pgbseq)-inter_min+1 ):
                    sm_start=difflib.SequenceMatcher(a=ntcc,b=str(pgbseq)[strtcand:strtcand+5])
                    if sm_start.ratio()>0.59:
                        cstarts.append(strtcand)
                for stopcand in range(inter_min-1,len(pgbseq)):
                    sm_stop=difflib.SequenceMatcher(a=ctcc,b=str(pgbseq)[stopcand-4:stopcand+1])
                    if sm_stop.ratio()>0.59:
                        cstops.append(stopcand)
                foundit=False
                for cstart in cstarts:
                    if foundit: break
                    for cstop in cstops:
                        if foundit: break
                        if ((cstop-cstart)<(curatedrow['ctccpos']-curatedrow['ntccpos']+5) and\
                            (cstop-cstart)>(curatedrow['ctccpos']-curatedrow['ntccpos']-5)):
                            for xpos in range(0,len(pgbseq)-len(xtalseqstr)+5):
                                sm=difflib.SequenceMatcher(a=xtalseqstr,b=str(pgbseq)[xpos:xpos+len(xtalseqstr)])
                                if sm.ratio()>0.95:
                                    foundit=True
                                    core_start=cstart
                                    core_stop=cstop
                                    print(f'---->Fuzzy match worked for xtal {curatedrow["pdbid"]} - {curatedrow["acc"]}:\
                                         ({core_start},{core_stop})')
                                if foundit:
                                    break
            else:
                print('******************* * * * *  * * *  *  *   *   * * *    *   *   *  **  * *    *********************')
                print(f'-------could not match. Consider setting enable_fuzzy column to yes for {curatedrow["pdbid"]}-------')
                print('************  *  * *  *  * *  *  * *    *   *   * *   * *     *    *   *  * *  *******************')
        if core_start is not None and core_stop is not None:
            #do a fix that nudges to start or stop if w.in 3 of first/last position-
            if core_start<3:
                core_start=0
            if core_stop>len(pgbseq)-3:
                core_stop=len(pgbseq)-1
            update_tuple=(core_start,core_stop,pgbrow['seq_checksum'],curatedrow['acc'])
            c.execute('''UPDATE CURATEDXTALS SET pgbsr_ccstart=(?), pgbsr_ccstop=(?), pgbsr_seqchecksum=(?) WHERE acc=(?)''',\
                       update_tuple)
        else:
            print(f'+++++------FAILED TO FIND CORE START,STOP FOR {curatedrow["pdbid"]} ({curatedrow["acc"]})------++++++')
    xdbconn.close()




def build_curatedxtaltable(dbpathstr,xtalxlfpath,sheet_name="Sheet1",xtalsdb_relpath='up1',xtalsdb_name='EDATA.sql'):
    xtaldf=pd.read_excel(xtalxlfpath,sheet_name=sheet_name,dtype={'flagged':str})
    xtaldf['enable_fuzzy'].fillna('no',inplace=True)
    xtaldf['flagged'].fillna('no',inplace=True)
    xtaldf.replace(to_replace={'enable_fuzzy':{'yes':1,'no':0},'flagged':{'yes':1,'no':0}},inplace=True)
    xtaldf['enable_fuzzy']=xtaldf['enable_fuzzy'].astype(int)
    xtaldf['flagged']=xtaldf['flagged'].astype(int)
    xtalfname=os.path.basename(xtalxlfpath)
    xtalfile_mtime=os.stat(xtalxlfpath).st_mtime
    dbpath=Path(dbpathstr)
    conn=seqdbutils.gracefuldbopen(dbpath)
    c=conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS CURATEDXTALS (acc text, xtalfname text, xtalf_mtime int, \
                cazypdblist text, pdbid text, pdbchain text, ntccpos int, ctccpos int, ntccseq text, ctccseq text,\
                enable_fuzzy int, pgbsr_ccstart int, pgbsr_ccstop int, pgbsr_seqchecksum text)''')
    for _,xrow in xtaldf.dropna(subset=['curated_pdb']).iterrows():
        assert(type(xrow['cc_startnum']) in (float,int)),"cc_startnum type not numeric"
        entry_status=check_cxstatus(c,xrow,xtalfile_mtime)
        if entry_status=='skip':
            continue
        if entry_status in ['update_required','delete']: #drop 
            c.execute('''DELETE FROM CURATEDXTALS WHERE pdbid=(?)''',(xrow['curated_pdb'],))
        if entry_status in ['absent','update_required']: #adding or replacing
            new_tuple=(xrow['gb_accession'],xtalfname,xtalfile_mtime,xrow['PDBS'],xrow['curated_pdb'],\
                 xrow['chain'],int(xrow['cc_startnum']),int(xrow['cc_stopnum']),xrow['cc_startseq'],xrow['cc_stopseq'],\
                 xrow['enable_fuzzy'],None,None,None)
            c.execute('''INSERT INTO CURATEDXTALS VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',new_tuple)
#        else:
#            assert(entry_status=='exists'),'existing entry_status unclear' #otherwise skip since ok already
    #now let's match up start/stops with sequences in database-
    xtalsdbpath=None
    if xtalsdb_relpath=='up1':
        xtalsdbpath= dbpath.resolve().parent.parent / xtalsdb_name
    align_proteingbs_seq(conn,xtalsdbpath)
    conn.commit()
    conn.close()


def mergein_predictedccdata(dbpathstr):
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    seqdbutils.check_tables_exist(conn,['CURATEDXTALS','CCDATA'])
    c=conn.cursor()

    #c.execute('''SELECT CCDATA.acc,CCDATA.ccb,CCDATA.ccstart,CCDATA.ccstop,CURATEDXTALS.ntccpos,CURATEDXTALS.ctccpos FROM CCDATA INNER JOIN CURATEDXTALS ON CCDATA.acc = CURATEDXTALS.acc''')
    c.execute('''SELECT * FROM CCDATA INNER JOIN CURATEDXTALS ON CCDATA.acc = CURATEDXTALS.acc''')
    cxccs=c.fetchall()
    for cxcc in cxccs:
#        curacc=cxcc['acc']
#        c.execute('''SELECT pklgbsr FROM PROTEINGBS WHERE acc=(?)''',(curacc,))
#        gbsr=pickle.loads(c.fetchone())
#        myseq=gbsr.seq
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