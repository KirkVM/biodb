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
            prevrow['ntfullpos']==xrow['full_startnum'] and prevrow['ctfullpos']==xrow['full_stopnum'] and
            prevrow['ntfullseq']==xrow['full_startseq'] and prevrow['ctfullseq']==xrow['full_stopseq'] and
            prevrow['acc']==xrow['gb_accession'] and prevrow['cazypdblist']==xrow['PDBS'] and
            prevrow['enable_fuzzy']==xrow['enable_fuzzy']):
            return 'exists'
        else:
            assert (xtalfile_mtime>prevrow['xtalf_mtime']),"replacing with val from older excel file?"
            return 'update_required'
    return 'prev_status_unclear'

class AlnRegister:
    def __init__(self,enzyme_format,acc,pdbid):
        self.enzyme_format=enzyme_format
        self.acc=acc
        self.pdbid=pdbid
        self.start=None
        self.stop=None
        self.message=''
        self.success=False

def do_register_aln(pgbseq,curatedrow,eformat,alnreg,mode,xtalseqstr=None):
    if eformat=='cc':
        curated_startseq=curatedrow['ntccseq']
        curated_stopseq=curatedrow['ctccseq']
        curated_startresnum=curatedrow['ntccpos']
        curated_stopresnum=curatedrow['ctccpos']
    elif eformat=='full':
        curated_startseq=curatedrow['ntfullseq']
        curated_stopseq=curatedrow['ctfullseq']
        curated_startresnum=curatedrow['ntfullpos']
        curated_stopresnum=curatedrow['ctfullpos']
    
    inter_seqlen=curated_stopresnum-len(curated_stopseq) - (curated_startresnum+len(curated_startseq)) + 1
    inter_max=inter_seqlen+2
    inter_min=inter_seqlen-2
    seqrestring=f'({curated_startseq}.{{{inter_min},{inter_max}}}{curated_stopseq})'
    seqRE=re.compile(seqrestring) 
    if mode=='easy': 
        pgbmatch=seqRE.search(str(pgbseq))
        if pgbmatch:
            alnreg.start,alnreg.stop=pgbmatch.span()
            alnreg.stop-=1
            alnreg.message+=f'found {curatedrow["pdbid"]} {alnreg.enzyme_format} the easy way by matching genbank sequence'
            alnreg.message+=f' {curatedrow["acc"]}:({alnreg.start},{alnreg.stop})\n'
            alnreg.success=True
        else:
            alnreg.message+=f'-->match to genbank seq failed for {curatedrow["pdbid"]} for eformat {alnreg.enzyme_format}\n'
    elif mode=='xtal':
        xstrmatch=seqRE.search(xtalseqstr)
        if xstrmatch:
            max_ratio=0
            xtalcoreseqstr=xstrmatch.group(1)
            for strtcand in range(0,len(pgbseq)-inter_min ):
                for stopcand in range(strtcand+inter_min,strtcand+inter_max):
                    sm=difflib.SequenceMatcher(a=xtalcoreseqstr,b=str(pgbseq),autojunk=False)
                    if sm.ratio()>max_ratio:
                        alnreg.start=strtcand
                        alnreg.stop=stopcand-1
                        alnreg.success=True #may continually overwrite best register
        if alnreg.success:
            alnreg.message+= f'-->found {curatedrow["pdbid"]} by matching xtalpeptide to genbank sequence'
            alnreg.message+= f' {curatedrow["acc"]}: ({alnreg.start},{alnreg.stop})\n'
        else:
            alnreg.message+=f'---->match to xtalseq failed for {curatedrow["pdbid"]} for eformat {alnreg.enzyme_format}\n'
    elif mode=='fuzzy':
        substarts=[]
        substops=[]
        for strtcand in range(0,len(pgbseq)-inter_min+1 ):
            sm_start=difflib.SequenceMatcher(a=curated_startseq,b=str(pgbseq)[strtcand:strtcand+5])
            if sm_start.ratio()>0.59:
                substarts.append(strtcand)
        for stopcand in range(inter_min-1,len(pgbseq)):
            sm_stop=difflib.SequenceMatcher(a=curated_stopseq,b=str(pgbseq)[stopcand-4:stopcand+1])
            if sm_stop.ratio()>0.59:
                substops.append(stopcand)
        foundit=False
        for substart in substarts:
            if alnreg.success: 
                break
            for substop in substops:
                if foundit: break
                if ((substop-substart)<(curated_stopresnum-curated_startresnum+5) and\
                    (substop-substart)>(curated_stopresnum-curated_startresnum-5)):
                    for xpos in range(0,len(pgbseq)-len(xtalseqstr)+5):
                        sm=difflib.SequenceMatcher(a=xtalseqstr,b=str(pgbseq)[xpos:xpos+len(xtalseqstr)])
                        if sm.ratio()>0.95:
                            foundit=True
                            alnreg.start=substart
                            alnreg.stop=substop-1
                            alnreg.success=True
                        if alnreg.success:
                            break
        if alnreg.success:
            alnreg.message+=f'---->Fuzzy match worked for xtal {curatedrow["pdbid"]} - {curatedrow["acc"]},\
                            eformat {alnreg.enzyme_format}: ({alnreg.start},{alnreg.stop})'
        else:
            alnreg.message+=f'|||||||||FAILED to find register for {curatedrow["pdbid"]} despite enabling fuzzy matching|||||||||||'
    #do a final clean-up to round to start or stop of sequence--
    if alnreg.success:
        if alnreg.start<3:
            alnreg.start=0
        if alnreg.stop>len(pgbseq)-3:
            alnreg.stop=len(pgbseq)-1
    return alnreg
 


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

    aln_registers=[] 
    c.execute('''SELECT * FROM CURATEDXTALS''')
    curatedrows=c.fetchall()
    for curatedrow in curatedrows:

        c.execute('''SELECT * FROM PROTEINGBS WHERE acc=(?)''',(curatedrow['acc'],))
        pgbrow=c.fetchone()
        pgbseq=pickle.loads(pgbrow['pklgbsr']).seq

        if (pgbrow['seq_checksum']==curatedrow['pgbsr_seqchecksum'] and \
            curatedrow['pgbsr_ccstart'] is not None and curatedrow['pgbsr_ccstop'] is not None and
            curatedrow['pgbsr_fullstart'] is not None and curatedrow['pgbsr_fullstop'] is not None):
            print(f'skipping previously-completed {curatedrow["pdbid"]}')
            continue #skip this one!

        print(f'now examining {curatedrow["pdbid"]}')
        for enzyme_format in ['cc','full']:   
            aln_register=AlnRegister(enzyme_format,curatedrow['acc'],curatedrow['pdbid']) 
            #just copy over cc values to full if possible---
            if (enzyme_format=='full' and curatedrow['ntccpos']==curatedrow['ctccpos'] and \
                curatedrow['ntfullpos']==curatedrow['ctfullpos']): 
                aln_register.start=aln_registers[-1].start
                aln_register.stop=aln_registers[-1].stop
                aln_register.message='used cc values for full'
                aln_register.success=aln_registers[-1].success
                continue
            
            else:
                aln_register=do_register_aln(pgbseq,curatedrow,enzyme_format,aln_register,'easy')
                if not aln_register.success:
                    #get sequence from pdb file
                    xc.execute('''SELECT * FROM XTALS WHERE pdbid=(?)''',(curatedrow['pdbid'],))
                    xtalrow=xc.fetchone()
                    xtalfpath= xtalsdbpath.parent / xtalrow['relpath']
                    xtalstruct=mmcparser.get_structure(curatedrow['pdbid'],xtalfpath)
                    xtalchain=xtalstruct
                    xtalppchain=ppb.build_peptides(xtalchain)
                    xtalseqstr=str(xtalppchain[0].get_sequence())
                    aln_register=do_register_aln(pgbseq,curatedrow,enzyme_format,aln_register,'xtal',xtalseqstr=xtalseqstr)
                    if not aln_register.success:
                        if curatedrow['enable_fuzzy']:
                            aln_register=do_register_aln(pgbseq,curatedrow,enzyme_format,aln_register,'fuzzy',xtalseqstr=xtalseqstr)
                        else:
                            aln_register.message='************ * * * *  * * *  *  *   *   * * *    *   *   *  **  * *    *********************\n'
                            aln_register.message+=f'-------could not match. Consider setting enable_fuzzy column to yes for {curatedrow["pdbid"]}-------\n'
                            aln_register.message+='************  *  * *  *  * *  *  * *    *   *   * *   * *  * *   *  * *  *******************\n'
                aln_registers.append(aln_register)

    for aln_register in aln_registers:
        print(aln_register.message)
        if aln_register.success:
            update_tuple=(aln_register.start,aln_register.stop,pgbrow['seq_checksum'],aln_register.acc)
            if aln_register.enzyme_format=='cc':
                c.execute('''UPDATE CURATEDXTALS SET pgbsr_ccstart=(?), pgbsr_ccstop=(?), pgbsr_seqchecksum=(?) \
                            WHERE acc=(?)''', update_tuple)
            if aln_register.enzyme_format=='full':
                c.execute('''UPDATE CURATEDXTALS SET pgbsr_fullstart=(?), pgbsr_fullstop=(?), pgbsr_seqchecksum=(?) \
                            WHERE acc=(?)''', update_tuple)
        else:
            print(f'+++++------FAILED TO FIND {aln_register.enzyme_format} START,STOP FOR {aln_register.pdbid} ({aln_register.acc})------++++++\n')
    xdbconn.close()




def build_curatedxtaltable(dbpathstr,xtalxlfpath,sheet_name="Sheet1",xtalsdb_relpath='up1',xtalsdb_name='EDATA.sql'):
    xtaldf=pd.read_excel(xtalxlfpath,sheet_name=sheet_name,dtype={'flagged':str})
    xtaldf['enable_fuzzy'].fillna('no',inplace=True)
    xtaldf['flagged'].fillna('no',inplace=True)
    xtaldf.replace(to_replace={'enable_fuzzy':{'yes':1,'no':0},'flagged':{'yes':1,'no':0}},inplace=True)
    xtaldf['enable_fuzzy']=xtaldf['enable_fuzzy'].astype(int)
    xtaldf['flagged']=xtaldf['flagged'].astype(int)
    xtaldf.full_startnum.fillna(xtaldf.cc_startnum,inplace=True)
    xtaldf.full_startseq.fillna(xtaldf.cc_startseq,inplace=True)
    xtaldf.full_stopnum.fillna(xtaldf.cc_stopnum,inplace=True)
    xtaldf.full_stopseq.fillna(xtaldf.cc_stopseq,inplace=True)
    xtalfname=os.path.basename(xtalxlfpath)
    xtalfile_mtime=os.stat(xtalxlfpath).st_mtime
    dbpath=Path(dbpathstr)
#    return xtaldf
    conn=seqdbutils.gracefuldbopen(dbpath)
    c=conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS CURATEDXTALS\
                (acc text, xtalfname text, xtalf_mtime int, cazypdblist text, pdbid text, pdbchain text,\
                ntccpos int, ctccpos int, ntccseq text, ctccseq text,\
                ntfullpos int, ctfullpos int, ntfullseq text, ctfullseq text,\
                pgbsr_ccstart int, pgbsr_ccstop int, pgbsr_fullstart int, pgbsr_fullstop int,
                enable_fuzzy int, pgbsr_seqchecksum text)''')
    for _,xrow in xtaldf.dropna(subset=['curated_pdb']).iterrows():
        assert(type(xrow['cc_startnum']) in (float,int)),"cc_startnum type not numeric"
        entry_status=check_cxstatus(c,xrow,xtalfile_mtime)
        if entry_status=='skip':
            continue
        if entry_status in ['update_required','delete']: #drop 
            c.execute('''DELETE FROM CURATEDXTALS WHERE pdbid=(?)''',(xrow['curated_pdb'],))
        if entry_status in ['absent','update_required']: #adding or replacing
            new_tuple=(xrow['gb_accession'],xtalfname,xtalfile_mtime,xrow['PDBS'],xrow['curated_pdb'],xrow['chain'],\
                int(xrow['cc_startnum']),int(xrow['cc_stopnum']),xrow['cc_startseq'],xrow['cc_stopseq'],\
                int(xrow['full_startnum']),int(xrow['full_stopnum']),xrow['full_startseq'],xrow['full_stopseq'],\
                 None,None,None,None,xrow['enable_fuzzy'],None)
            c.execute('''INSERT INTO CURATEDXTALS VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',new_tuple)
#        else:
#            assert(entry_status=='exists'),'existing entry_status unclear' #otherwise skip since ok already
    #now let's match up start/stops with sequences in database-
    xtalsdbpath=None
    if xtalsdb_relpath=='up1':
        xtalsdbpath= dbpath.resolve().parent.parent / xtalsdb_name
    align_proteingbs_seq(conn,xtalsdbpath)
    conn.commit()
    conn.close()

def extract_curatedsrs(dbpathstr,form='cc',resetid_db_acc='False'):
    """returns a list of type SeqRecord containing all the sequence files in CURATEDXTALS, OTHERS? table"""
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    seqdbutils.check_tables_exist(conn,['CURATEDXTALS','PROTEINGBS'])
    c=conn.cursor()

    c.execute('''SELECT * FROM CURATEDXTALS INNER JOIN PROTEINGBS \
                 ON CURATEDXTALS.acc=PROTEINGBS.acc WHERE CURATEDXTALS.pdbid NOT NULL''') #is WHERE needed?
    allrows=c.fetchall()
    print(f'found {len(allrows)} non-null entries from CURATEDXTALS')

    srs=[]
    for row in allrows:
        pklgbsr=row['pklgbsr']
        if pklgbsr is not None:
            newsr=pickle.loads(row['pklgbsr'])
            newsr.seq=newsr.seq[row['pgbsr_ccstart']:row['pgbsr_ccstop']]
            if resetid_db_acc:
                newsr.id=row['acc']
            srs.append(newsr)
    conn.close()
    return srs

#def eval_ccdata(dbpathstr):
#    conn=seqdbutils.gracefuldbopen(dbpathstr)
#    seqdbutils.check_tables_exist(conn,['CURATEDXTALS','CCDATA'])
#    c=conn.cursor()
#
#    #c.execute('''SELECT CCDATA.acc,CCDATA.ccb,CCDATA.ccstart,CCDATA.ccstop,CURATEDXTALS.ntccpos,CURATEDXTALS.ctccpos FROM CCDATA INNER JOIN CURATEDXTALS ON CCDATA.acc = CURATEDXTALS.acc''')
#    c.execute('''SELECT * FROM CCDATA INNER JOIN CURATEDXTALS ON CCDATA.acc = CURATEDXTALS.acc''')
#    cxccs=c.fetchall()
#    for cxcc in cxccs:
#        try:
#            print(cxcc['acc'],cxcc['ccstart'],cxcc['ccstop'],cxcc['ntccpos'],cxcc['ctccpos'])
#            print(cxcc['ntccseq'],cxcc['ctccseq'])
#        except:
#            conn.close()
#            sys.exit()
#    conn.close()

#tests-numtypes coming out of db
#unique list of accs in PROTEINGBS