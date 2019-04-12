import re,os,datetime
from Bio.PDB import PDBList
from pathlib import Path
from . import seqdbutils

def build_xtalstable(dbpathstr,sourcedbpathstrs,pdbformat='cif'):
    dbpath=Path(dbpathstr)
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month:02d}-{today.day:02d}'
    conn=seqdbutils.gracefuldbopen(dbpath) 
    xtaldirpath=Path(dbpath).parent / 'Xtals'
    if not xtaldirpath.exists():
        os.mkdir(xtaldirpath)
    c=conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS XTALS (pdbid text, acc text, srcdb text,dldate text, 
              relpath text, pdbformat text, dlsuccess int, obsolete int)''')
    pdbl=PDBList()
    cxRE=re.compile('([A-Za-z0-9]{4})\[')
    #by default folder containing xtals is named 'xtals' and in same directory as xtalsdb
    for srcdbpathstr in sourcedbpathstrs:
        srcdbpath=Path(srcdbpathstr)
        srcdbstr=srcdbpath.name
        src_conn=seqdbutils.gracefuldbopen(srcdbpath)    
        seqdbutils.check_tables_exist(src_conn,['CAZYSEQDATA'])
        src_c=src_conn.cursor()
        src_c.execute('SELECT acc,pdbids FROM CAZYSEQDATA WHERE pdbids NOT NULL')
        pdbrows=src_c.fetchall()
        dbpdbs=[]
        dbaccs=[]
        for pdbrow in pdbrows:
            pdbentry=pdbrow['pdbids']
            accentry=pdbrow['acc']
            pdbs=cxRE.findall(pdbentry)
            dbaccs.extend([accentry for _ in range(len(pdbs))]) 
            dbpdbs.extend(pdbs)
        #dbpdbs=['4IM4']
        pdbl.download_pdb_files(dbpdbs,pdir=xtaldirpath)#,obsolete=True)
        for acc,pdb in zip(dbaccs,dbpdbs):
            rel_pdbpath=xtaldirpath / f'{pdb}.cif'
            download_success=os.path.exists(rel_pdbpath)
            str_relpath=str(rel_pdbpath)
            c.execute('''SELECT COUNT(*) FROM XTALS WHERE pdbid=(?) AND dlsuccess=(?)''',(pdb,1))
            already_downloaded=c.fetchone()[0]
            if already_downloaded:
                continue
            c.execute('''SELECT COUNT(*) FROM XTALS WHERE pdbid=(?) AND dlsuccess=(?)''',(pdb,0))
            previously_failed=c.fetchone()[0]
            if previously_failed:
                if download_success:
                    print(f'new download of previously failed {pdb}')
                    c.execute('''UPDATE XTALS SET dldate = (?), dlsuccess = (?) WHERE pdbid=(?)''',(todaystr,1,pdb))
                continue
            c.execute('''INSERT INTO XTALS VALUES (?,?,?,?,?,?,?,?)''',\
                      (pdb,acc,srcdbstr,todaystr,str_relpath,pdbformat,download_success,None))
        conn.commit()
        src_conn.close()
    conn.close()

def add_obsoletevalues(dbpathstr):
    dbpath=Path(dbpathstr)
    conn=seqdbutils.gracefuldbopen(dbpath) 
    c=conn.cursor()
    pdbl=PDBList()
    obsolete_list=pdbl.get_all_obsolete()
    for ocode in obsolete_list:
        c.execute('''SELECT * FROM XTALS WHERE pdbid=(?)''',(ocode,))
        obsrow=c.fetchone()
        if obsrow is not None:
            srcdb=obsrow['srcdb']
            srcacc=obsrow['acc']
            print(f'PDB {ocode} (acc-{srcacc}) from {srcdb} marked obsolete')
            c.execute('''UPDATE XTALS SET obsolete=(?) WHERE pdbid=(?)''',(1,ocode))
    conn.commit()
    conn.close()


#how to get sequence from a structure-
#from Bio.PDB import MMCIFParser,PPBuilder
#parser=MMCIFParser()
#structure=parser.get_structure('coolname',pdbloc)
#chain=structure[0]['A']
#ppb=PPBuilder()
#ppchain=ppb(chain)
#seq=ppchain[0].get_sequence()
