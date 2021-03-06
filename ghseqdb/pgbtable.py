import time,datetime,pickle
from Bio import Entrez
from Bio.SeqUtils.CheckSum import seguid
from . import seqdbutils
from requesters import entrez_requests

def get_taxid(sr):
    try:
        for sf in sr.features:
            if sf.type=='source':
                for k in sf.qualifiers:
                    if k=='db_xref':
                        dbxrs=sf.qualifiers[k]
                        for dbxr in dbxrs:
                            svals=dbxr.split(':')
                            if svals[0]=='taxon':
                                taxid=int(svals[1])
        return taxid
    except:
        return None
 
def updatedb(c,acc,failcount):
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month:02d}-{today.day:02d}'
    sr=entrez_requests.getgbpsr(acc)
    goodsr=True
    if sr is None: 
        print(f'could not download {acc}')
        goodsr=False
    elif sr.name!=acc:
        print(f"sr.name {sr.name} does not match acc {acc}. Checking sr.id...")
        if acc in sr.id:
            print(f"forgiven... acc= {acc}, sr.id= {sr.id}, seqlen= {len(sr.seq)}.")
        else:
            print(f"no good. sr.id= {sr.id}, seqlen= {len(sr.seq)}. skipping...")
            goodsr=False
    if goodsr:
        try:
            accvrsn=sr.annotations['sequence_version']
            accvrsn=int(accvrsn)
        except:
            print(f'cannot get version # for {acc}')
            accvrsn=None
        seq_checksum=seguid(sr.seq)
        sr_checksum=seguid(sr)
        if failcount==0:
            new_tuple=(acc,accvrsn,get_taxid(sr),todaystr,todaystr,seq_checksum,sr_checksum,pickle.dumps(sr),0)
            c.execute('''INSERT INTO PROTEINGBS VALUES (?,?,?,?,?,?,?,?,?)''',new_tuple)
        else:
            update_tuple=(accvrsn,get_taxid(sr),todaystr,todaystr,seq_checksum,sr_checksum,pickle.dumps(sr),0,acc)
            c.execute('''UPDATE PROTEINGBS SET version = (?), taxid=(?), ncbiscan_date = (?), modify_date = (?), seq_checksum = (?), \
                            sr_checksum = (?), pklgbsr = (?), failcount = (?) WHERE acc = (?)''',update_tuple)
    else:
        if failcount==0:
            new_tuple=(acc,None,None,todaystr,None,None,None,None,1)
            c.execute('''INSERT INTO PROTEINGBS VALUES (?,?,?,?,?,?,?,?,?)''',new_tuple)
        else:
            update_tuple=(failcount+1,todaystr,acc)
            c.execute('''UPDATE PROTEINGBS SET failcount = (?), ncbiscan_date = (?) WHERE acc = (?)''',update_tuple)

def get_tblaccs(cursor):
    accs=[]
    try:
        cursor.execute('''SELECT acc FROM CAZYSEQDATA''')
        accs2add=[x['acc'] for x in cursor.fetchall()]
        accs.extend(accs2add)
    except:
        print('no CAZYSEQDATA table')
    try:
        cursor.execute('''SELECT acc FROM FROZENSEQS''')
        accs.extend([x['acc'] for x in cursor.fetchall()])
    except:
        print('no FROZENSEQS table')
    return accs

def build_proteingbtable(dbpath,email,api_key,refresh=False,retry_fails=False,stopat=None):
    Entrez.email=email
    Entrez.api_key=api_key
    conn=seqdbutils.gracefuldbopen(dbpath) 
    c=conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS PROTEINGBS (acc text, version text, taxid int, ncbiscan_date text,\
                 modify_date text, seq_checksum text, sr_checksum, pklgbsr glob, failcount int)''')
    c.execute('''SELECT acc FROM PROTEINGBS''')
    cur_pgbaccs=[x['acc'] for x in c.fetchall()]
    accs2find=get_tblaccs(c)  #replace this with seqdbutils.check_table_exists and then a select
    print(f'{len(cur_pgbaccs)} existing (successful) entries in PROTEINGBS table')
    dlaccs=list(set(accs2find).difference(cur_pgbaccs))

    print(f'{len(dlaccs)} remaining seqs to pull')
    for snum,dlacc in enumerate(dlaccs):
        updatedb(c,dlacc,0)
        if snum>0 and snum%25==0:
            print(f'through {snum} sequences')
            conn.commit()
            time.sleep(5)
        if stopat is not None and snum==stopat:
            break
    if retry_fails:
        c.execute('''SELECT * FROM PROTEINGBS WHERE failcount!=0''')
        rtrows=c.fetchall()
        for snum,rtrow in enumerate(rtrows):
            updatedb(c,rtrow['acc'],rtrow['failcount'])
            if snum>0 and snum%25==0:
                print(f'through {snum} retry sequences')
                conn.commit()
                time.sleep(5)

    conn.commit()
    conn.close()


def extractsrs(dbpathstr,format='fasta',resetid_db_acc=False):
    """returns a file containing all the sequence files in PROTEINGBS table"""
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    seqdbutils.check_tables_exist(conn,["PROTEINGBS"])
    c=conn.cursor()
    c.execute('''SELECT * FROM PROTEINGBS WHERE acc NOT NULL''')
    allrows=c.fetchall()
    print(f'found {len(allrows)} non-null entries')
    missingsr_accs=[]
    srs=[]
    for row in allrows:
        pklgbsr=row['pklgbsr']
        if pklgbsr is not None:
            newsr=pickle.loads(row['pklgbsr'])
            if resetid_db_acc:
                newsr.id=row['acc']
            srs.append(newsr)
        else:
            missingsr_accs.append(row['acc'])
    conn.close()
    return srs



def update_table_with_taxid(dbpathstr):
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    seqdbutils.check_tables_exist(conn,["PROTEINGBS"])
    c=conn.cursor()
    #see if we need to add taxid row
    c.execute('''SELECT * FROM PROTEINGBS''')
    row=c.fetchone()
    if 'taxid' not in row.keys():
        print('adding column taxid')
        c.execute('''ALTER TABLE PROTEINGBS ADD COLUMN taxid int''')
    
    #now go through the list
    c.execute('''SELECT acc,pklgbsr FROM PROTEINGBS WHERE pklgbsr NOT NULL AND taxid IS NULL''')
    rows=c.fetchall()
    print(f'examining {len(rows)} rows for taxids')
    for row in rows:
        sr=pickle.loads(row['pklgbsr'])
        taxid=get_taxid(sr)
        if taxid is not None:
            update_tuple=(taxid,row['acc'])
            c.execute('''UPDATE PROTEINGBS SET taxid = (?) WHERE acc = (?)''',update_tuple)
    conn.commit()
    conn.close()
