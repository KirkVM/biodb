import time,datetime,pickle
from Bio import Entrez
from Bio.SeqUtils.CheckSum import seguid
from . import seqdbutils
from requesters import entrez_requests

def updatedb(c,acc,failcount):
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month:02d}-{today.day:02d}'
    sr=entrez_requests.getgbpsr(acc)
    if sr is not None:
        assert(sr.name==acc),"why doesn't sr.name==dlacc?"
        try:
            accvrsn=sr.annotations['sequence_version']
            accvrsn=int(accvrsn)
        except:
            print(f'cannot get version # for {acc}')
            accvrsn=None
        seq_checksum=seguid(sr.seq)
        sr_checksum=seguid(sr)
        if failcount==0:
            new_tuple=(acc,accvrsn,todaystr,todaystr,seq_checksum,sr_checksum,pickle.dumps(sr),0)
            c.execute('''INSERT INTO PROTEINGBS VALUES (?,?,?,?,?,?,?,?)''',new_tuple)
        else:
            update_tuple=(accvrsn,todaystr,todaystr,seq_checksum,sr_checksum,pickle.dumps(sr),0,acc)
            c.execute('''UPDATE PROTEINGBS SET version = (?), ncbiscan_date = (?), modify_date = (?), seq_checksum = (?), \
                            sr_checksum = (?), pklgbsr = (?), failcount = (?) WHERE acc = (?)''',update_tuple)
    else:
        if failcount==0:
            new_tuple=(acc,None,todaystr,None,None,None,None,1)
            c.execute('''INSERT INTO PROTEINGBS VALUES (?,?,?,?,?,?,?,?)''',new_tuple)
        else:
            update_tuple=(failcount+1,todaystr,acc)
            c.execute('''UPDATE PROTEINGBS SET failcount = (?), ncbiscan_date = (?) WHERE acc = (?)''',update_tuple)
        print(f'could not download {acc}')

def get_tblaccs(cursor):
    accs=[]
    try:
        cursor.execute('''SELECT acc FROM CAZYSEQDATA''')
        accs2add=[x['acc'] for x in cursor.fetchall()]
        accs.extend(accs2add)
    except:
        print('no CAZYSEQDATA table')
    return accs

def build_proteingbtable(dbpath,email,api_key,refresh=False,retry_fails=False,stopat=None):
    Entrez.email=email
    Entrez.api_key=api_key
    conn=seqdbutils.gracefuldbopen(dbpath) 
    c=conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS PROTEINGBS (acc text, version text, ncbiscan_date text,\
                 modify_date text, seq_checksum text, sr_checksum, pklgbsr glob, failcount int)''')
    c.execute('''SELECT acc FROM PROTEINGBS WHERE failcount=0''')
    cur_pgbaccs=[x['acc'] for x in c.fetchall()]
    accs2find=get_tblaccs(c) 
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


def extractsrs(dbpath,format='fasta'):
    """returns a file containing all the sequence files in PROTEINGBS table"""
    conn=seqdbutils.gracefuldbopen(dbpath)
    c=conn.cursor()
    c.execute('''SELECT * FROM PROTEINGBS WHERE acc NOT NULL''')
    allrows=c.fetchall()
    print(f'found {len(allrows)} non-null entries')
    missingsr_accs=[]
    srs=[]
    for row in allrows:
        pklgbsr=row['pklgbsr']
        if pklgbsr is not None:
            srs.append(pickle.loads(row['pklgbsr']))
        else:
            missingsr_accs.append(row['acc'])
    conn.close()
    return srs