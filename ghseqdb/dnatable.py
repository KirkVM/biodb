import re,os,datetime,pickle,time
from Bio import Entrez,SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils.CheckSum import seguid
#from Bio.SeqUtils.CheckSum import seguid
#from . import seqdbutils
from requesters import entrez_requests
from ghseqdb import seqdbutils

def getdnaseqs(acclist,email,api_key):
#    Entrez.email=email
#    Entrez.api_key=api_key
    srs=entrez_requests.getdnasrs(acclist,email,api_key)
    for sr in srs:
        print(sr,len(sr.seq),sr.id)

def update_entry(c,srid,sr):
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month:02d}-{today.day:02d}'
    try:
        accvrsn=sr.annotations['sequence_version']
        accvrsn=int(accvrsn)
    except:
        print(f'cannot get version # for {srid}')
        accvrsn=None
    seq_checksum=seguid(sr.seq)
    pklseq=pickle.dumps(sr.seq)
    update_tuple=(accvrsn,todaystr,todaystr,seq_checksum,pklseq,0,srid)
    c.execute('''UPDATE DNATABLE SET version = (?), ncbiscan_date = (?), modify_date = (?), seq_checksum = (?),\
                pklseq = (?), failcount = (?) WHERE acc = (?)''',update_tuple)

def update_failure(c,srid):
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month:02d}-{today.day:02d}'
    c.execute('''SELECT * FROM DNATABLE WHERE acc=(?)''',(srid,))
    fcount=c.fetchone()['failcount']+1
    update_tuple=(todaystr,fcount,srid)
    c.execute('''UPDATE DNATABLE SET ncbiscan_date = (?), failcount = (?) WHERE acc = (?)''',update_tuple)

def update_if_needed(c,srid,sr):
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month:02d}-{today.day:02d}'
    if sr is None:
        update_failure(c,srid)
    c.execute('''SELECT * FROM DNATABLE WHERE acc=(?)''',(srid,))
    row=c.fetchone()
    if sr is not None:
        newseq_checksum=seguid(sr.seq)
        prevseq_checksum=row['seq_checksum']
        if newseq_checksum==prevseq_checksum:
            c.execute('''UPDATE DNATABLE SET ncbiscan_date = (?) WHERE acc = (?)''',(todaystr,srid))
        else:
            print(f'overwriting {srid} due to diff in checksum')
            update_entry(c,srid,sr)



def build_dnatable(dbpathstr,email,api_key,newdb=False,customlistfile=None,refresh=False,retry_fails=False,stopat=None):
    Entrez.email=email
    Entrez.api_key=api_key
    conn=seqdbutils.gracefuldbopen(dbpathstr) 
    c=conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS DNATABLE (acc text, version text, ncbiscan_date text, modify_date text,\
                seq_checksum text,pklseq glob, failcount int)''')
    assert(os.path.exists(customlistfile))
    if customlistfile:
        with open(customlistfile,'r') as f:
            accs2find=[x.strip() for x in f.readlines()]
    c.execute('''SELECT * FROM DNATABLE''')
    rows=c.fetchall()
    existing_accs=[x['acc'] for x in rows]
    newacc2find=list(set(accs2find).difference(existing_accs))
    print(f'{len(newacc2find)} new accession codes to check')
    print(newacc2find)
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month:02d}-{today.day:02d}'
    for srcnt,srid in enumerate(newacc2find):
        print(srid)
        new_tuple=(srid,None,todaystr,None,None,None,None)
        c.execute('''INSERT INTO DNATABLE VALUES (?,?,?,?,?,?,?)''',new_tuple)
        sr=entrez_requests.getdnasr(srid)
        if sr is not None:
            update_entry(c,srid,sr)
        else:
            update_failure(c,srid)
        if srcnt>0 and srcnt%10==0:
            time.sleep(5)
            conn.commit()
    #now look at old ones
    old2check=list(set(accs2find).union(existing_accs))
    print(f'now checking a bunch of previous: {len(old2check)} many')
    for srid in old2check:
        c.execute('''SELECT * FROM DNATABLE WHERE acc=(?)''',(srid,))
        row=c.fetchone()
        existingdate=datetime.date(*[int(x) for x in row['ncbiscan_date'].split('-')])
        days_since_update=(today-existingdate).days
        if days_since_update>30:
            sr=entrez_requests.getdnasr(srid)
            update_if_needed(c,srid,sr)
    conn.commit()