import sqlite3,re,datetime
from . import seqdbutils
from scrapers import cazydbscrapers

def build_cazytable(ghfam,dbpath):
    """scrapes CAZY db for accession codes/annotations info, then downloads seqs through NCBI Entrez

    Arguments:
        ghfam: shorthand name of GH family of interest (GH5, GH43, etc)
        email: email to use in registering with Entrez eutil API
        outfolder: path to folder to output sequence files (creates by default if needed)

    Returns:
        sqlite db (also writes out pseq fasta file)
    """    
    conn=seqdbutils.gracefuldbopen(dbpath) 
    c=conn.cursor()
#    try:
#        c.execute('''SELECT COUNT (*) FROM CAZYSEQDATA''')
#        cazydbsize=c.fetchone()[0]

    c.execute('''CREATE TABLE IF NOT EXISTS CAZYSEQDATA (acc text, version text, scrapedate text, \
                subfam text, extragbs text, ecs text, pdbids text, uniprotids text)''')
    print('starting CAZY scrape')
    czes_=cazydbscrapers.scrape_cazyfam(f'{ghfam}')
    assert(len(czes_)>0), f"Unable to scrape CAZY for selection {ghfam}"
    print(f'found {len(czes_)} entries. building DB')
    add_count=0
    update_count=0
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month}-{today.day}'
    accRE=re.compile("(.+)\.(\d+)")
    for cze in czes_:
        maingbacc=cze.gbids_[0]
        try:
            acc,accvrsn=accRE.match(maingbacc).groups()
        except:
            print(f'no sequence version for {maingbacc}')
            acc=maingbacc
            accvrsn=None
        c.execute('''SELECT * FROM CAZYSEQDATA WHERE acc = (?)''',(acc,))
        existingentries=c.fetchall()
        assert(len(existingentries))<=1, f"more than 1 entry exists for {acc}"
        subfam=None
        extragbs=None
        ecs=None
        pdbids=None
        uniprotids=None
        if cze.family!=None:
            subfam=cze.family
        if len(cze.gbids_)>1:
            extragbs=''
            for egb in cze.gbids_[1:]:
                extragbs+=f'{egb},'
            extragbs=extragbs[:-1]
        if len(cze.ecs_)>0:
            ecs=''
            for ec in cze.ecs_:
                ecs+=f'{ec},'
            ecs=ecs[:-1]
        if len(cze.pdbids_)>0:
            pdbids=''
            for pdbid in cze.pdbids_:
                pdbids+=f'{pdbid},'
            pdbids=pdbids[:-1]
        if len(cze.uniprotids_)>0:
            uniprotids=''
            for uniprotid in cze.uniprotids_:
                uniprotids+=f'{uniprotid},'
            uniprotids=uniprotids[:-1]

        #now update entry if it's been over 1 month
        if len(existingentries)==0:
            new_tuple=(acc,accvrsn,todaystr,subfam,extragbs,ecs,pdbids,uniprotids)
            c.execute('''INSERT INTO CAZYSEQDATA VALUES (?,?,?,?,?,?,?,?)''',new_tuple)
            add_count+=1
        else:
            update_tuple=(accvrsn,todaystr,subfam,extragbs,ecs,pdbids,uniprotids,acc)
            existingdate=datetime.date(*[int(x) for x in existingentries[0]['scrapedate'].split('-')])
            days_since_update=(today-existingdate).days
            if days_since_update>30:
                c.execute('''UPDATE CAZYSEQDATA SET version = (?), scrapedate = (?), subfam = (?), \
                            extragbs = (?), ecs = (?), pdbids = (?), uniprotids = (?) WHERE acc = (?)''',update_tuple)
                update_count+=1
    conn.commit()
    conn.close()
    print(f'added {add_count} entries, updated {update_count} entries')
