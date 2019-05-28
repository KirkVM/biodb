import sqlite3,re,datetime
from . import seqdbutils
from scrapers import cazydbscrapers

def build_cazytable(ghfam,dbpathstr,drop_old=False):
    """scrapes CAZY db for accession codes/annotations info, then downloads seqs through NCBI Entrez

    Arguments:
        ghfam: shorthand name of GH family of interest (GH5, GH43, etc)
        email: email to use in registering with Entrez eutil API
        outfolder: path to folder to output sequence files (creates by default if needed)

    Returns:
        sqlite db (also writes out pseq fasta file)
    """    
    conn=seqdbutils.gracefuldbopen(dbpathstr,create_new=True) 
    c=conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS CAZYSEQDATA (acc text, version text, scrapedate text, \
                subfam text, extragbs text, ecs text, pdbids text, uniprotids text)''')
    print('starting CAZY scrape')
    czes_=cazydbscrapers.scrape_cazyfam(f'{ghfam}')
    assert(len(czes_)>0), f"Unable to scrape CAZY for selection {ghfam}"
    print(f'found {len(czes_)} entries. building DB')
    add_count=0
    update_count=0
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month:02d}-{today.day:02d}'
    accRE=re.compile("(.+)\.(\d+)")
    czegbs=[]
    for cze in czes_:
        maingbacc=cze.gbids_[0]
        try:
            acc,accvrsn=accRE.match(maingbacc).groups()
        except:
            print(f'no sequence version for {maingbacc}')
            acc=maingbacc
            accvrsn=None
        czegbs.append(acc)
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
            prev_entry=existingentries[0]
            existingdate=datetime.date(*[int(x) for x in prev_entry['scrapedate'].split('-')])
            days_since_update=(today-existingdate).days
            if days_since_update>15: #compare values
                cazy_changes=False
                if subfam!=prev_entry['subfam'] or extragbs!=prev_entry['extragbs'] or \
                               ecs!=prev_entry['ecs'] or pdbids!=prev_entry['pdbids']:
                    print(f'change to {acc}. Updating')
                    cazy_changes=True
                if subfam!=prev_entry['subfam']:
                    print(f'++++ subfam change from {prev_entry["subfam"]} to {subfam}')
                if extragbs!=prev_entry['extragbs']:
                    print(f'++++ extragbs change from {prev_entry["extragbs"]} to {extragbs}')
                if ecs!=prev_entry['ecs']:
                    print(f'+___ECECECECECECEC___ ecs change from {prev_entry["ecs"]} to {ecs}___ECECECECECECEC___+')
                if pdbids!=prev_entry['pdbids']:
                    print(f'+___PDBPDBPDBPDBPDB____ pdbids change from {prev_entry["pdbids"]} to {pdbids} ___PDBPDBPDBPDBPDB____+')
                if accvrsn!=prev_entry['version']:
                    print(f'++++ accvrsn change for {acc}** not updating if this is the only cHange **')
                if cazy_changes:
                    update_tuple=(accvrsn,todaystr,subfam,extragbs,ecs,pdbids,uniprotids,acc)
                    c.execute('''UPDATE CAZYSEQDATA SET version = (?), scrapedate = (?), subfam = (?), \
                            extragbs = (?), ecs = (?), pdbids = (?), uniprotids = (?) WHERE acc = (?)''',update_tuple)
                    update_count+=1
        ###still need logic to drop old entries, if no longer there or re-classified as an alt id of a difft entry
    conn.commit()
    c.execute('''SELECT * FROM CAZYSEQDATA''')
    dbrows=c.fetchall()
    dbgbs=[x['acc'] for x in dbrows]
    missing_gbs=list(set(dbgbs).difference(czegbs))
    missing_count=len(missing_gbs)
    if missing_count>0:
        print(f'{missing_count} CAZYSEQDATA rows could be dropped based on cazy.org. Consider setting drop_old=True')
    for mgb in missing_gbs:
        print(f'{mgb} in CAZYSEQDATA but no longer on cazy.org for this family')
        #TODO iterate through CAZYSEQDDATA WHERE extragbs is not null and compare
    conn.close()
    print(f'added {add_count} entries, updated {update_count} entries')
