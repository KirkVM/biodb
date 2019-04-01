import argparse,os,sqlite3,datetime# sqlite3,os,sys,argparse,yaml
from Bio import SeqIO,Entrez
from . import cazydbscrapers
from . import entrez_requests
import sqlite3,atexit,re
import pickle,time

def updatedb(c,acc,failcount):
    today=datetime.date.today()
    todaystr=f'{today.year}-{today.month}-{today.day}'
    try:
        sr=entrez_requests.getgbpsr(acc)
        assert(sr.name==acc),"why doesn't sr.name==dlacc?"
        accvrsn=sr.annotations['sequence_version']
        try:
            accvrsn=int(accvrsn)
        except:
            print(f'cannot get version # for {acc}')
            accvrsn=None
        checksum=seguid(sr.seq)
        if failcount==0:
            new_tuple=(acc,accvrsn,todaystr,checksum,pickle.dumps(sr),0)
            c.execute('''INSERT INTO PROTEINGBS VALUES (?,?,?,?,?,?)''',new_tuple)
        else:
            update_tuple=(accvrsn,todaystr,checksum,pickle.dumps(sr),0,acc)
            c.execute('''UPDATE PROTEINGBS SET version = (?), dldate = (?), checksum = (?), \
                            pklgbsr = (?), failcount = (?) WHERE acc = (?)''',update_tuple)
    except:
        if failcount==0:
            new_tuple=(acc,None,None,None,None,1)
            c.execute('''INSERT INTO PROTEINGBS VALUES (?,?,?,?,?,?)''',new_tuple)
        else:
            update_tuple=(failcount+1,acc)
            c.execute('''UPDATE PROTEINGBS SET failcount = (?) WHERE acc = (?)''',update_tuple)

        print(f'could not download {acc}')

def getprotein_gbsrs(dbfpath,email,api_key,refresh=False,retry_fails=False):
    Entrez.email=email
    Entrez.api_key=api_key
    conn=sqlite3.connect(dbfpath)
    atexit.register(conn.close)
    conn.row_factory=sqlite3.Row
    c=conn.cursor()
    try:
        c.execute('''SELECT COUNT (*) FROM CAZYDATA''')
        cazydbsize=c.fetchone()[0]
    except:
        return "no db found"
    c.execute('''SELECT acc FROM CAZYDATA''')
    cazyaccs=[x['acc'] for x in c.fetchall()]
    print(f'{len(cazyaccs)} entries in CAZYDATA table')

    c.execute('''CREATE TABLE IF NOT EXISTS PROTEINGBS (acc text, version text, \
                 dldate text, checksum text, pklgbsr glob, failcount int)''')
    c.execute('''SELECT acc FROM PROTEINGBS WHERE failcount=0''')
    cur_pgbaccs=[x['acc'] for x in c.fetchall()]
    print(f'{len(cur_pgbaccs)} existing (successful) entries in CAZYDATA table')

    dlaccs=list(set(cazyaccs).difference(cur_pgbaccs))
    print(f'{len(dlaccs)} remaining CAZY seqs to pull')
    for snum,dlacc in enumerate(dlaccs):
        updatedb(c,dlacc,0)
        if snum>0 and snum%25==0:
            print(f'through {snum} sequences')
            conn.commit()
            time.sleep(5)
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


#    print('Now downloading fasta protein sequences through Biopython-implementation of Entrez eutil API')
#
#    aHT=entrez_requests.pullgb_fromcazyobjs(czes_,email)



#    for subfam in aHT.keys():
#        outfpath=os.path.join(outfolder,f"GH{ghfam}_{subfam}seqs.fasta")
#        SeqIO.write(aHT[subfam],outfpath,"fasta")

if __name__=="__main__":
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('ghfam',help='gh family number (e.g., 5)')
    parser.add_argument('email',help='email address to provide Entrez e-utilities')
    parser.add_argument('--outfolder','-o',help='location of output',default='Seqs_GH<ghfam>')
    parser.add_argument('--force','-f',help='force write of output into existing folder',action='store_true')
    args=parser.parse_args()
    if args.outfolder=='Seqs_GH<ghfam>':
        outfldrpath=f'Seqs_GH{args.ghfam}'
    else:
        outfldrpath=args.outfolder
    if os.path.exists(outfldrpath) and not args.force:
        exit(f'folder {outfldrpath} exists. Use --outfolder=<fpath> to specify new location \
                 or -f option to force overwrite in existing location')
    grab_cazyseqs(args.ghfam,args.email,outfldrpath)


