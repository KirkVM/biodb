import argparse,os# sqlite3,os,sys,argparse,yaml
from Bio import SeqIO
import kazykools# import pullfamily
import entrez_requests

def grab_cazyseqs(ghfam,email,outfolder):
    if not os.path.exists(outfolder):
        #os.mkdir(os.path.join(os.getcwd(),outfolder))
        os.mkdir(outfolder)
    print('starting CAZY scrape')
    czes_=kazykools.scrape_cazyfam(f'GH{ghfam}')
    print('Now downloading fasta protein sequences through Biopython-implementation of Entrez eutil API')
    aHT=entrez_requests.pullgb_fromcazyobjs(czes_,email)
    for subfam in aHT.keys():
        outfpath=os.path.join(outfolder,f"GH{ghfam}_{subfam}seqs.fasta")
        SeqIO.write(aHT[subfam],outfpath,"fasta")

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


