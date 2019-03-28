import os,time
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sqlite3,atexit

def _getgbseqrec(gbid):
	handle=Entrez.efetch(db="protein",id=gbid,rettype="fasta",retmode="xml")
	a=Entrez.read(handle)
	sr=SeqRecord(Seq(a[0]['TSeq_sequence'],IUPAC.protein),id=gbid,description=a[0]['TSeq_defline'])
	return sr

def getgbpsr(gbacc):
	handle=Entrez.efetch(db="protein",id=gbacc,rettype="gb",retmode="text")
#    print(handle)
	sr=SeqIO.read(handle,"genbank")
	return sr

def getgbsrs(gbaccs,email,api_key,pause_scheme='default'):
    """public call to ensure playing nice"""
    Entrez.email=email
    if api_key is not None:
        Entrez.api_key=api_key
    srs=[]
    for x,gbacc in enumerate(gbaccs):
        sr=_getgbpsr(gbacc)
        srs.append(sr)
        if pause_scheme=='default':
            if x>0 and x%10==0:
                time.sleep(5)
    return srs