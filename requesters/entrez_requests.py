import os,time,re
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sqlite3,atexit

def _getgbseqrec(gbid):
	handle=Entrez.efetch(db="protein",id=gbid,rettype="fasta",retmode="xml")
	a=Entrez.read(handle)
	sr=SeqRecord(Seq(a[0]['TSeq_sequence'],IUPAC.protein),id=gbid,description=a[0]['TSeq_defline'])
	return sr

def getgbpsr(gbacc):
    try:
        handle=Entrez.efetch(db="protein",id=gbacc,rettype="gb",retmode="text")
        sr=SeqIO.read(handle,"genbank")
    except:
        sr=None
    return sr


def getgbsrs(gbaccs,email,api_key,pause_scheme='default'):
    """public call to ensure playing nice"""
    Entrez.email=email
    Entrez.api_key=api_key
    srs=[]
    for x,gbacc in enumerate(gbaccs):
        sr=_getgbpsr(gbacc)
        srs.append(sr)
        if pause_scheme=='default':
            if x>0 and x%10==0:
                time.sleep(5)
    return srs


def getdnasrs(gbaccs,email,api_key,pause_scheme='default'):
    Entrez.email=email
    Entrez.api_key=api_key
    srs=[]
    for x,gbacc in enumerate(gbaccs):
        sr=getdnasr(gbacc)
        if sr is not None:
            srs.append(sr)
        if pause_scheme=='default':
            if x>0 and x%10==0:
                time.sleep(2)
    return srs

def getdnasr(srid):
    nucaccRE=re.compile('([A-Z]{1,3}_{0,1}\d{5,9}(\.\d)*)') #very basic,just to get rid of 'accession' or other extra txt
    nargnRE=re.compile('(complement\()*(join\()*((.+):(.+))+') #greedy!
    rgnRE=re.compile('([A-Z\d+\._]+):<*(\d+)\.\.>*(\d+)') #include complement thingy

    gbhandle=Entrez.efetch(db="protein",id=srid,rettype="gp",retmode="xml")
    gbrec=Entrez.read(gbhandle)
    try:
        gbrecHT=gbrec[0]
    except: #total failure...
        print('couldnt get gbrecHT')
        return None

    assert (len(gbrec)==1), print(f'{len(gbrec)} values for {srid}')
    nucacc=None
    if 'GBSeq_source-db' in gbrec[0].keys():
        nucacc_entry=gbrec[0]['GBSeq_source-db']
        naobj=nucaccRE.search(nucacc_entry)
        if naobj:
            nucacc=naobj.group(1)
        else:
            print('couldnt find nuc acc in string {}'.format(nucacc_entry))
            return None

    cdsHT=None
    if 'GBSeq_feature-table' in gbrec[0].keys():
        for ftrindex,ftrel in enumerate(gbrec[0]['GBSeq_feature-table']):
            if 'GBFeature_key' in ftrel.keys():
                if ftrel['GBFeature_key']=='CDS':
                    if cdsHT is not None:
                        print('2 cds elements? for {}'.format(srid))
                    cdsHT=ftrel

    #now find the coding region:
    coding_data=None
    transl_data=None
    explicit_codonstart=None
    if cdsHT is not None:
        #print(cdsHT)
        if 'GBFeature_quals' in cdsHT.keys():
            for kvHT in cdsHT['GBFeature_quals']:
                if kvHT['GBQualifier_name']=='coded_by':
                    coding_data=kvHT['GBQualifier_value']
                if kvHT['GBQualifier_name']=='transl_table':
                    transl_data=kvHT['GBQualifier_value']
                if kvHT['GBQualifier_name']=='codon_start':
                    explicit_codonstart=int(kvHT['GBQualifier_value'])-1

    complement=False
    join=False
    cdstarts_=[]
    cdstops_=[]
    if coding_data is not None:
        accrgnobj=nargnRE.match(coding_data)
        if accrgnobj:
            if accrgnobj.group(1) is not None:
                complement=True
            if accrgnobj.group(2) is not None:
                join=True
            rgnstr=accrgnobj.group(3)
            rgnobjs_=rgnRE.findall(rgnstr)
            for rgnobj in rgnobjs_:
                new_nucacc=rgnobj[0]
                cdstarts_.append(int(rgnobj[1])-1)
                cdstops_.append(int(rgnobj[2]))
                if nucacc is not None and nucacc!=new_nucacc:
                    nucacc_trimmed=nucacc.split('.')[0]
                    newacc_trimmed=new_nucacc.split('.')[0]
                    if nucacc_trimmed!=newacc_trimmed:
                        print('which na acc: {0} or {1}'.format(nucacc,new_nucacc))
                    print('using {0}'.format(new_nucacc))
                    nucacc=new_nucacc

    #if explicit_codonstart is not None and (complement or join): #figuring out how to handle this...
    #if explicit_codonstart and join and complement: #figuring out how to handle this...
    #    print('what do we do here? quitting for now....<<<-----...')
    #    break

    fulldnastr=None
    if nucacc is not None and len(cdstarts_)>0:
        dnahandle=Entrez.efetch(db="nucleotide",id=nucacc,rettype="gb",retmode="xml")
        dnarec=Entrez.read(dnahandle)
        if len(dnarec)!=1:
            print("{0} DNA values in dnarec for {1}({2})".format(len(dnarec),srid,nucacc))
        if len(dnarec)>0:
            if 'GBSeq_sequence' in dnarec[0].keys():
                fulldnastr=dnarec[0]['GBSeq_sequence']
        else:
            print('...cant find dnaseq.....skipping........')
            return None


    codingseq=None
    if fulldnastr:
        merge_region=0
        dnaseq=''
        for start,stop in zip(cdstarts_,cdstops_):
            if complement: #treat differently in case of explicit_codonstart
                if explicit_codonstart is not None and merge_region==0: #override!
                    start=int(start)
                    stop=int(stop)-(explicit_codonstart) #b/c apparently codon_start refers to AFTER revcomp()
                if merge_region==0:
                    start=max(0,int(start)-3)
                if merge_region==len(cdstarts_)-1:
                    stop=min(int(stop)+3,len(fulldnastr))
            else:
                if explicit_codonstart is not None: #override!
                    if merge_region==0:
                        start=start+int(explicit_codonstart)
                    stop=int(stop)
                if merge_region==0:
                    start=max(0,int(start)-3)
                if merge_region==len(cdstarts_)-1:
                    stop=min(int(stop)+3,len(fulldnastr))
            dnaseq+=fulldnastr[start:stop].upper()
            #print(len(dnaseq))
            merge_region+=1
        codingseq=Seq(dnaseq,alphabet=generic_dna)

    good_dna=False
    if codingseq:
        dnastr=str(codingseq)
        gcount=dnastr.count('G')
        ccount=dnastr.count('C')
        tcount=dnastr.count('T')
        acount=dnastr.count('A')
        if gcount+ccount+tcount+acount!=len(dnaseq):
            print('some weirdos in there. skipping.........')
            dnaseq=dnaseq.replace('Y','T') #t or c
            dnaseq=dnaseq.replace('S','G') #g or c
            dnaseq=dnaseq.replace('M','A') #a or c
            dnaseq=dnaseq.replace('W','G') #g or c
        else:
            good_dna=True

    if good_dna:
        if complement:
            codingseq=codingseq.reverse_complement()
        dnasr=SeqRecord(codingseq,id=srid)
        return dnasr
    else:
        return None