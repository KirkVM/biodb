import requests,re,difflib,bs4
#from lxml import html,etree
from Bio import Entrez,SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import codecs
def scrape_ncbi_genomepage(localhtmlfpath):
    catHT={}
#    page=open(localhtmlfpath,'r')
    page=codecs.open(localhtmlfpath,'r','utf-8')
#    page=page.encode('utf-8')
    soup=bs4.BeautifulSoup(page,'lxml')#,from_encoding="utf-8")#,ignore=True)
    body=str(soup.body)

    trangeRE=re.compile("TemperatureRange:(\w+)(<|\s|,)")
    shapeRE=re.compile("Shape:(\w+)(<|\s|,)")
    motilityRE=re.compile("Motility:(\w+)(<|\s|,)")
    o2RE=re.compile("OxygenReq:(\w+)(<|\s|,)")
    gramRE=re.compile("Gram:(\w+)(<|\s|,)")
    trobj=trangeRE.search(body)
    if trobj:
        catHT['TRange']=trobj.group(1)
    srobj=shapeRE.search(body)
    if srobj:
        catHT['Shape']=srobj.group(1)
    mrobj=motilityRE.search(body)
    if mrobj:
        catHT['Motility']=mrobj.group(1)
    orobj=o2RE.search(body)
    if orobj:
        catHT['Oxygen']=orobj.group(1)
    grobj=gramRE.search(body)
    if grobj:
        catHT['Gram']=grobj.group(1)
    return catHT


def match_to_protein(dnaseq,protseq,accid,matchreq="perfect"):
    dnaseq0=dnaseq
    dstarts_=[0,1,2,3,4,5,6,7,8]
    dstops_=[-6,-5,-4,-3,-2,-1,0]
    relstart=None
    relstop=None
    protstr=str(protseq)
    bestdiff=2
    sloppymatch=False
    longestmatch=len(protstr)-2
    for dstart in dstarts_:
        for dstop in dstops_:
            if dstop==0:
                newseq=dnaseq[dstart:]
            else:
                newseq=dnaseq[dstart:dstop]
            if len(newseq)%3==0:
                translated=newseq.translate()
                tstr=str(translated)
                if tstr[0]=='V' and newseq[0:3]=='GTG': #probably alt. bacterial start codon
                    if protstr[0]=='M':
                        tstr='M'+tstr[1:]
                smatcher=difflib.SequenceMatcher(a=tstr,b=protstr,autojunk=False)
                longestmatch=smatcher.find_longest_match(0,len(tstr),0,len(protstr))
                sizediff=len(protstr)-longestmatch.size
                if sizediff<bestdiff:
                    bestdiff=sizediff
                    relstart=dstart+(longestmatch.a*3)
                    relstop=relstart+longestmatch.size*3
                if matchreq=="sloppy":
                    gdnufs_=difflib.get_close_matches(protstr,[tstr],cutoff=float(len(protstr)-2)/float(len(protstr)))
                    if len(gdnufs_)>0:
                        relstart=dstart
                        relstop=dstart+len(protstr)*3
                        sloppymatch=True
    if bestdiff<2:#relstart is not None and relstop is not None:
        dnaseq=dnaseq[relstart:relstop]
        dnasr=SeqRecord(dnaseq,id=accid)
        if bestdiff>0:
            print(bestdiff)
        return dnasr
    elif sloppymatch:
        dnaseq=dnaseq[relstart:relstop]
        dnasr=SeqRecord(dnaseq,id=accid)
        print('slop but ok')
        return dnasr
    elif matchreq=="perfect": #one more try, allow 1 missing
        dnasr=match_to_protein(dnaseq0,protseq,accid,matchreq="sloppy")
        return dnasr
    else:
        return None



nucaccRE=re.compile('([A-Z]{1,3}_{0,1}\d{5,9}(\.\d)*)') #very basic,just to get rid of 'accession' or other extra txt
nargnRE=re.compile('(complement\()*(join\()*((.+):(.+))+') #greedy!
rgnRE=re.compile('([A-Z\d+\._]+):<*(\d+)\.\.>*(\d+)') #include complement thingy
def getfullseqs(srid):
    Entrez.email="kirk.vandermeulen@wisc.edu"
    #myProteinSequence=Seq(str(sr.seq),IUPAC.protein)

    gbhandle=Entrez.efetch(db="protein",id=srid,rettype="gp",retmode="text")
    ncbigbsr=list(SeqIO.parse(gbhandle,"genbank"))[0]

    gbhandle_forparse=Entrez.efetch(db="protein",id=srid,rettype="gp",retmode="xml")
    gbrec=Entrez.read(gbhandle_forparse)
    try:
        gbrecHT=gbrec[0]
    except: #total failure...
        return None, None, None

    if len(gbrec)!=1:print("{0} values in gbrec for {1}".format(len(gbrec),srid))
    nucacc=None
    if 'GBSeq_source-db' in gbrec[0].keys():
        nucacc_entry=gbrec[0]['GBSeq_source-db']
        naobj=nucaccRE.search(nucacc_entry)
        if naobj:
            nucacc=naobj.group(1)
        else:
            print('couldnt find nuc acc in string {}'.format(nucacc_entry))

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
        codingseq=Seq(dnaseq,IUPAC.unambiguous_dna)

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
        dnasr=match_to_protein(codingseq,ncbigbsr.seq,srid)
        if dnasr is None:
            print('no matchy')
            psr=None
        else:
            psr=SeqRecord(dnasr.seq.translate(),id=dnasr.id)
        return dnasr,psr,ncbigbsr
    else:
        return None,None,ncbigbsr
