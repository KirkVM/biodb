import requests
from lxml import html,etree
import re,time
#from Bio import Entrez
#from Bio.SeqRecord import SeqRecord
#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
#from Bio import SeqIO
from io import open
import bs4

class CAZYEntry:
    def __init__(self):
        self.gbids_=[]
        self.gblinks_=[]
        self.name=None
        self.name_raw=None
        self.pdbids_=[]
        self.pdblinks_=[]
        self.uniprotids_=[]
        self.uniprotlinks_=[]
        self.organism=None
        self.orglink=None
        self.family=None
        self.kingdom=None
        self.ecs_=[]

def parse_fam(data):
    famnum=None
    if data.string is not None:
        famnum=data.string
    return famnum

def parse_ec(data):
    ecs_=[]
    for el in data.find_all('a'):
        ecs_.append(el.string)
    return ecs_

def parse_name(data):
    names_=list(data.strings)
    fullname=names_[0]
    if len(names_)>1:
        nickname=names_[1].strip()[1:-1] #pulls parentheses
    else:
        nickname=None
    return fullname,nickname

def parse_org(data):
    hlink=None
    textval=None
    hlink=data.a.get('href')
    if data.b is not None:
        textval=data.b.string
    elif data.a.string is not None:
        textval=data.a.string
    return hlink,textval

def parse_gb(data):
    hlink=None
    gbtextvals_=None
    hlink=data.a.get('href')
    gbtextvals_=list(data.strings)
    return hlink,gbtextvals_
	
def parse_uniprot(data):
    hlink=None
    textval=None
    if data.a is not None:
        hlink=data.a.get('href')
        textval=data.a.string
    return hlink,textval

def parse_pdb(data):
    hlinks_=[]
    textvals_=[]
    for pdbel in data.find_all('a'):
        hlinks_.append(pdbel.get('href'))
        textvals_.append(pdbel.string)
    return hlinks_,textvals_

def parse_row(row,kingdom=None):
    cze=CAZYEntry()	
    cze.kingdom=kingdom
    tds_=row.find_all('td')
    name,nickname=parse_name(tds_[0])
    ecs_=parse_ec(tds_[1])
    orglink,orgtext=parse_org(tds_[2])
    gblink,gbtexts_=parse_gb(tds_[3])
    uplink,uptext=parse_uniprot(tds_[4])
    pdblinks_,pdbtexts_=parse_pdb(tds_[5])
    famnum=parse_fam(tds_[6])

    cze.name=name
    cze.organism=orgtext
    cze.orglink=orglink
    cze.ecs_=ecs_
    cze.gbids_=gbtexts_
    cze.gblinks_.append(gblink)
    cze.uniprotids_.append(uptext)
    cze.uniprotlinks_.append(uplink)
    cze.pdbids_=pdbtexts_
    cze.pdblinks_=pdblinks_
    cze.family=famnum
    return cze



def find_pagins(rowel,addlpagesRE):
    pagins_=[]
    for linkel in rowel.find_all('a'):
        linkstr=linkel.get('href')
        if linkstr is not None:
            pobj=addlpagesRE.search(linkstr)
            if pobj:
                pagins_.append(pobj.group(1))
    return pagins_


def parsepage(plinkstr,kingdom=None):
    czes_=[]
    accRE=re.compile("protein\&.*val=(.+)")
    page=requests.get(plinkstr)
    soup=bs4.BeautifulSoup(page.content,"lxml")
    
    ###for testing local page#####
    #page=open('CAZy - GH5.html','r')
    #soup=bs4.BeautifulSoup(page,"lxml")
    ########################################
    kings_=['Archaea','Bacteria','Eukaryota','Viruses','unclassified']
    body=soup.body
    for row in body.find_all('tr'):#.children:
#        if True not in [
        tds_=row.find_all('td')
        #lets make sure this is a row with an accession code (do special stuff here?...)
        if True not in [accRE.search(str(td)) is not None for td in tds_]:
            #lets see if the row reflects new kingdom:
            #print(row.prettify())
            s_=row.find_all('span')
            if len(s_)>0:
                if s_[0].string in kings_:
                    kingdom=s_[0].string
            continue
        czes_.append(parse_row(row,kingdom=kingdom))
        #print(row.prettify())
    return czes_


def scrape_cazyfam(famname):
    accRE=re.compile("protein\&.*val=(.+)")
    famselect=famname+"_all"
    addlpagesRE=re.compile("{0}\.html\?debut_PRINC=(\d+)\#pagination_PRINC".format(famselect))
    mainpagelink="http://www.cazy.org/{0}.html".format(famselect)
    page=requests.get(mainpagelink)
    soup=bs4.BeautifulSoup(page.content,"lxml")
    ###for testing local page#####
    #page=open('CAZy - GH5.html','r')
    #soup=bs4.BeautifulSoup(page,"lxml")
    ################################
    
    body=soup.body   
    pagins_=[]
    
    for row in body.find_all('tr'):#.children:
        tds_=row.find_all('td')
#        #lets make sure this is a row with an accession code (do special stuff here?...)
        if True not in [accRE.search(str(td)) is not None for td in tds_]:
            #promising...:
            pagins_.extend(find_pagins(row,addlpagesRE))
    #this relies on last page number being visible and steps being 1000 per page
    if len(pagins_)>0:
        pageset=set(pagins_)
        orderedps_=sorted([int(p) for p in pageset])
        if len(orderedps_)>1:
#        if ps2add>1:
            ps2add=int((orderedps_[-1]-orderedps_[-2])/1000)
            newps_=[]
            for p in range(1,ps2add):
                p2add=p*1000+orderedps_[-2]
                newps_.append(p2add)
            orderedps_.extend(newps_)
        pagins_=sorted(orderedps_)
    #read-in first page
    czes_=[]
    czes_.extend(parsepage(mainpagelink))
    #print(czes_[0].gbids_)
    #read-in addl pages
    for p in pagins_:
        mplink="http://www.cazy.org/{0}.html?debut_PRINC={1}#pagination_PRINC".format(famselect,p)
        czes_.extend(parsepage(mplink))
        print("read in {0} ids".format(len(czes_)))
    czes_=list(set(czes_))
    print("after duplicate removal: read in {0} ids".format(len(czes_)))
    return czes_


if __name__=="__main__":
#    parsepage("http://www.cazy.org/GH5_all.html")
    scrape_cazyfam('GH5')
