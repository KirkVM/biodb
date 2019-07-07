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