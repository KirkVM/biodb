import requests,os,re,codecs
#from lxml import html,etree
import bs4
from io import BytesIO,StringIO
import pdb
#if straight from web:############
#pagelink="https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=2513237176"
#page=requests.get(pagelink)
#tree=html.fromstring(page.content)
#######################################


metadatafields_=['Temperature Range','Oxygen Requirement','Isolation','Isolation Country']
metadatafields_.extend(['Ecosystem','Ecosystem Category','Ecosystem Type','Habitat','Host Name','Biotic Relationships'])
metadatafields_.extend(['Specific Ecosystem','Sporulation','Sample Body Site','Sample Body Subsite'])

def scrape_jgi_imgpage(localhtmlfpath):
#    infoRE=re.compile("NCBI Taxon ID")
#    myparser=etree.HTMLParser()#encoding="UTF-8")
#    tree=html.fromstring(open(localhtmlfpath,'r').read(),parser=myparser)
    catHT={}
    overlaps_=[]
#    page=open(localhtmlfpath,'r')
    page=codecs.open(localhtmlfpath,'r','utf-8')

    soup=bs4.BeautifulSoup(page,'lxml')#'html.parser'
    #body=soup.body
    tables_=soup.find_all('table')
    for table in tables_:
        trs_=table.find_all('tr')
        for tr in trs_:
            th=tr.find('th')
            if tr.find('th') and tr.find('td'):
                key=tr.find('th').string.strip()
                ditem=tr.find('td')
                if ditem.string is not None:
                    val=ditem.string.strip()
                else:
                    linkout=ditem.find('a')
                    if linkout is not None:
                        if linkout.string is not None:
                            val=linkout.string.strip()
                        else:
                            val=None
                    else:
                        if key!='Project Geographical Map':
                            print('cant find info for {0} in {1}'
                            .format(key,localhtmlfpath))
                        val=None
                if val is not None and len(val)>0:
                    catHT[key]=val
    return catHT
