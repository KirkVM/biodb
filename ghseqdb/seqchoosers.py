import re,pickle,sqlite3
from ete3 import Tree,NCBITaxa
import numpy as np
import xarray as xr

from . import seqdbutils
#class SeqFamData:
#    def __init__(self):
#        self.seqxds=None #sequence dataset
#        self.tree=None
#        self.db=None



def get_taxdict(sr,ncbitaxa):
    rankinfodict={}
    for sf in sr.features:
        if sf.type=='source':
            for k in sf.qualifiers:
                if k=='db_xref':
                    dbxrs=sf.qualifiers[k]
                    for dbxr in dbxrs:
                        svals=dbxr.split(':')
                        if svals[0]=='taxon':
                            taxid=int(svals[1])
                            taxlineage=ncbitaxa.get_lineage(taxid)
                            rankdict=ncbitaxa.get_rank(taxlineage)
                            valdict=ncbitaxa.get_taxid_translator(taxlineage)
                            rev_rankdict={rankdict[k]:k for k in rankdict}
                            rev_valdict={valdict[k]:k for k in valdict}
                            rankinfodict={rankdict[k]:[k,valdict[k]] for k in rankdict}
    return rankinfodict

def build_seqxds(vpwxra_pklfpathstr,dbpathstr):
    '''builds an xarray dataset from vscurate xarray, adding '''
    vpwxra=None
    with open(vpwxra_pklfpathstr,'rb') as f:
        vpwxra=pickle.load(f)
    for x in vpwxra.coords['curateseqs']:
        vpwxra.loc[:,x,'normscore']= \
            (vpwxra.loc[:,x,'score']-vpwxra.loc[:,x,'score'].min())/(vpwxra.loc[:,x,'score'].max()-vpwxra.loc[:,x,'score'].min())
    
    taxra=xr.DataArray(  np.full((len(vpwxra.dbseqs),7),np.nan), \
        coords=[vpwxra.dbseqs,['superkingdom','phylum','class','order','family','genus','species']], \
        dims=['dbseqs','ranks'])

    ncbitaxa=NCBITaxa()
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    conn.row_factory=sqlite3.Row
    c=conn.cursor()
    tonamedict={}
    for accentry in vpwxra.dbseqs:
        acc=accentry.data.item(0)
        c.execute('''SELECT * FROM PROTEINGBS WHERE acc=(?)''',(acc,))
        sr=pickle.loads(c.fetchone()['pklgbsr'])
        taxdict=get_taxdict(sr,ncbitaxa)
        for k in taxdict:
            if k in taxra.ranks:
                taxra.loc[acc,k]=taxdict[k][0]
                tonamedict[taxdict[k][1]]=taxdict[k][0]

    conn.close()
    mergeds=xr.Dataset(data_vars={'vpwxra':vpwxra,'taxra':taxra})
    return mergeds

def build_etetree(tfpathstr,ghdb):
    '''builds an xarray dataset from vscurate xarray, adding annotation data'''
    t=Tree(tfpathstr)

    ncbitaxa=NCBITaxa()
