import re,pickle,sqlite3
from ete3 import Tree,NCBITaxa
import numpy as np
import xarray as xr
import pandas as pd

from . import seqdbutils
#class SeqFamData:
#    def __init__(self):
#        self.seqxds=None #sequence dataset
#        self.tree=None
#        self.db=None

def get_intpwid_df(pwfpathstr,format='clustal'):
    pwdf=pd.read_csv(pwfpathstr,sep='\s+',index_col=0,skiprows=1,header=None)
    
    accvrsnRE=re.compile('(\w+)(\..+)*')
    idxrename={}
    for idxname in pwdf.index:
        reobj=accvrsnRE.match(idxname)
        if reobj:
            idxrename[idxname]=reobj.group(1)
        else:
            idxrename[idxname]=idxname#reobj.group(1)
    pwdf.rename(idxrename,axis='index',inplace=True)
    pwdf.columns=list(pwdf.index)
    pwdf.astype(int)
    return pwdf



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

def build_seqxds(vpwxracc_fpathstr,dbpathstr,vpwxrafull_fpathstr=None):
    '''builds an xarray dataset from vscurate xarray, adding '''
    vpwxra_cc=xr.open_dataarray(vpwxracc_fpathstr)
    for x in vpwxra_cc.curateseq.values:#().coords[:,'curateseqs'].data:
        vpwxra_cc.loc[:,x,'normscore'] =  (vpwxra_cc.loc[:,x,'score']- vpwxra_cc.loc[:,x,'score'].min()) /    \
             (vpwxra_cc.loc[:,x,'score'].max() - vpwxra_cc.loc[:,x,'score'].min())
    mergeds=xr.Dataset(data_vars={'vpwxra_cc':vpwxra_cc})
    if vpwxrafull_fpathstr is not None:
        vpwxra_full=xr.open_dataarray(vpwxrafull_fpathstr)
        for x in vpwxra_full.curateseq.values:#().coords[:,'curateseqs'].data:
            vpwxra_full.loc[:,x,'normscore'] =  (vpwxra_full.loc[:,x,'score']- vpwxra_full.loc[:,x,'score'].min()) /    \
                (vpwxra_full.loc[:,x,'score'].max() - vpwxra_full.loc[:,x,'score'].min())
        mergeds['vpwxra_full']=vpwxra_full#xr.Dataset(data_vars={'vpwxra_cc':vpwxracc})

    taxra=xr.DataArray(  np.full((len(vpwxra_cc.dbseq),7),np.nan), \
        coords=[vpwxra_cc.dbseq,['superkingdom','phylum','class','order','family','genus','species']], \
        dims=['dbseq','ranks'])

    ncbitaxa=NCBITaxa()
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    conn.row_factory=sqlite3.Row
    c=conn.cursor()
    tonamedict={}
    for accentry in vpwxra_cc.dbseq:
        acc=accentry.data.item(0)
        c.execute('''SELECT * FROM PROTEINGBS WHERE acc=(?)''',(acc,))
        sr=pickle.loads(c.fetchone()['pklgbsr'])
        taxdict=get_taxdict(sr,ncbitaxa)
        for k in taxdict:
            if k in taxra.ranks:
                taxra.loc[acc,k]=taxdict[k][0]
                tonamedict[taxdict[k][1]]=taxdict[k][0]

    conn.close()
#    mergeds=xr.Dataset(data_vars={'vpwxra':vpwxra,'taxra':taxra})
    #mergeds=xr.Dataset(data_vars={'vpwxra':vpwxracc})
    mergeds['taxra']=taxra
    return mergeds
