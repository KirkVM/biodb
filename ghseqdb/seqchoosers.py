import re,pickle,sqlite3
from ete3 import Tree,NCBITaxa
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

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
#        sr=pickle.loads(c.fetchone()['pklgbsr'])
        #taxdict=get_taxdict(sr,ncbitaxa)
        row=c.fetchone()
        taxid=row['taxid']
        if taxid is not None:
            taxlineage=ncbitaxa.get_lineage(taxid)
            rankdict=ncbitaxa.get_rank(taxlineage)
            valdict=ncbitaxa.get_taxid_translator(taxlineage)
            rankinfodict={rankdict[k]:[k,valdict[k]] for k in rankdict}
            for k in rankinfodict:
                if k in taxra.ranks:
                    taxra.loc[acc,k]=rankinfodict[k][0]
                    tonamedict[rankinfodict[k][1]]=rankinfodict[k][0]
    conn.close()
#    mergeds=xr.Dataset(data_vars={'vpwxra':vpwxra,'taxra':taxra})
    #mergeds=xr.Dataset(data_vars={'vpwxra':vpwxracc})
    mergeds['taxra']=taxra
    return mergeds

import matplotlib.pyplot as plt
def plot_ebounds(mds,dbfpath,dbseqidx,scorefield='normscore',start_stops='avg',ax=None):
    assert (scorefield in ['normscore','score'])
    assert (start_stops in ['avg','mdn'])
    start_yid=f'{start_stops}_start'
    stop_yid=f'{start_stops}_stop'
    if ax is None:
        _,ax=plt.subplots(1,1)
    seqid=mds.dbseq.values[dbseqidx]
    ax.plot(mds.vpwxra_cc.loc[seqid,:,scorefield],mds.vpwxra_cc.loc[seqid,:,start_yid],'g.',marker='o',ms=10)
    ax.plot(mds.vpwxra_cc.loc[seqid,:,scorefield],mds.vpwxra_cc.loc[seqid,:,stop_yid],'g.',marker='s',ms=10)
    ax.plot(mds.vpwxra_full.loc[seqid,:,scorefield],mds.vpwxra_full.loc[seqid,:,start_yid],'b.',marker='o',ms=7)
    ax.plot(mds.vpwxra_full.loc[seqid,:,scorefield],mds.vpwxra_full.loc[seqid,:,stop_yid],'b.',marker='s',ms=7)

    fycoord=.55
    ax.text(0.3,fycoord,mds.dbseq.values[dbseqidx],size=20,transform=ax.transAxes)
    fycoord-=.06
    ax.grid(True,axis='y')
    conn=seqdbutils.gracefuldbopen(dbfpath)
    c=conn.cursor()

    c.execute('''SELECT * FROM PROTEINGBS WHERE acc=(?)''',(seqid,))
    pgbrow=c.fetchone()
    pgb=pickle.loads(pgbrow['pklgbsr'])
    ax.text(0.3,fycoord,f'seqlen={len(pgb.seq)}',size=14,transform=ax.transAxes)
    c.execute('''SELECT * FROM HMMERSEQDATA WHERE acc=(?)''',(seqid,))
    hrows=c.fetchall()
    for hrow in hrows:
        fycoord-=.06
        alstartenv=hrow['align_start']-hrow['hmm_start'] #switches to zero-based
        alstopenv=hrow['align_stop']-(hrow['motifsize']-hrow['hmm_stop']+1) #switches to zero-based
        #htext=f"{hrow['motifacc']}: {hrow['align_start']}-{hrow['align_stop']} ({hrow['env_start']}-{hrow['env_stop']} )"
        htext=f"{hrow['motifacc']}: {hrow['align_start']}-{hrow['align_stop']} ({alstartenv}-{alstopenv} )"
        ax.text(0.3,fycoord,htext,size=14,transform=ax.transAxes)
    conn.close()
    
