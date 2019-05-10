import re,pickle,sqlite3,math
from ete3 import Tree,NCBITaxa
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import t as student

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

def plot_ebounds(mds,dbpathstr,dbseqidx,scorefield='normscore',start_stops='avg',ax=None):
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
    conn=seqdbutils.gracefuldbopen(dbpathstr)
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

def ebound_xtaleval(mds,dbpathstr,scorefield='normscore',start_stops='avg'):
    mdsxtal=mds.drop([x for x in mds.dbseq.values if x not in mds.curateseq.values],dim='dbseq')
    numseqs=len(mdsxtal.dbseq)
    numrows=math.ceil(numseqs/3)
    fig,axs=plt.subplots(numrows,3,figsize=(30,7*numrows))
    for seqnum,ax in zip(range(numseqs),np.ravel(axs)):
        plot_ebounds(mds,dbpathstr,seqnum,ax=ax)

def get_maxesize(dbpathstr):
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    c=conn.cursor()
    c.execute('''SELECT * FROM PROTEINGBS''')
    maxlength=0
    rows=c.fetchall()
    for row in rows:
        gbsr=pickle.loads(row['pklgbsr'])
        maxlength=max(maxlength,len(gbsr))
    conn.close()
    return maxlength

def get_esizeprior(dbpathstr,max_esize=5000):
    if max_esize is None: #find it ourselves
        maxlength=get_maxesize(dbpathstr)        
    else:
        maxlength=max_esize
    
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    c=conn.cursor()
    c.execute('''SELECT * FROM CURATEDXTALS''')
    rows=c.fetchall()
    cc_sizes=[]
    for row in rows:
        start=row['pgbsr_ccstart']
        stop=row['pgbsr_ccstop']
        cc_sizes.append(stop-start+1)
    ccmed=int(round(np.median(cc_sizes),0))
    conn.close()

#    maxset=student(v1,60,loc=0,scale=40)
    lseg=np.array([1,ccmed+1])-np.array([0,ccmed])
    probs=[]
    coolarr=np.zeros((max_esize,max_esize))
    probs=[]
    for yoffset in range(121):
        v1=np.cross(lseg,np.array([0,ccmed+yoffset])-np.array([0,ccmed]))/np.linalg.norm(lseg)
        v2=student.pdf(v1,60,loc=0,scale=40)
        probs.append(v2)
        for xpos in range(max_esize):
            if (ccmed-yoffset+xpos >=0) and (ccmed-yoffset+xpos<max_esize):
                coolarr[ccmed-yoffset+xpos,xpos]=v2
            if ccmed+yoffset+xpos<max_esize:
                coolarr[ccmed+yoffset+xpos,xpos]=v2
    minset=np.min(probs)
    unders=coolarr<minset
    coolarr[unders]=minset
    return coolarr

def calc_prob_grid(acc,esize_prior,dbpathstr,mds,motif_dict={'PF00150':[16,32]}):
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    c=conn.cursor()
    c.execute('''SELECT * FROM PROTEINGBS WHERE acc=(?)''',(acc,))
    gbsr=pickle.loads(c.fetchone()['pklgbsr'])
    genelength=len(gbsr.seq)
    p_array=np.zeros((genelength,genelength,8))
    xyidx=[str(x) for x in range(genelength)]
    probcoords=['x','y','esize_prior','pfam_prior','sum_prior','prob_data','data_likelihood','posterior']
    probdr=xr.DataArray(p_array,coords=[xyidx,xyidx,probcoords],dims=['start','stop','prob'])
    spgrid=np.meshgrid(np.arange(0,genelength,1),np.arange(0,genelength,1))
    probdr.loc[:,:,'x']=spgrid[0]
    probdr.loc[:,:,'y']=spgrid[1]
    probdr.loc[:,:,'esize_prior']=esize_prior[:genelength,:genelength]
    
    maxset=student.pdf(0,20,loc=0,scale=30)
    c.execute('''SELECT * FROM HMMERSEQDATA WHERE acc=(?)''',(acc,))
    hrows=c.fetchall()
    for hrow in hrows:
        if hrow['motifacc'] in motif_dict.keys():
            pfam_begin=hrow['align_start']-hrow['hmm_start']
            pfam_end=hrow['align_stop']+(hrow['motifsize']-hrow['hmm_stop'])-1
            exp_next=motif_dict[hrow['motifacc']][0]
            exp_cext=motif_dict[hrow['motifacc']][1]
            exp_estart=pfam_begin-exp_next
            exp_estop=pfam_end+exp_cext
            print(pfam_begin,exp_estart)
            print(pfam_end,exp_estop)
            student_dict={}
            temparr=np.zeros((spgrid[0].shape))
            for xoffset in range(60):#xg in np.arange(exp_estart-50,exp_estart+50,1):
                for yoffset in range(60):
                    dist2peak=int(round(np.linalg.norm(np.array([exp_estart-xoffset,exp_estop-yoffset])-np.array([exp_estart,exp_estop])),0))
#                    v2=student.pdf(dist2peak,10,loc=0,scale=20)
                    if dist2peak in student_dict.keys():
                        v2=student_dict[dist2peak]
                    else:
                        v2=student.pdf(dist2peak,20,loc=0,scale=30)
                        student_dict[dist2peak]=v2
                    exp_estarts=[exp_estart-xoffset,exp_estart+xoffset]
                    exp_estops=[exp_estop-yoffset,exp_estop+yoffset]
                    for x0 in exp_estarts:
                        for y0 in exp_estops:
                            if x0>=0 and x0<genelength and y0>=0 and y0<genelength:
                                temparr[y0,x0]=v2

            minset=min(1e-3*maxset,np.min([student_dict[x] for x in student_dict.keys()]))
            junk=temparr<minset
            temparr[junk]=minset
            probdr.loc[:,:,'pfam_prior']+=temparr

    conn.close()
#    probdr.loc[:,:,'data_likelihood']+=1e-3
    student_dict={}
    for cs in mds.curateseq:
        temparr=np.zeros((spgrid[0].shape))
        x=int(round(mds.vpwxra_cc.loc[acc,cs,'avg_start'].values.item(0),0))
        y=int(round(mds.vpwxra_cc.loc[acc,cs,'avg_stop'].values.item(0),0))
        ns=mds.vpwxra_cc.loc[acc,cs,'normscore']
        for xg in np.arange(x-10,x+10,1):
            for yg in np.arange(y-10,y+10,1):#,seqgrid[1].ravel()):
                if xg<0 or yg<0 or xg>=spgrid[0].shape[0] or yg>=spgrid[0].shape[0]:
                    continue
                dist2peak=np.linalg.norm(np.array([xg,yg])-np.array([x,y]))
                if dist2peak in student_dict.keys():
                    temparr[yg,xg]=1e3*student_dict[dist2peak]*1e3**(ns-0.1)
                else:
                    studval=student.pdf(dist2peak,5,loc=0,scale=5)
                    student_dict[dist2peak]=studval
                    temparr[yg,xg]=1e3*student_dict[dist2peak]*1e3**(ns-0.1)
                    #temparr[yg,xg]=studval
#        temparr=1e3*studval*1e3**(ns-0.1)
#                studval=student.pdf(np.linalg.norm(np.array([xg,yg])-np.array([x,y])),5,loc=0,scale=5)
#                probdr.loc[:,:,'data_likelihood'][yg,xg]+=1e3*studval*1e3**(ns-0.1)
        probdr.loc[:,:,'prob_data']+=temparr#1e3*temparr*1e3**(ns-0.1)
    probdr.loc[:,:,'sum_prior']=probdr.loc[:,:,'esize_prior']+probdr.loc[:,:,'pfam_prior']
    
    probdr.loc[:,:,'data_likelihood']=probdr.loc[:,:,'prob_data']*probdr.loc[:,:,'sum_prior']#1e3*temparr*1e3**(ns-0.1)
    probdr.loc[:,:,'posterior']=probdr.loc[:,:,'data_likelihood']/np.sum(probdr.loc[:,:,'data_likelihood'])
    return probdr
#        hs.append(v2)   