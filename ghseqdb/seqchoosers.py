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
        alstopenv=hrow['align_stop']+(hrow['motifsize']-hrow['hmm_stop']+1) #switches to zero-based
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
        plot_ebounds(mdsxtal,dbpathstr,seqnum,ax=ax)

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

def build_prior_array(median,shape,scale,size):
    coolarr=np.zeros((size,size))
    probs=[]
    lseg=np.array([1,median+1])-np.array([0,median])
    for yoffset in range(2*shape):
        v1=np.cross(lseg,np.array([0,median+yoffset])-np.array([0,median]))/np.linalg.norm(lseg)
        v2=student.pdf(v1,shape,loc=0,scale=scale)
        probs.append(v2)
        for xpos in range(size):
            if (median-yoffset+xpos >=0) and (median-yoffset+xpos<size):
                coolarr[median-yoffset+xpos,xpos]=v2
            if median+yoffset+xpos<size:
                coolarr[median+yoffset+xpos,xpos]=v2
    minset=np.min(probs)
    unders=coolarr<minset
    coolarr[unders]=minset
    return coolarr
 

def get_esizeprior(dbpathstr,max_esize=5000,eform='cc'):
    assert (eform in ['cc','full','both']),'eform argument invalid ("cc","full", or "both")'
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    c=conn.cursor()
    c.execute('''SELECT * FROM CURATEDXTALS''')
    rows=c.fetchall()
    ccsizes=[]
    fullsizes=[]
    for row in rows:
        ccstart=row['pgbsr_ccstart']
        ccstop=row['pgbsr_ccstop']
        ccsizes.append(ccstop-ccstart+1)
        if eform in ['full','both']:
#            try:
            fullstart=row['pgbsr_fullstart']
            fullstop=row['pgbsr_fullstop']
            if abs(fullstart-ccstart)>5 or abs(fullstop-ccstop)>5:
                fullsizes.append(fullstop-fullstart+1)
#            except:
#                breakpoint()
    ccmedian=int(round(np.median(ccsizes),0))
    cc_shape=int((5/6)*max(0.15*max(ccsizes),max(ccsizes)-ccmedian,ccmedian-min(ccsizes)))
    cc_scale=int(2*cc_shape/3)
    if eform in ['full','both']:
        fullmedian=int(round(np.median(fullsizes),0))
        full_shape=int((5/6)*max(0.15*max(fullsizes),max(fullsizes)-fullmedian,fullmedian-min(fullsizes)))
        full_scale=int(2*full_shape/3)
    conn.close()

    cc_arr=build_prior_array(ccmedian,cc_shape,cc_scale,max_esize)
    if eform in ['full','both']:
        full_arr=build_prior_array(fullmedian,full_shape,full_scale,max_esize)
    if eform=='cc':
        return cc_arr
    elif eform=='full':
        return full_arr
    return cc_arr+full_arr

def build_pfam_expmatrix(parray,expstart,expstop):
    student_dict={}
    maxset=student.pdf(0,20,loc=0,scale=30) #what is this????
    for xoffset in range(60):#xg in np.arange(exp_estart-50,exp_estart+50,1):
        for yoffset in range(60):
            dist2peak=int(round(np.linalg.norm(np.array([expstart-xoffset,expstop-yoffset])-np.array([expstart,expstop])),0))
#           v2=student.pdf(dist2peak,10,loc=0,scale=20)
            if dist2peak in student_dict.keys():
                v2=student_dict[dist2peak]
            else:
                v2=student.pdf(dist2peak,20,loc=0,scale=30)
                student_dict[dist2peak]=v2
            exp_estarts=[expstart-xoffset,expstart+xoffset]
            exp_estops=[expstop-yoffset,expstop+yoffset]
            for x0 in exp_estarts:
                for y0 in exp_estops:
                    if x0>=0 and x0<parray.shape[0] and y0>=0 and y0<parray.shape[0]:
                        parray[y0,x0]=v2
    minset=min(1e-3*maxset,np.min([student_dict[x] for x in student_dict.keys()]))
    junk=parray<minset
    parray[junk]=minset
    return parray

addl_motifs={}
core_motifs={'GH5':['PF00150',[16,32]],'GH43':['PF04616',[2,13]],'GH8':['PF01270',[11,8]]} #Glyco_hydro_5,Glyco_hydro_43,Glyco_hydro_8
core_motifs.update({'GH9':['PF00759',[3,9]],'GH48':['PF02011',[24,32]],'AA10':['PF03067',[0,3]]})#Glyco_hydro_9,Glyco_hydro_48
core_motifs.update({'AA9':['PF03443',[0,10]],'GH10':['PF00331',[3,2]],'GH26':['PF02156',[5,17]]})#Glyco_hydro_61,Glyco_hydro_10,Glyco_hydro_26
core_motifs.update({'GH30':['PF02055',[-17,67]],'GH6':['PF01341',[5,5]]})#Glyco_hydro_30,Glyco_hydro_6
addl_motifs.update({'GH5':{'PF18448':[5,5] }}) #CBM_X2, anything for GH5_5s?
addl_motifs.update({'GH43':{'PF16369':[5,5],'PF17851':[5,5],'PF07081':[5,5],'PF14200':[5,5],'PF05270':[5,5],'PF00754':[5,5]} }) 
#GH43_C,GH43_C2,DUF (always seen with 17851),RicinB_lectin_2,AbfB,F5_F8_type_C
addl_motifs.update({'GH8':{}})
addl_motifs.update({'GH9':{'PF02927':[5,5],'PF00942':[5,5]}})#CelD_N,CBM_3
addl_motifs.update({'GH48':{}})
addl_motifs.update({'AA10':{}})
addl_motifs.update({'AA9':{}})
addl_motifs.update({'GH10':{}})
addl_motifs.update({'GH26':{}})
addl_motifs.update({'GH30':{'PF17189':[5,5]}}) #Glyco_hydro_30C
addl_motifs.update({'GH6':{}})
def calc_prob_grid(acc,ghfam,esize_prior,dbpathstr,mds,eform='cc',exclude_self=False):
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
    
    c.execute('''SELECT * FROM HMMERSEQDATA WHERE acc=(?)''',(acc,))
    hrows=c.fetchall()
    conn.close()
    #first iterate through and find cc data...
    exp_core_begin=None
    exp_core_end=None
    mm_cnt=0
    for hrow in hrows:
        motifacc=hrow['motifacc']
        if motifacc == core_motifs[ghfam][0]:
            mm_cnt+=1
            pfam_core_begin=hrow['align_start']-hrow['hmm_start']
            pfam_core_end=hrow['align_stop']+(hrow['motifsize']-hrow['hmm_stop'])-1
            exp_next=core_motifs[ghfam][1][0] #nt-extension
            exp_cext=core_motifs[ghfam][1][1] #ct-extension
            exp_core_begin=pfam_core_begin-exp_next
            exp_core_end=pfam_core_end+exp_cext
            temparr=build_pfam_expmatrix(np.zeros((spgrid[0].shape)),exp_core_begin,exp_core_end)
            probdr.loc[:,:,'pfam_prior']+=temparr
    if mm_cnt>1:
        print('*possible double domain construct')
    addl_motif_accs=list(addl_motifs[ghfam].keys())
    if exp_core_begin is not None and exp_core_end is not None and eform=='full':
        for hrow in hrows:
            motifacc=hrow['motifacc']
            if motifacc in addl_motif_accs:
#                print(motifacc,hrow['align_start'],hrow['align_stop'],hrow['hmm_stop'])
                print(f'pfam aux motif {motifacc}')
                pfam_motif_begin=hrow['align_start']-hrow['hmm_start']
                pfam_motif_end=hrow['align_stop']+(hrow['motifsize']-hrow['hmm_stop'])-1
                if pfam_motif_begin<exp_core_begin and pfam_motif_end<exp_core_end:
                    exp_motif_begin=pfam_motif_begin-addl_motifs[ghfam][motifacc][0]
                    temparr=build_pfam_expmatrix(np.zeros((spgrid[0].shape)),exp_motif_begin,exp_core_end)
                    probdr.loc[:,:,'pfam_prior']+=temparr
                if pfam_motif_begin>exp_core_begin and pfam_motif_end>exp_core_end:
                    exp_motif_end=pfam_motif_end+addl_motifs[ghfam][motifacc][1]
                    temparr=build_pfam_expmatrix(np.zeros((spgrid[0].shape)),exp_core_begin,exp_motif_end)
                    probdr.loc[:,:,'pfam_prior']+=temparr

#    probdr.loc[:,:,'data_likelihood']+=1e-3
    student_dict={}
    if eform=='cc':
        efvpwxra=mds.vpwxra_cc
    elif eform=='full':
        efvpwxra=mds.vpwxra_full

    for cs in mds.curateseq.values:
        if exclude_self and cs==acc:
            continue
        temparr=np.zeros((spgrid[0].shape))
        x=int(round(efvpwxra.loc[acc,cs,'avg_start'].values.item(0),0))
        y=int(round(efvpwxra.loc[acc,cs,'avg_stop'].values.item(0),0))
        ns=efvpwxra.loc[acc,cs,'normscore']
        for xg in np.arange(x-10,x+10,1):
            for yg in np.arange(y-10,y+10,1):#,seqgrid[1].ravel()):
                if xg<0 or yg<0 or xg>=spgrid[0].shape[0] or yg>=spgrid[0].shape[0]:
                    continue
                dist2peak=np.linalg.norm(np.array([xg,yg])-np.array([x,y]))
                if dist2peak in student_dict.keys():
                    temparr[yg,xg]=1e3*student_dict[dist2peak]*1e3**(ns.values.item(0)*1.8-0.1)
                else:
                    studval=student.pdf(dist2peak,5,loc=0,scale=5)
                    student_dict[dist2peak]=studval
                    temparr[yg,xg]=1e3*student_dict[dist2peak]*1e3**(ns-0.1)
        probdr.loc[:,:,'prob_data']+=temparr#1e3*temparr*1e3**(ns-0.1)
    probdr.loc[:,:,'sum_prior']=probdr.loc[:,:,'esize_prior']+probdr.loc[:,:,'pfam_prior']
    
    probdr.loc[:,:,'data_likelihood']=probdr.loc[:,:,'prob_data']*probdr.loc[:,:,'sum_prior']#1e3*temparr*1e3**(ns-0.1)
    probdr.loc[:,:,'posterior']=probdr.loc[:,:,'data_likelihood']/np.sum(probdr.loc[:,:,'data_likelihood'])
    return probdr
#
import skimage
from skimage import morphology

#get number of regions with pmf>0.2
def get_threshold_pmf_count(image,lblimage,rprops):
    thresh_rprops=[]
    total_prob=0.0
    for rprop in rprops:
        cool=np.where(lblimage==rprop.label,image,0)
        sumregionprob=cool.sum()
        #print(sumregionprob,rprop.intensity_image.sum())
        assert(abs(sumregionprob-rprop.intensity_image.sum())<1e-8)
        total_prob+=sumregionprob
        if sumregionprob>0.2:
            thresh_rprops.append(rprop)      
    return total_prob,thresh_rprops

def get_startstop_rois(posteriorgridimg):
    for log_threshprob in np.linspace(-1,-5,50):
        threshprob=np.power(10,log_threshprob)
        markers =np.zeros_like(posteriorgridimg)
        markers[posteriorgridimg>threshprob]=1
        lblimage,lcnt=skimage.morphology.label(markers,connectivity=2,return_num=True)
        rps=skimage.measure.regionprops(lblimage,intensity_image=posteriorgridimg)
        sumprob,thresh_rprops=get_threshold_pmf_count(posteriorgridimg,lblimage,rps)
        if sumprob>0.9:
            break
    return thresh_rprops

def get_start_stop(probxr):
    glength=len(probxr.start)
    posteriorgridimg=probxr.loc[:,:,'posterior'].values.copy()
    rps=get_startstop_rois(posteriorgridimg)
    if len(rps)>0:
        maxprob=max([np.sum(x.intensity_image) for x in rps])
        if maxprob<0.3:
            print('***warning-low max prob...***')
    ssrp=None
    if len(rps)==1:
        ssrp=rps[0]
    elif len(rps)==2:
        rp1_wtdcentroid=rps[0].weighted_centroid
        rp1_wtdlength=rp1_wtdcentroid[0]-rp1_wtdcentroid[1]
        rp2_wtdcentroid=rps[1].weighted_centroid
        rp2_wtdlength=rp2_wtdcentroid[0]-rp2_wtdcentroid[1]
        ssrp=rps[0]
        if abs(rp2_wtdlength-rp1_wtdlength)>40:
            print('selecting based on significant length diff indicating aux module')
            if rp2_wtdlength>rp1_wtdlength:
                ssrp=rps[1]
#                rp2wtdscore=min(3,3*(rp2_wtdlength-rp1_wtdlength)/100)*np.sum(rps[1].intensity_image)
#                if rp2wtdscore>np.sum(rps[0].intensity_image):
#                    print('sufficient score/lngth diff to choose longer')
#                    ssrp=rps[1]
#                else:
#                    ssrp=rps[0]
#            if rp1_wtdlength>rp2_wtdlength:
#                rp1wtdscore=min(3,3*(rp1_wtdlength-rp2_wtdlength)/100)*np.sum(rps[0].intensity_image)
#                if rp1wtdscore>np.sum(rps[1].intensity_image):
#                    print('sufficient score/lngth diff to choose longer')
#                    ssrp=rps[0]
#                else:
#                    ssrp=rps[1]
##                ssrp=rps[1]
        else: 
            print('choosing based on prob intensity')
            if np.sum(rps[1].intensity_image)>np.sum(rps[0].intensity_image):
                ssrp=rps[1]
#            elif np.max(rps[1].intensity_image)>2*np.max(rps[0].intensity_image):
#                ssrp=rps[1]

    elif len(rps)>2:
        print('>2 possibilities. curious...')
        return None
    else:
        print('no probable location')
        return None

    construct_maxlength=ssrp.bbox[2]-ssrp.bbox[1]
    construct_minlength=ssrp.bbox[0]-ssrp.bbox[3]
    lwtdrp_intensity_image=ssrp.intensity_image.copy()#posteriorgrid.copy()

    for row in range(lwtdrp_intensity_image.shape[0]):
        for col in range(lwtdrp_intensity_image.shape[1]):
            orig_rval=row+ssrp.bbox[0]
            orig_cval=col+ssrp.bbox[1] 
            sfactor=((orig_rval-orig_cval)-0.9*construct_minlength)/(construct_maxlength-construct_minlength)
            lwtdrp_intensity_image[row][col]*=sfactor#((orig_rval-orig_cval)-0.9*construct_minlength)/(construct_maxlength-construct_minlength)
    maxlwtd=np.argmax(lwtdrp_intensity_image)
    max_lwtd_loc=np.unravel_index(maxlwtd,lwtdrp_intensity_image.shape)
    bestrow,bestcol=max_lwtd_loc[0]+ssrp.bbox[0],max_lwtd_loc[1]+ssrp.bbox[1]
    if bestcol<=10:
        bestcol=0
    if abs(bestrow-glength)<=16:
        bestrow=glength
    return bestcol,bestrow

def calc_prob_density(probxr):
    stacked=probxr.loc[:,:,'posterior'].stack(coords=('xvals','yvals'))
    return stacked
import matplotlib.gridspec as gridspec
import skimage
def talign_plot(mplt1,mplt2,pwidf):
#from itertools import combinations_with_replacement
    numleaves_t1=len(mplt1.ordered_leaves)
    numleaves_t2=len(mplt2.ordered_leaves)
    inter_array=np.zeros((numleaves_t1,numleaves_t2))
    for a in range(numleaves_t1):
        for b in range(numleaves_t2):
            a_accs=mplt1.ordered_leaves[a].accs
            b_accs=mplt2.ordered_leaves[b].accs
            inter_array[a,b]=np.log(pwidf.loc[a_accs[0],b_accs[0]])
#            subdf=pwidf.loc[a_accs,b_accs]
#            avgval=np.log(subdf.mean().mean())
#            inter_array[a,b]=avgval
    startimg=skimage.img_as_float(inter_array)
#    cool2=skimage.exposure.rescale_intensity(startimg,in_range=(startimg.min(),3.95))
    cool2=skimage.exposure.rescale_intensity(startimg,in_range=(3.2,4.6))
    print(startimg.min(),startimg.max())
    print(cool2.min(),cool2.max())
#    prin
#    cool2=startimg

    fig=plt.figure(figsize=(20,20))
    plot_grid=gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0)
    ax1=plt.Subplot(fig,plot_grid[1])
    ax2=plt.Subplot(fig,plot_grid[2])
    ax3=plt.Subplot(fig,plot_grid[3])
#mpltree.orientation='top'
#mpltree.orientation='left'
    fig.add_subplot(ax1)
    fig.add_subplot(ax3)
    fig.add_subplot(ax2)
    mplt2.render(ax=ax1,orientation='top',tree_lw=1.0)#,create_leaf_names=True,draw_leaf_names=True);

    mplt1.render(ax=ax2,orientation='left',tree_lw=1.0)#,create_leaf_names=True,draw_leaf_names=True);
    ax3.imshow(cool2,aspect='auto')
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.axis('off')
    ax3.spines['left'].set_visible(False)
    ax1.axis('off')
    ax2.axis('off');
    plt.tight_layout()



def bigplotter(accs,ghfam,mds,eform,savefigname=None):
    #for a list of gbacc ids, calculates probgrid, then start_stops, then plot_ebounds
    dbpathstr=f'{ghfam}/{ghfam}DB.sql'
    esize_prior=get_esizeprior(dbpathstr,eform=eform)
    fig,axs=plt.subplots(nrows=-(-len(accs)//3),ncols=3,figsize=(22, -(-len(accs)//3)*6 )    )
    for accnum,acc in enumerate(accs):
        ax=np.ravel(axs)[accnum]
        if eform in ['both','full']:
            ef='full'
        else:
            ef='cc'
        probxr=calc_prob_grid(acc,ghfam,esize_prior,dbpathstr,mds,eform=ef,exclude_self=False)
        startstop=get_start_stop(probxr)
        dbseqidx=list(mds.dbseq.values).index(acc)
        plot_ebounds(mds,dbpathstr,dbseqidx,ax=ax) 
        if startstop is None:
            startstop=['nd','nd']
        sortedpdbscores=mds.vpwxra_cc.loc[acc,:,'normscore'].sortby(mds.vpwxra_cc.loc[acc,:,'normscore'],ascending=False)
        ax.text(0.55,0.9,f'{startstop[0]} - {startstop[1]}',size=22,transform=ax.transAxes)
        for score_rank in range(3):
            ax.text(0.45,0.8-score_rank*0.06,
                f'{sortedpdbscores.curateseq[score_rank].values.item(0)}: {sortedpdbscores[score_rank].values.item(0):.2f}',
                size=20,transform=ax.transAxes)
        if savefigname is not None:
            plt.savefig(savefigname,dpi=100)

def get_avg_exts(ghfam,dbpath):
    motifselection=core_motifs[ghfam][0]
    cc_sizes=[]
    n_exts=[]
    c_exts=[]
    conn=sqlite3.connect(dbpath)
    conn.row_factory=sqlite3.Row
    c=conn.cursor()
    c.execute('''SELECT * FROM CURATEDXTALS''')
    rows=c.fetchall()
    for row in rows:
        acc=row['acc']
        start=row['pgbsr_ccstart']
        stop=row['pgbsr_ccstop']
        cc_sizes.append(stop-start+1)
        c.execute('''SELECT * FROM HMMERSEQDATA WHERE acc=(?)''',(acc,))
        hrows=c.fetchall()
        for hrow in hrows:
            if hrow['motifacc']==motifselection:
                halstart=hrow['align_start']-1
                halstop=hrow['align_stop']-2
                if halstart>=start and halstop<=stop:
                    pfam_begin=hrow['align_start']-hrow['hmm_start']
                    pfam_end=hrow['align_stop']+(hrow['motifsize']-hrow['hmm_stop'])-1
                    n_exts.append(pfam_begin-start)
                    c_exts.append(stop-pfam_end)    
    conn.close()
    return np.median(cc_sizes),np.median(n_exts),np.median(c_exts)