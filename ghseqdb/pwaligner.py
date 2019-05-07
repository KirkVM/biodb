#import pwaligner
import pickle,sys,tarfile,re,os,dask,copyreg
import xarray as xr
from pathlib import Path
from Bio import SeqIO,pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import numpy as np

class VPWAInfo:
    def __init__(self,dbacc,curateacc,bpalns,gap_penalty,ext_penalty,penalize_end_gaps,chtc_mode=False):
        ''''sets up aln info

        Arguments:
            dbacc- accession code of sequence to be investigated
            curateacc- acc code of curated sequence (length truncd to known domain(s) size)
            bpalns: a biopython alignment object (list of alignments-
                    -alignment ordering should be performed according to pw(dbseq,curateseq) (can I verify this??)
        '''
        self.dbacc=dbacc
        self.curateacc=curateacc
        self.alns=bpalns #a biopython pw aligns result
        self.gap_penalty=gap_penalty
        self.ext_penalty=ext_penalty
        self.penalize_end_gaps=penalize_end_gaps
        self.dbseqlen=None
        self.curateseqlen=None
        self.seqaln_density=None
        self.vseqaln_curatedensity=None
        self.score=None
        self.starts=[]
        self.stops=[]
        self.avg_start=None
        self.avg_stop=None
        self.median_start=None
        self.median_stop=None
        self.initialize_info()
        if len(self.alns)>0:
            self.buildinfo()
        if chtc_mode:
            del self.alns
    def initialize_info(self):
        if len(self.alns)==0:
            self.score=np.nan
        else: #use aln[0] to set up
            aln=self.alns[0]
            self.score=aln[2]
            self.dbseqlen=len(aln[0])-aln[0].count('-')
            self.curateseqlen=len(aln[1])-aln[1].count('-')
            self.seqaln_density=np.zeros((self.dbseqlen))
            self.vseqaln_curatedensity=np.zeros((self.dbseqlen))
            self.avg_start=0.0
            self.avg_stop=0.0
    def buildinfo(self):
        for aln in self.alns:
            curateseq_ended=False
            dbalnseq=aln[0]
            curatealnseq=aln[1]
            db_pos=-1
            curate_pos=-1
            aln_pos=-1
            aln_length=len(aln[0])
            for aln_pos in range(aln_length):
                dbres=dbalnseq[aln_pos]
                curateres=curatealnseq[aln_pos]
                if dbres!='-':
                    db_pos+=1
                    aln_nbrs=dbalnseq[max(0,aln_pos-3):min(aln_length,aln_pos+4)]
                    aln_nbrsize=len(aln_nbrs)-1 #don't include current position
                    resdensity=(aln_nbrsize-aln_nbrs.count('-'))/aln_nbrsize
                    self.seqaln_density[db_pos]+=resdensity
                    
                    valn_nbrs=curatealnseq[max(0,aln_pos-3):min(aln_length,aln_pos+4)]
                    valn_nbrsize=len(valn_nbrs) #do include current position
                    vresdensity=(valn_nbrsize-valn_nbrs.count('-'))/valn_nbrsize
                    self.vseqaln_curatedensity[db_pos]+=vresdensity
                if curateres!='-':
                    curate_pos+=1
                if curate_pos==0:
                    self.starts.append(db_pos)
                if curate_pos==self.curateseqlen-1 and not curateseq_ended:
                    self.stops.append(db_pos)
                    curateseq_ended=True
        self.avg_start=np.mean(self.starts)
        self.avg_stop=np.mean(self.stops)
        self.median_start=np.median(self.starts)
        self.median_stop=np.median(self.stops)
        self.seqaln_density/=len(self.alns)
        self.vseqaln_curatedensity/=len(self.alns)

def biopython_vpwaln(dbsr,curatesr,gap_penalty=-10,ext_penalty=-0.5,penalize_end_gaps=(True,False),chtc_mode=False):
    alns=pairwise2.align.globalds(dbsr.seq,curatesr.seq,blosum62,gap_penalty,ext_penalty,penalize_end_gaps=penalize_end_gaps)
    vpwainfo=VPWAInfo(dbsr.id,curatesr.id,alns,gap_penalty=gap_penalty,ext_penalty=ext_penalty,penalize_end_gaps=penalize_end_gaps,chtc_mode=chtc_mode)
    return vpwainfo

def get_curate_alstats(dbgbsr,cgbsrs,gap_penalty=-10,ext_penalty=-0.5,penalize_end_gaps=(True,False),chtc_mode=False):
    vpwainfos=[]
    for cgbsr in cgbsrs:
        curvpwainfo=biopython_vpwaln(dbgbsr,cgbsr,gap_penalty=gap_penalty,ext_penalty=ext_penalty,penalize_end_gaps=penalize_end_gaps,chtc_mode=chtc_mode)
        vpwainfos.append(curvpwainfo)
    return vpwainfos

def multi_curate_alstats(dbgbsrs,cgbsrs,gap_penalty=-10,ext_penalty=-0.5,penalize_end_gaps=(True,False),chtc_mode=False):
    allvpwainfos=[]
    for dbgbsr in dbgbsrs:
        vpwainfos=get_curate_alstats(dbgbsr,cgbsrs,gap_penalty=gap_penalty,ext_penalty=ext_penalty,penalize_end_gaps=penalize_end_gaps,chtc_mode=chtc_mode)
        allvpwainfos.append(vpwainfos)
    return allvpwainfos
    

####### HERE IS STUFF TO GET PW ALIGNMENT FROM DIRECTORY OF *.tgz FILES###########
def build_xarray_fromvpwais(vpwais):
    #for pwai in pwais:
    data_array=np.zeros((len(vpwais),len(vpwais[0]),10))
    dbaccs=[x[0].dbacc for x in vpwais]
    caccs=[x.curateacc for x in vpwais[0]]
    for dbidx,dbpw in enumerate(vpwais):
        for cidx,vpwai in enumerate(dbpw):
            #print(pwai.median_start)
            data_array[dbidx][cidx][0]=vpwai.score
            data_array[dbidx][cidx][1]=np.nan
            data_array[dbidx][cidx][2]=np.nan
            data_array[dbidx][cidx][3]=np.nan
            data_array[dbidx][cidx][4]=vpwai.median_start
            data_array[dbidx][cidx][5]=vpwai.avg_start
            data_array[dbidx][cidx][6]=vpwai.median_stop
            data_array[dbidx][cidx][7]=vpwai.avg_stop
            data_array[dbidx][cidx][8]=np.nan
            data_array[dbidx][cidx][9]=np.nan
    params=['score','normscore','lognormscore','log10normscore','mdn_start','avg_start','mdn_stop','avg_stop','start_conf','stop_conf']
    vpwdr=xr.DataArray(data_array,coords=[dbaccs,caccs,params],dims=['dbseq','curateseq','parameter'])
    return vpwdr

class VPWAInfoUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "__main__":
            return VPWAInfo
        return super().find_class(module, name)

def xra_from_file(fpath,tarball=True,pickled=True):
    if tarball:
        tf=tarfile.open(fpath,mode="r:gz")
        ti=tf.getmembers()[0]
        vpwafile=tf.extractfile(ti)
        if pickled: #not sure why this would not be case
            #vpwais=pickle.load(vpwafile)
            Unpickler=VPWAInfoUnpickler(vpwafile)
            vpwais=Unpickler.load() 
        tf.close()
    return build_xarray_fromvpwais(vpwais)

def build_vpwxra_from_directory(pathstr,tarballs=True,pickled=True):    
    fpaths=[]
    if tarballs:
        tRE=re.compile('.+\.tgz|.+\.tar\.gz')
        tfnames=[x for x in os.listdir(pathstr) if tRE.match(x) is not None]
        for tfname in tfnames:#[:5]:
            tfpath=Path(pathstr) / tfname
            fpaths.append(tfpath.resolve())
    else:
        print('tarball method only one currently implemented')
        return None
    xras=[]
    appendtasks=[]
    for fpath in fpaths:
        appendtasks.append( dask.delayed( xra_from_file,pure=True)(fpath))
    xras=dask.compute(*appendtasks)
    vpwxra=xr.concat(xras,dim='dbseq')
    return vpwxra
####### HERE IS STUFF TO GET PW ALIGNMENT FROM DIRECTORY OF *.tgz FILES###########

def do_chtc_function(argv=None):
    if argv is None:
        argv=sys.argv
    assert (len(argv)==5), "require 4 args: dbsrfile,cursrfile,numsplits,splitnum"
    dbsrfile=argv[1]
    cursrfile=argv[2]
    numsplits=int(argv[3])
    splitnum=int(argv[4])
    dbsrs=list(SeqIO.parse(dbsrfile,'fasta'))
    cursrs=list(SeqIO.parse(cursrfile,'fasta'))
    idxsplits=np.array_split(range(len(dbsrs)),numsplits)
    split_dbsrs=dbsrs[idxsplits[splitnum][0]:idxsplits[splitnum][-1]+1]
    allvpwainfos=multi_curate_alstats(split_dbsrs,cursrs,chtc_mode=True)
    
    #copyreg.pickle(VPWAInfo, pickle_VPWAInfo)
    with open(f'vpwaset{splitnum}.pkl','wb') as f:
        pickle.dump(allvpwainfos,f)

if __name__=="__main__":
    do_chtc_function()
