import operator,pdb
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class dMatch:
        def __init__(self,theline):
                entries=list(map(lambda x:x.strip(),theline.split()))
                self.dname=entries[0].split('.')[0]
                self.dsize=int(entries[2])
                self.accession=entries[3]
                self.psize=int(entries[5])
                self.fevalue=float(entries[6])
                self.fscore=float(entries[7])
                self.fbias=float(entries[8])
                self.matchnum=int(entries[9])
                self.matchtot=int(entries[10])
                self.dcevalue=float(entries[11])
                self.dievalue=float(entries[12])
                self.dscore=float(entries[13])
                self.dbias=float(entries[14])
                self.hmmstart=int(entries[15])
                self.hmmstop=int(entries[16])
                self.alstart=int(entries[17])
                self.alstop=int(entries[18])
                self.envstart=int(entries[19])
                self.envstop=int(entries[20])
                self.acc=float(entries[21])
                self.picdstart=None
                self.picdstop=None
                self.modded_start=False #has picdstart position been modified to accommodate adjacent domain?
                self.modded_stop=False #has picdstop position been modified to accommodate adjacent domain?
                


def parse_hmmfile(domtbloutfpath):
        """returns a dictionary with a list for each domain annoation:
        [d_name,hmm_domain_start,hmm_domain_end,protein_domain_start,protein_domain_end,hmm_domain_size,coverage,cvalue]
        lists are sorted according to protein_domain_start"""
        dbcanfile=open(domtbloutfpath,'r')
        dHT={}
        for line in dbcanfile.readlines():
                if line[0]=="#":continue
                hmmermatch=dMatch(line)
#                if not dHT.has_key(hmmermatch.accession):dHT[hmmermatch.accession]=[]
                if not hmmermatch.accession in dHT.keys():dHT[hmmermatch.accession]=[]
                dHT[hmmermatch.accession].append(hmmermatch)
        #now sort each by dievalue:
        for acc in dHT:
                dHT[acc]=sorted(dHT[acc],key=lambda x:x.dievalue)
        return dHT


def cleanannotations(dmatches_):
        topmatches_=[]
        for dmatch in dmatches_:
                dmatch.picdstart=max(1,dmatch.alstart-dmatch.hmmstart+1)
                dmatch.picdstop=min(dmatch.psize,dmatch.alstop+(dmatch.dsize-dmatch.hmmstop))

                addable=True
                modded=False 
                tmindex=0
                while addable and tmindex<len(topmatches_):
                        topmatch=topmatches_[tmindex]
                        #for readability:
                        #if topmatch.accession=='AIQ32060.1':pdb.set_trace()
                        curname=dmatch.dname
                        curstart=dmatch.picdstart
                        curstop=dmatch.picdstop
                        prevname=topmatch.dname
                        prevstart=topmatch.picdstart
                        prevstop=topmatch.picdstop
                
                        mergestatus=None
                        if (curstart>=prevstart and curstop<prevstop) or (curstart>prevstart and curstop<=prevstop):
                                mergestatus="within"
                        elif curstart<=prevstart and curstop>=prevstop:
                                mergestatus="span"
                        elif curstart<=prevstop and curstop>prevstop:
                                mergestatus="ohright"
                        elif curstart<prevstart and curstop>=prevstart:
                                mergestatus="ohleft"
                        elif curstart>prevstop or curstop<prevstart:
                                mergestatus="clean"
                        else:
                                exit("this shouldn't happen!")
        
                        if mergestatus in ['within','span']:
                                addable=False
                        elif mergestatus=="ohright":
                                overlapsize=prevstop-curstart+1
                                freesize=curstop-prevstop
                                #trying to add CBM46 to right of GH5
                                if curname=="CBM46" and prevname=="GH5":
                                        if overlapsize<20:
                                                dmatch.picdstart=prevstop+1
                                                dmatch.modded_start=True
                                        else:
                                                addable=False
                                #trying to add CBM to right of GH
                                elif curname[0:3]=="CBM" and prevname[0:2]=="GH":
                                        if overlapsize<10:
                                                dmatch.picdstart=prevstop+1
                                                dmatch.modded_start=True
                                        else:
                                                addable=False
                                #trying to add GH to right of CBM
                                elif curname[0:2]=="GH" and prevname[0:3]=="CBM":
                                        if overlapsize<10:
                                                topmatch.picdstop=curstart-1
                                                topmatch.modded_stop=True
                                        else:
                                                addable=False
                        elif mergestatus=="ohleft":
                                overlapsize=curstop-prevstart+1
                                freesize=prevstart-curstart
                                #trying to add GH5 to left of CBM46
                                if curname=="GH5" and prevname=="CBM46":
                                        if overlapsize<20:
                                                topmatch.picdstart=curstop+1
                                                topmatch.modded_start=True
                                        else:
                                                addable=False
                                #trying to add GH to left of CBM
                                elif curname[0:2]=="GH" and prevname[0:3]=="CBM":
                                        if overlapsize<10:
                                                topmatch.picdstart=curstop+1
                                                topmatch.modded_start=True
                                        else:
                                                addable=False
                                #trying to add CBM to left of GH
                                elif curname[0:3]=="CBM" and prevname[0:2]=="GH":
                                        if overlapsize<10:
                                                dmatch.picdstop=prevstart-1
                                                dmatch.modded_stop=True
                                        else:
                                                addable=False
                        
                        #now, if overlap is <4 or <4% for both domains AND not already modded
                        #AND good e-values,truncate the larger domain and allow
                        if (mergestatus=='ohleft' and dmatch.modded_stop==False and topmatch.modded_start==False) or\
                                (mergestatus=='ohright' and dmatch.modded_start==False and topmatch.modded_stop==False):
                                if topmatch.dievalue<0.01 and dmatch.dievalue<0.01:
                                        #sort the domains so longer is first
                                        doms_=sorted([topmatch,dmatch],key=lambda x:x.dsize,reverse=True)
                                        if (overlapsize<4) or (float(overlapsize)/doms_[1].dsize<0.04):
                                                if doms_[0].picdstart<doms_[1].picdstart: #long one is first
                                                        doms_[0].picdstop=doms_[1].picdstart-1
                                                        doms_[0].modded_stop=True
                                                else:
                                                        doms_[0].picdstart=doms_[1].picdstop+1
                                                        doms_[0].modded_start=True
                                        else:
                                                addable=False
                                else:
                                        addable=False
                        tmindex+=1
                        #now see if necessary truncations make it implausible:
                        if (float(dmatch.picdstop-dmatch.picdstart))/(float(dmatch.dsize))<0.5:
                                addable=False
                if addable:
                        topmatches_.append(dmatch)
        return topmatches_


#returns a list of domains for each accession
def parse_dbcan(domtbloutfpath):
        domHT=parse_hmmfile(domtbloutfpath)
        for acc in domHT:
                domHT[acc]=cleanannotations(domHT[acc])
        return domHT




#this is written for dbcankools_old....not sure if it's functional?
#def getcbm46seqs(domtbloutfpath,fullseqfpath):
#       domdict=parse_hmmfile(domtbloutfpath)
#       for proteinacc in domdict.keys():
#               domerges(domdict[proteinacc])
#               killsmalls(domdict[proteinacc])
#               #killoverlaps(domdict[proteinacc])
#       cbm46count=0
#       fullseqrecs=list(SeqIO.parse(fullseqfpath,'fasta'))
#       fullseqids=map(lambda x:x.id,fullseqrecs)
#       nterm46s=[]
#       cterm46s=[]
#       for proteinacc in domdict.keys():
#               dentries=domdict[proteinacc]
#               cbm46entries=filter(lambda x:x[0]=='CBM46',dentries)
#               numcbm46=len(cbm46entries)
##              if numcbm46!=2 and numcbm46>0:
##                      print proteinacc
#               if numcbm46!=2:
#                       #currently this is 15 sequences. Check into these
#                       continue
#               #get the seq:
#               curprotseq=fullseqrecs[fullseqids.index(proteinacc)].seq
#               for x,cbm46 in enumerate(cbm46entries):
#                       hmm_start=int(cbm46[1])
#                       hmm_stop=int(cbm46[2])
#                       prot_start=int(cbm46[3])
#                       prot_stop=int(cbm46[4])
#                       fulldom_pstart=prot_start-(hmm_start-1)
#                       fulldom_pstop=prot_stop+(87-hmm_stop)
#                       cursubseq=curprotseq[fulldom_pstart:fulldom_pstop]
#                       if x==0:
#                               nterm46s.append(SeqRecord(cursubseq,id=proteinacc))
#                       elif x==1:
#                               cterm46s.append(SeqRecord(cursubseq,id=proteinacc))
#       SeqIO.write(nterm46s,"cbm46nt.fasta","fasta")
#       SeqIO.write(cterm46s,"cbm46ct.fasta","fasta")
#getcbm46seqs("superfreshall.txt","newseqannotations_81516/fullcodedprotein4_allplus.fasta")
