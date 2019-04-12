import pickle
from pathlib import Path
from ete3 import Tree
import pandas as pd
from . import seqdbutils
def curatedrefine(dbpathstr,tfpathstr,pkl_intpwiddfpathstr,acclist=['AAA23224'],edatadbloc='up1'):
    pwid_df=pd.read_pickle(pkl_intpwiddfpathstr)
    conn=seqdbutils.gracefuldbopen(dbpathstr)
    seqdbutils.check_tables_exist(conn,['CURATEDXTALS','PROTEINGBS','CCDATA'])
    c=conn.cursor()
    #c.execute('''SELECT CURATEDXTALS.acc,PROTEINGBS.pklgbsr FROM CURATEDXTALS INNER JOIN PROTEINGBS WHERE CURATEDXTALS.acc=PROTEINGBS.acc''')
    c.execute('''SELECT acc FROM CURATEDXTALS''')
    rows=c.fetchall()
    curatedaccs=[x['acc'] for x in rows]
    removedaccs=list(set(curatedaccs).difference(pwid_df.columns))
    for ra in removedaccs:
        c.execute('''SELECT * FROM CCDATA WHERE acc=(?)''',(ra,))
        row=c.fetchone()
        print(ra)
        for k in row.keys():
            print(k,row[k])
    curatedaccs=list(set(curatedaccs).intersection(pwid_df.columns))
    tree=Tree(tfpathstr)
    tree.prune(curatedaccs)
#    print(tree)
#    taccs=list(   set(accs).intersection(tree.get_leaf_names()) )
#    tree.prune(taccs)
    for acc in []:# ['AAA23231']:
        curaccs=curatedaccs[:]
        curaccs.append(acc)
        curaccs=list(set(curaccs).intersection(pwid_df.columns))
        myintdf=pwid_df.loc[curaccs,curaccs]
        print(myintdf.shape)
        top_comp_value=0
        top_comp=None
        for curatedacc in curatedaccs:
            d2acc=myintdf[acc].loc[curatedacc]
            if d2acc>top_comp_value:
                top_comp=curatedacc
                top_comp_value=d2acc

        c.execute('''SELECT * FROM CCDATA WHERE acc=(?)''',(acc,))
        ccrow=c.fetchone()
        c.execute('''SELECT * FROM CURATEDXTALS WHERE acc=(?)''',(top_comp,))
        cxrow=c.fetchone()
        c.execute('''SELECT * FROM PROTEINGBS WHERE acc=(?)''',(top_comp,))
        pgbrow=c.fetchone()
        curated_ntccpos=cxrow['ntccpos']
        curated_ntccseq=cxrow['ntccseq']

        curated_ctccpos=cxrow['ctccpos']
        curated_ctccseq=cxrow['ctccseq']
        #now we see if this matches the sequence-
        pgbsr=pickle.loads(pgbrow['pklgbsr'])
        print(curated_ntccpos,curated_ctccpos)
        print(curated_ntccseq,pgbsr.seq[max(0,curated_ntccpos-10):curated_ntccpos+10])
        print(curated_ctccseq,pgbsr.seq[curated_ctccpos-10:min(len(pgbsr.seq),curated_ctccpos+10)])
        print(str(pgbsr.seq))
    tree.set_outgroup(top_comp)
#    print(tree)
#        curated_ntccoffset=cxrow['']

#            print(curatedacc,d2acc)
#    print(tree)
#    print(accs)
#    print(versions)
      