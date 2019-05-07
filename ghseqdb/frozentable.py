import sqlite3,pickle#,re,datetime,pickle
from Bio import SeqIO
from . import seqdbutils

def build_frozentable(dbpathstr,seqfpath,notes=None,seqformat='fasta'):
    srsdict=SeqIO.to_dict(SeqIO.parse(seqfpath,seqformat))
    assert (len(srsdict.keys())>0), "no seqs"
    conn=seqdbutils.gracefuldbopen(dbpathstr) 
    c=conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS FROZENSEQS (acc text, pklfsr glob, notes text)''')
    c.execute('''SELECT acc FROM FROZENSEQS''')
    cur_fraccs=[x['acc'] for x in c.fetchall()]
    fraccs2add=list(set(srsdict.keys()).difference(cur_fraccs))
    print(f'attempting to add {len(fraccs2add)} of {len(srsdict.keys())} seqs in seqfile') 
    for fracc in fraccs2add:
        pklfsr=pickle.dumps(srsdict[fracc])
        frnotes=None
        if type(notes)==str:
            frnotes=notes
        newtuple=(fracc,pklfsr,frnotes)
        c.execute('''INSERT INTO FROZENSEQS VALUES (?,?,?)''',newtuple)
    conn.commit()
    conn.close()