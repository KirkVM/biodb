import os
import pandas as pd
from . import seqdbutils


def build_curatedxtaltable(dbpath,xtalxlfpath):
    xtaldf=pd.read_excel(xtalfpath,....)
    xtalfile_mtime=os.stat(xtalxlfpath).st_mtime

    conn=seqdbutils.gracefuldbopen(dbpath)
    c.execute('''CREATE TABLE IF NOT EXISTS CURATEDXTALS (acc text, chain text, xtalfname text,\
                xtalf_mtime int, ntccpos int, ctccpos int, ntccseq text, ctccseq text)''')
    pass