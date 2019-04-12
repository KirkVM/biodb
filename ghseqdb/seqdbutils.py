#import os
from pathlib import Path
import sqlite3,atexit,sys

def gracefuldbopen(dbpathstr):
    ''' returns a connection with after doing some boilerplate'''
        #make path if needed....
        #if os.path.exists...
    dbpath=Path(dbpathstr)
    if not dbpath.exists():
        if not dbpath.parent.exists():
            sys.exit(f'target folder {dbpath.parent} does not exist. Is dbpath correct?')
        else:
            print(f'creating new database {dbpath.name} at {dbpath.parent}')
    try:
        conn=sqlite3.connect(dbpath)
        atexit.register(conn.close)
        conn.row_factory=sqlite3.Row
        return conn
    except:
        raise Exception(f"could not connect to {dbpath.name}")

def check_tables_exist(conn,tblnames):
    c=conn.cursor()
    for tblname in tblnames:
        #import pdb;pdb.set_trace()
        if tblname=='CCDATA':
            c.execute('''SELECT COUNT (*) FROM CCDATA''')
        elif tblname=='PROTEINGBS':
            c.execute('''SELECT COUNT (*) FROM PROTEINGBS''')
        elif tblname=='CURATEDXTALS':
            c.execute('''SELECT COUNT (*) FROM CURATEDXTALS''')
        elif tblname=='CAZYSEQDATA':
            c.execute('''SELECT COUNT (*) FROM CAZYSEQDATA''')
        elif tblname=='XTALS':
            c.execute('''SELECT COUNT (*) FROM XTALS''')
        numrows=c.fetchone()[0]
        if numrows<1:
            conn.close()
            sys.exit(f'no table named {tblname}')







#this is something to run to periodically check status of dbase
#should write to a per-dbfolder logfile
def dbstatus_report(dbpath,dblogfpath):
    pass