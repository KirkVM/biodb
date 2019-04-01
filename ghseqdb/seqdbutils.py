#from pathlib import Path
#import os
import sqlite3,atexit

def gracefuldbopen(dbpath):
    ''' returns a connection with after doing some boilerplate'''
        #make path if needed....
        #if os.path.exists...
    try:
        conn=sqlite3.connect(dbpath)
        atexit.register(conn.close)
        conn.row_factory=sqlite3.Row
        return conn
    except:
        raise Exception("could not connect to db")


def table_exists(tblname,conn,fieldnames=None):
    try:
        cursor=conn.execute('''FROM (?) SELECT *''',(tblname,))
#        cursor=conn.execute('''SELECT * FROM CAZYSEQDATA''')
        colnames=[desc[0] for desc in cursor.description]
        print(colnames)
        if fieldnames is not None:
            if len(set(colnames).union(fieldnames))!=len(fieldnames):
                print('table exists but num columns is different than expected...')
        return True
    except:
        print('table does not exist yet')
        return False
#        table does not exist yet
 

#    try:
#        c.execute('''SELECT COUNT (*) FROM CAZYDATA''')
#        cazydbsize=c.fetchone()[0]
#        print(f'adding/updating existing CAZYDATA table, size {cazydbsize}')
#    except:
#        #table does not exist yet
#        c.execute('''CREATE TABLE CAZYDATA (acc text, version text, scrapedate text, \
#            subfam text, extragbs text, ecs text, pdbids text, uniprotids text)''')
 