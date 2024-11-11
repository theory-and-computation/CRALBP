import sqlite3
import uuid

def readFile(_filename):
    _entries = []
    with open(_filename, 'r') as fileobj:
        lines = fileobj.readlines()

        try:
            lines.remove('\n')
        except:
            pass

        for line in lines:
            split_line = line.strip().split(',')
            _entries.append((split_line[0], split_line[1], split_line[2], split_line[3], split_line[4]))
        return _entries


def addAllEntries(_allentries):
    for row in _allentries:
        addEntry(int(row[0]), row[1], row[2], int(row[3]))



def createTable():
    conn = sqlite3.connect('RSCM_contacts.db')
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS protein (
        resid INTEGER, 
        restype TEXT,
        lipid_restypes NOT NULL, 
        frame INTEGER)
    '''
    )
    conn.commit()
    conn.close()

# Function to add a new entry to the database
def addEntry(_resid: int, _restype: str, _lipid_restypes: tuple, _frame: int):
    conn = sqlite3.connect('RSCM_contacts.db')
    c = conn.cursor()
    c.execute("INSERT INTO protein (resid, restype, lipid_restypes, frame) VALUES (?, ?, ?, ?)", (_resid, _restype, _lipid_restypes, _frame))
    conn.commit()
    conn.close()


if __name__ == '__main__':

    #createTable()
    entries = readFile("./RSCM_contacts.txt")
    addAllEntries(entries)

