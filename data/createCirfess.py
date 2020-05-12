#! /usr/bin/python3

import sqlite3
from sqlite3 import Error

createPepSQL =  '''CREATE TABLE pep (
    ID TEXT PRIMARY KEY,
    pepSeq TEXT,
    Accession TEXT,
    range TEXT,
    numMissedCleavages INTEGER,
    OKforMS TEXT,
    numMotifsNXS INTEGER,
    numMotifsNXT INTEGER,
    numMotifsNXC INTEGER,
    numMotifsNXV INTEGER,
    numMotifsC INTEGER,
    numMotifsK INTEGER,
    motifLocNXS TEXT,
    motifLocNXT TEXT,
    motifLocNXC TEXT,
    motifLocNXV TEXT,
    motifLocC TEXT,
    motifLocK TEXT,
    numMotifsPhobiusNXS INTEGER,
    numMotifsPhobiusNXT INTEGER,
    numMotifsPhobiusNXC INTEGER,
    numMotifsPhobiusNXV INTEGER,
    numMotifsPhobiusC INTEGER,
    numMotifsPhobiusK INTEGER,
    motifLocPhobiusNXS TEXT,
    motifLocPhobiusNXT TEXT,
    motifLocPhobiusNXC TEXT,
    motifLocPhobiusNXV TEXT,
    motifLocPhobiusC TEXT,
    motifLocPhobiusK TEXT,
    numMotifsTMHMMNXS INTEGER,
    numMotifsTMHMMNXT INTEGER,
    numMotifsTMHMMNXC INTEGER,
    numMotifsTMHMMNXV INTEGER,
    numMotifsTMHMMC INTEGER,
    numMotifsTMHMMK INTEGER,
    motifLocTMHMMNXS TEXT,
    motifLocTMHMMNXT TEXT,
    motifLocTMHMMNXC TEXT,
    motifLocTMHMMNXV TEXT,
    motifLocTMHMMC TEXT,
    motifLocTMHMMK TEXT
);'''

createProtSQL =  '''CREATE TABLE prot (
    Accession TEXT PRIMARY KEY,
    Seq TEXT,
    Length INTEGER,
    StringOutPhobius TEXT,
    numTMPhobius INTEGER,
    numICPhobius INTEGER,
    numECPhobius INTEGER,
    StringOutTMHMM TEXT,
    numTMTMHMM INTEGER,
    numICTMHMM INTEGER,
    numECTMHMM INTEGER,
    SigPepPhobius TEXT,
    ScorePhobius TEXT,
    SigPepSignalP TEXT,
    ScoreSignalP TEXT,
    SigPepPrediSI TEXT,
    ScorePrediSI TEXT,
    numSPpredictions INTEGER,
    SPC INTEGER,
    SPCstring TEXT,
    numMotifsNXS INTEGER,
    motifsLocNXS TEXT,
    numMotifsNXT INTEGER,
    motifsLocNXT TEXT,
    numMotifsNXC INTEGER,
    motifsLocNXC TEXT,
    numMotifsNXV INTEGER,
    motifsLocNXV TEXT,
    numMotifsC INTEGER,
    motifsLocC TEXT,
    numMotifsK INTEGER,
    motifsLocK TEXT,
    numMotifsPhobiusNXS INTEGER,
    motifsLocPhobiusNXS TEXT,
    numMotifsPhobiusNXT INTEGER,
    motifsLocPhobiusNXT TEXT,
    numMotifsPhobiusNXC INTEGER,
    motifsLocPhobiusNXC TEXT,
    numMotifsPhobiusNXV INTEGER,
    motifsLocPhobiusNXV TEXT,
    numMotifsPhobiusC INTEGER,
    motifsLocPhobiusC TEXT,
    numMotifsPhobiusK INTEGER,
    motifsLocPhobiusK TEXT,
    numMotifsTMHMMNXS INTEGER,
    motifsLocTMHMMNXS TEXT,
    numMotifsTMHMMNXT INTEGER,
    motifsLocTMHMMNXT TEXT,
    numMotifsTMHMMNXC INTEGER,
    motifsLocTMHMMNXC TEXT,
    numMotifsTMHMMNXV INTEGER,
    motifsLocTMHMMNXV TEXT,
    numMotifsTMHMMC INTEGER,
    motifsLocTMHMMC TEXT,
    numMotifsTMHMMK INTEGER,
    motifsLocTMHMMK TEXT,
    numPepMC0 INTEGER,
    numPepMC1 INTEGER,
    numPepMC2 INTEGER,
    numPepOKMC0 INTEGER,
    numPepOKMC1 INTEGER,
    numPepOKMC2 INTEGER,
    numPepMC0NXS INTEGER,
    numPepMC0NXT INTEGER,
    numPepMC0NXC INTEGER,
    numPepMC0NXV INTEGER,
    numPepMC0C INTEGER,
    numPepMC0K INTEGER,
    numPepMC1NXS INTEGER,
    numPepMC1NXT INTEGER,
    numPepMC1NXC INTEGER,
    numPepMC1NXV INTEGER,
    numPepMC1C INTEGER,
    numPepMC1K INTEGER,
    numPepMC2NXS INTEGER,
    numPepMC2NXT INTEGER,
    numPepMC2NXC INTEGER,
    numPepMC2NXV INTEGER,
    numPepMC2C INTEGER,
    numPepMC2K INTEGER,
    numPepMC0PhobiusNXS INTEGER,
    numPepMC0PhobiusNXT INTEGER,
    numPepMC0PhobiusNXC INTEGER,
    numPepMC0PhobiusNXV INTEGER,
    numPepMC0PhobiusC INTEGER,
    numPepMC0PhobiusK INTEGER,
    numPepMC1PhobiusNXS INTEGER,
    numPepMC1PhobiusNXT INTEGER,
    numPepMC1PhobiusNXC INTEGER,
    numPepMC1PhobiusNXV INTEGER,
    numPepMC1PhobiusC INTEGER,
    numPepMC1PhobiusK INTEGER,
    numPepMC2PhobiusNXS INTEGER,
    numPepMC2PhobiusNXT INTEGER,
    numPepMC2PhobiusNXC INTEGER,
    numPepMC2PhobiusNXV INTEGER,
    numPepMC2PhobiusC INTEGER,
    numPepMC2PhobiusK INTEGER,
    numPepMC0TMHMMNXS INTEGER,
    numPepMC0TMHMMNXT INTEGER,
    numPepMC0TMHMMNXC INTEGER,
    numPepMC0TMHMMNXV INTEGER,
    numPepMC0TMHMMC INTEGER,
    numPepMC0TMHMMK INTEGER,
    numPepMC1TMHMMNXS INTEGER,
    numPepMC1TMHMMNXT INTEGER,
    numPepMC1TMHMMNXC INTEGER,
    numPepMC1TMHMMNXV INTEGER,
    numPepMC1TMHMMC INTEGER,
    numPepMC1TMHMMK INTEGER,
    numPepMC2TMHMMNXS INTEGER,
    numPepMC2TMHMMNXT INTEGER,
    numPepMC2TMHMMNXC INTEGER,
    numPepMC2TMHMMNXV INTEGER,
    numPepMC2TMHMMC INTEGER,
    numPepMC2TMHMMK INTEGER
);'''

# createIndex1 = 'CREATE INDEX pep_acc_idx ON pep (Accession);'
# createIndex2 = 'CREATE INDEX pep_nmc_idx ON pep (numMissedCleavages);'
# createIndex3 = 'CREATE INDEX pep_ok_idx ON pep (OKforMS);'
# createIndex4 = 'CREATE INDEX prot_nspp_idx ON prot ( numSPpredictions );'
# createIndex5 = 'CREATE INDEX prot_npOK0_idx ON prot ( numPepOKMC0 );'
# createIndex6 = 'CREATE INDEX prot_npOK1_idx ON prot ( numPepOKMC1 );'
# createIndex7 = 'CREATE INDEX prot_npOK2_idx ON prot ( numPepOKMC2 );'

createIndex1 = 'CREATE INDEX prot_Acc_idx on prot(Accession);'
createIndex2 = 'CREATE INDEX prot_npOK0_idx ON prot ( numPepOKMC0 );'
createIndex3 = 'CREATE INDEX prot_npOK1_idx ON prot ( numPepOKMC1 );'
createIndex4 = 'CREATE INDEX prot_npOK2_idx ON prot ( numPepOKMC2 );'
createIndex5 = 'CREATE INDEX pep_ok_idx on pep(OKforMS);'
createIndex6 = 'CREATE INDEX pep_MC_idx on pep(numMissedCleavages);'
createIndex7 = 'CREATE INDEX pep_acc_idx on pep(Accession);'
createIndex8 = 'CREATE INDEX prot_AccOKMC on pep (Accession, OKforMS, numMissedCleavages);'
createIndex9 = 'CREATE INDEX pep_nmnxs_idx on pep (numMotifsNXS);'
createIndex10 = 'CREATE INDEX pep_accMCnmnxs_idx on pep (Accession, numMissedCleavages, numMotifsNXS);'
createIndex11 = 'CREATE INDEX pep_accMCnmnxt_idx on pep (Accession, numMissedCleavages, numMotifsNXT);'
createIndex12 = 'CREATE INDEX pep_accMCnmnxc_idx on pep (Accession, numMissedCleavages, numMotifsNXC);'
createIndex13 = 'CREATE INDEX pep_accMCnmnxv_idx on pep (Accession, numMissedCleavages, numMotifsNXV);'
createIndex14 = 'CREATE INDEX pep_accMCnmc_idx on pep (Accession, numMissedCleavages, numMotifsC);'
createIndex15 = 'CREATE INDEX pep_accMCnmk_idx on pep (Accession, numMissedCleavages, numMotifsK);'

def create_cirfess_structure() :
    try:
        conn = sqlite3.connect('../cirfess.db')
    except:
        print( 'Connection Error: ' + str(Error) )
    cursor = conn.cursor()
    try:
        cursor.execute('DROP TABLE IF EXISTS pep;')
        cursor.execute('DROP TABLE IF EXISTS prot;')
    except:
        print( 'Error dropping tables: ' + str(error) )
    try:
        cursor.execute(createProtSQL)
        cursor.execute(createPepSQL)
    except:
        print( 'Table creation error: ' + str(Error) )
    try:
        cursor.execute(createIndex1)
        cursor.execute(createIndex2)
        cursor.execute(createIndex3)
        cursor.execute(createIndex4)
        cursor.execute(createIndex5)
        cursor.execute(createIndex6)
        cursor.execute(createIndex7)
        cursor.execute(createIndex8)
        cursor.execute(createIndex9)
        cursor.execute(createIndex10)
        cursor.execute(createIndex11)
        cursor.execute(createIndex12)
        cursor.execute(createIndex13)
        cursor.execute(createIndex14)
        cursor.execute(createIndex15)
    except:
        print( 'Index creation error: ' + str(Error) )
    conn.commit()

create_cirfess_structure()
