import sqlite3
import argparse

def gen_seqs(db, out, seq_type):
    # TODO: set up so seq_type actually does something
    conn = sqlite3.connect(db)
    c = conn.cursor()
    with open(out, 'w') as f:
        for i in range(10500):
            t = (str(i),)
            exec_str = 'SELECT sequences.%s FROM sequences, homology WHERE homology.component_id=? AND sequences.catalog_id="ILLUMINA_0171-6-AMPHIPLICA" AND homology.run_id=416 AND homology.sequence_id=sequences.sequence_id LIMIT 1;' % (seq_type)
            c.execute(exec_str, t)
            seq = c.fetchone()
            if seq != None:
                f.write('%s|%s\n' % (i, seq[0]))

if __name__ == '__main__':

    desc = "generates root sequences from an agalma homologize run"
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-db', '--database', help='agalma database')
    parser.add_argument('-o', '--out', help='out file')
    parser.add_argument('-t', '--type', default='nucleotide', help='nucleotide or protein')

    opts = parser.parse_args()
    db = opts.database
    out = opts.out
    seq_type = opts.type

    if seq_type == "nucleotide":
        seq_type = "nucleotide_seq"
    elif seq_type == "protein":
        seq_type = "protein_seq"
    else:
        print "seq_type is an error"

    gen_seqs(db,out,seq_type)
