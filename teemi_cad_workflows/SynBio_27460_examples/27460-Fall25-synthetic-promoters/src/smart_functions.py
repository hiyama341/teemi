#rasmus' yndlingsfunktioner

from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord

def read_fasta_to_dseqrecords(fasta_path):
    """Return (list_of_Dseqrecord, list_of_names)"""
    recs = []
    names = []
    for r in SeqIO.parse(fasta_path, "fasta"):
        seq = str(r.seq).upper().replace("\n", "")
        dr = Dseqrecord(seq)
        try:
            dr.name = r.id
        except Exception:
            pass
        recs.append(dr)
        names.append(r.id)
    return recs, names