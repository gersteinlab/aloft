#maketrans is not a part of string module in Python3
#handle python 2.x and 3.x differently
try:
    #2.x
    import string
    comptable = string.maketrans("ATGC", "TACG")
except AttributeError:
    #3.x
    comptable = str.maketrans("ATGC", "TACG")

code={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',\
      'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',\
      'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',\
      'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',\
      'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',\
      'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',\
      'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',\
      'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',\
      'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',\
      'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',\
      'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',\
      'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',\
      'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',\
      'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',\
      'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',\
      'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

#reverses and complements string, eg from ACT to AGT
def compstr(strand):
    return strand.translate(comptable)[::-1]

#translates string, eg from ACT to T
def translate_aa(seq):
    return "".join(code[seq[i:i+3]] for i in range(0, len(seq)-2, 3))