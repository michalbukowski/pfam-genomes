import gzip, io
import pandas as pd
from os import linesep as eol
from collections import defaultdict

dna_alphabet = {
    'A' : 'T', 'G' : 'C', 'C' : 'G', 'T' : 'A',
    'Y' : 'R', 'R' : 'Y', 'W' : 'W', 'S' : 'S',
    'K' : 'M', 'M' : 'K', 'D' : 'H', 'V' : 'B',
    'H' : 'D', 'B' : 'V', 'N' : 'N', 'X' : 'X',
    '-': '-'
}

codon_tab = defaultdict(lambda: 'X', {
    'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
    'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
    'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
    'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

    'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R',
    'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R',
    'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R',
    'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R',

    'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S',
    'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S',
    'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R',
    'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R',

    'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G',
    'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G',
    'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G',
    'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'
})

starts = 'ATG TTG CTG GTG'.split()
stops = 'TAA TAG TGA'.split()

class Seq:
    def __init__(self, seqid, title, seq):
        self.seqid = seqid
        self.title = title
        self.seq   = seq
        items = [ field.split('=', 1) for field in title.split(' ') ]
        meta = dict(items)
        self.__dict__.update(meta)
    
    def __iter__(self):
        return self.seq
    
    def __getitem__(self, key):
        return self.seq[key]
    
    def __len__(self):
        return len(self.seq)
    
    def __repr__(self):
        if len(self.seq) <= 18:
            return self.seq
        else:
            return self.seq[:9] + '...' + self.seq[-9:]
    
    def __str__(self):
        return self.seq
    
    def fasta(self, lw=60):
        lines = eol.join(self[i:i+lw] for i in range(0, len(self), lw))
        title = ' ' + self.title if self.title != '' else ''
        fasta = f'>{self.seqid}{title}{eol}{lines}{eol}'
        return fasta

class Orf:
    def __init__(self, start, end, seq, trans):
        self.start = start
        self.end       = end
        self.seq       = seq
        self.trans     = trans

def fasta_meta(header, meta_keys):
    header = 'seqid=' + header
    meta = {}
    for key in meta_keys:
        start = header.find(f'{key}=')
        if start == -1:
            raise KeyError(f'Key "{key}" not present in metadata: {header}')
        end = header.find(' ', start)
        end = end if end != -1 else None
        value = header[start+len(key)+1:end]
        meta[key] = value
    return meta

def fasta(seqid, title, seq, lw=60):
    lines = eol.join(seq[i:i+lw] for i in range(0, len(seq), lw))
    title = ' ' + title if title != '' else ''
    fasta = f'>{seqid}{title}{eol}{lines}{eol}'
    return fasta

def revcmpl(seq):
    seq = list(seq)
    seq.reverse()
    seq_revcmpl = ''.join( [ dna_alphabet[letter] for letter in seq ] )
    return seq_revcmpl

def translate(seq):
    seq_len = len(seq)
    assert seq_len % 3 == 0, f'Cannot translate a sequence of {seq_len} length'
    codons = [ seq[i:i+3] for i in range(0, seq_len, 3) ]
    trans  = ''.join( [ codon_tab[codon] for codon in codons ] )
    return trans

def gb_location(seq, start, end=None, downstream=None, upstream=None):
    if end is not None:
        if type(start) is not int or type(end) is not int:
            raise TypeError('Start and end must be integers')
        if start <= 0 or end <= 0:
            raise KeyError(
                f'Minimal value for start and end is 1 ({start}..{end})'
            )
        
        strand = 1 if start <= end else -1
        start, end = sorted([start, end])
        start -= 1
        
        if strand == -1:
            upstream, downstream = downstream, upstream
        
        if upstream is not None:
            if start < upstream:
                upstream = start
            start -= upstream
        if downstream is not None:
            if end + downstream > len(seq):
                downstream = len(seq) - end
            end += downstream
        
        sub_seq = seq[start:end]
        if strand == -1:
            sub_seq = revcmpl(sub_seq)
            upstream, downstream = downstream, upstream
        
        if downstream is not None or upstream is not None:
            return sub_seq, start, end
        else:
            return sub_seq
    
    else:
        if type(end) is not int:
            raise TypeError('Index must be an integer')
        if end <= 0:
            raise KeyError('Minimal value of index is 1')
        return seq[key - 1]

def fetch_seqids(fpath):
    seqids = []
    f = gzip.open(fpath, 'rt', encoding='utf-8') if fpath.endswith('.gz') else open(fpath)
    for line in f:
        if line.startswith('>'):
            seqid = line[1:-1].split(' ', 1)[0]
            seqids.append(seqid)
    f.close()
    return seqids

def read_fasta(fpath, seqids=None):
    seqid = None
    seqs = {}

    f = gzip.open(fpath, 'rt', encoding='utf-8') if fpath.endswith('.gz') else open(fpath)
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            seqid, title = line[1:].split(' ', 1)
            if seqids is None or seqid in seqids:
                seqs[seqid] = Seq(seqid, title, '')
            else:
                seqid = None
        elif seqid is not None:
            seqs[seqid].seq += line
    f.close()
    
    if seqids is not None and len(seqs) != len(seqids):
        seqids = ', '.join(seqids - set(seqs.keys()))
        raise Exception(f'Cannot find {seqids} in {fpath}')
    
    return seqs

def read_gff3(fpath):
    cols = 'seqid source ftype start end score strand phase attrs'.replace(' ', '\t') + eol
    f = open(fpath)
    seqid = None
    assert f.readline().startswith('##gff-version')
    regions = {}
    for line in f:
        if line.startswith('##'):
            assert line.startswith('##sequence-region')
            seqid = line.rstrip().split(' ')[1]
            regions[seqid] = cols
        else:
            regions[seqid] += line
    regions = { seqid: pd.read_csv(io.StringIO(regions[seqid]), sep='\t')
                for seqid in regions }
    return regions

def find_orfs(seq, minlen):
    orfs = []
    for strand, seq in (1, seq), (-1, revcmpl(seq)):
        for shift in range(3):
            codons = [ seq[i:i+3] for i in range(shift, len(seq)-2, 3) ]
            start = None
            for i, codon in enumerate(codons):
                if codon in stops and start is not None:
                    if (i - start + 1)*3 >= minlen:
                        orfseq   = ''.join(codons[start:i+1])
                        trans = translate(orfseq)
                        start, end = start*3+shift + 1, (i+1)*3+shift
                        if strand == -1:
                            start = len(seq) - start + 1
                            end = len(seq) - end + 1
                        orfs.append(Orf(start, end, orfseq, trans))
                    start = None
                elif codon in starts and start is None:
                    start = i
    return orfs

