# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

import gzip
from os import linesep as eol

class Seq:
    '''A helper class for DNA sequence handling. Each Seq object is described by
       sequence id (seqid), one-line string with metadata (title) and
       the sequence itself (seq, wihout extra characters, e.g. newline characters).
       If in the header (title) metadata are present in a 'key=value' format,
       these are assigned to a Seq obejct as its properties.
    '''
    def __init__(self, seqid, title, seq):
        self.seqid = seqid
        self.title = title
        self.seq   = seq
        # Check for 'key=value' pairs in the header line, if present, save to
        # dictionary and update with it the self.__dict__ object to make
        # turn them into properties of a new Seq object.
        items = [ field.split('=', 1) for field in title.split(' ')
                  if '=' in field ]
        if len(items) > 0:
            meta = dict(items)
            self.__dict__.update(meta)
    
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
        '''Returns the sequence stored in the Seq object as a string compliant
           with FASTA format. Arguments:
           lw : line length in characters, default 60
        '''
        return fasta(self.seqid, self.title, self.seq, lw)

def fasta(seqid, title, seq, lw=60):
    '''Given a sequence id, metadata and a sequence, generates its representation
       in FASTA format. Arguments:
       seqid : squence id
       title : one-line string with metadata
       seq   : sequence wihout extra characters, e.g. newline characters
       lw    : line length in characters, default 60
       Returns:
       fasta : string with the sequence written in FASTA format
    '''
    lines = eol.join(seq[i:i+lw] for i in range(0, len(seq), lw))
    title = ' ' + title if title != '' else ''
    fasta = f'>{seqid}{title}{eol}{lines}{eol}'
    return fasta

def fasta_meta(header, meta_keys):
    '''Tries to read values for metadata keys from a given FASTA header. Raises
       a KeyError if any of the keys is not found. Arguments:
       header    : string with the sequence header
       meta_keys : a list of keys expected to be present in the header as 'key=value' pairs
       Returns:
       meta : a dictionary of {key : value} pairs
    '''
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

def read_fasta(fpath, seqids=None):
    '''Reads sequecnes from a FASTA file (or GZIPped FASTA file) to Seq objects.
       If particular sequence ids are requesed and any of them is not found
       in the file, a KeyError is raised. Arguments:
       fpath  : FASTA file path
       seqids : a collection of sequence ids that are suposed to be read from
                the file, dafault None (read all sequences)
       Returns:
       seqs : a dictionary of {seqid : Seq} pairs
    '''
    if seqids is not None:
        seqids = set(seqids)
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
        raise KeyError(f'Cannot find {seqids} in {fpath}')
    
    return seqs

