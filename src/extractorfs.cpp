// Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

/* Searches for all possible open reading frames (ORFs) in input FASTA format
   data in the context of given alphabet and translation table. Saves ORFs and
   their translations to FASTA files or prints them to stdout. Prints diagnostic
   messages to stderr.
   Required and [ optional ] command line arguments:
   --alph      : a path to a text file with sequence alphabet (first line - all
                 possible characters, second - complementary chracters)
   --tab       : a path to a text file with a translation table in a short format
                 as on https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
   --asmacc    : assembly accession version to be saved as metadata in sequence
                 headers
   --minlen    : minimal ORF's length in nucleotides
   [ --in ]    : optional, input FASTA file (uncompressed), if not given, data is
                 read from stdin
   [ --seqs ]  : optional, output FASTA file with ORF's sequeces, if not given,
                 data is written to stdout
   [ --trans ] : optional, output FASTA file with ORF's translations, if not given,
                  data is written to stdout
   USAGE:
   ./extractorfs --alph ALPHABET --tab TRANS_TAB --asmacc ASM_ACC --minlen MIN_LEN \
                 [ --in INPUT_FASTA --seqs ORF_FASTA --trans TRNAS_FASTA ]
*/
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string_view>
#include <iterator>
#include <algorithm>
#include <vector>
#include <set>
#include <map>

using namespace std;

/* USAGE message to be shown when required arguments are not provided.
*/
const char USAGE[] = "USAGE: extractorfs "
                     "--alph SEQ_ALPH_FILE "
                     "--tab TRANS_TAB_FILE "
                     "--asmacc ASSEMBLY_ACCESSION "
                     "--minlen ORF_MINLEN "
                     "[--in INPUT_FNA_FILE] "
                     "[--seqs OUTPUT_FNA_FILE] "
                     "[--trans OUTPUT_FAA_FILE]";

//------------------------------------------------------------------------------
/* Alphabet class used for reading an alphabet from a file and storing it in a
   map<char, char>, letter -> complementary letter.
*/
class Alphabet {
    public:
        map<char, char>* mapping;
        Alphabet();
        ~Alphabet();
        static Alphabet* fromFile(string);
};

/* Initialises a new object.
*/
Alphabet::Alphabet() {
    this->mapping = new map<char, char>();
}

/* Cleans up on destroy.
*/
Alphabet::~Alphabet() {
    delete this->mapping;
}

/* Reads an alphabet form a text file. Arguments:
   string fpath : a path to a text file with the alphabet
   Resturns:
   alph : Alphabet* object pointer
*/
Alphabet* Alphabet::fromFile(string fpath) {
    Alphabet* alph = new Alphabet();
    vector<string> lines;
    string line;
    
    /* Load the lines from a file to a vector of strings
    */
    ifstream inFile(fpath);
    while(getline(inFile, line)) {
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        lines.push_back(line);
    }
    inFile.close();
    
    /* The expected number of lines is 2 (alphabet characters and complementary
       characters).
    */
    if(lines.size() != 2) return NULL;
    
    /* Every charater from the second line is complementary to the character
       from the first line at the same position. Based on this assumption
       fill up the <char,char> map.
    */
    for(size_t i=0; i<lines[0].size(); i++) {
        char letter_one = lines[0][i];
        char letter_two = lines[1][i];
        (*alph->mapping)[letter_one] = letter_two;
    }
    
    return alph;
}

//------------------------------------------------------------------------------
/* CodonTable class used for reading a translation table from a file and storing
   it in a map<string, char>, three letter codon -> one letter amino acid residue.
   A CodonTable object also contains a set<string> of start and stop codons
   important for searching ORFs.
*/
class CodonTable {
    public:
        set<string>* starts;
        set<string>* stops;
        map<string, char>* mapping;
        CodonTable();
        ~CodonTable();
        static CodonTable* fromFile(string);
};

/* Initialises a new object.
*/
CodonTable::CodonTable() {
    this->starts = new set<string>();
    this->stops  = new set<string>();
    this->mapping = new map<string, char>();
}

/* Cleans up on destroy.
*/
CodonTable::~CodonTable() {
    delete this->starts;
    delete this->stops;
    delete this->mapping;
}

/* Loads codon table from a file. Arguments:
   string fpath -- a path to a text file with the translation table
   Returns:
   tab : CodonTable* object pointer
*/
CodonTable* CodonTable::fromFile(string fpath) {
    CodonTable* tab = new CodonTable();
    vector<string> lines;
    string line;
    
    /* Load the lines from a file to a vector of strings
    */
    ifstream inFile(fpath);
    while(getline(inFile, line)) {
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        lines.push_back(line);
    }
    inFile.close();
    
    /* The expected number of lines is 5 (first, second and third letter of a codon,
       whether the codon is a start or stop codon, encoded amio acid residue).
    */
    if(lines.size() != 5) return NULL;
    
    /* Characters at the same position of subsequent lines represents codon triplets,
       whether a codon is a start (M) or stop (*) codon and encoded amino acid residue.
       Based on these assumptions fill up a <string,char> map, and <string> vectors
       of start and stop codon.
    */
    for(size_t i=0; i<lines[0].size(); i++) {
        char residue = lines[0][i];
        char special = lines[1][i];
        string codon = string({lines[2][i], lines[3][i], lines[4][i]});
        if(special == 'M') (*tab->starts).insert(codon);
        else if(special == '*') (*tab->stops).insert(codon);
        (*tab->mapping)[codon] = residue;
    }
    
    return tab;
}

//------------------------------------------------------------------------------
/* Seq top-level class used for storing and manipulating biological sequecnes.
*/
class Seq {
    public:
        string seqid, title, seq, srcid;
        size_t start, end;
        Seq(string, string, string, string, size_t, size_t);
        string fasta(bool, size_t);
        static bool cmpSeqs(Seq*, Seq*);
        static vector<Seq*>* fromFile(string);
};

/* Initialises a Seq object. Arguments:
   string seqid : sequence id
   string title : sequence metadata in a one-line string
   string seq   : sequence, without any extra characters, e.g. new-line charaters
   string srcid : sequence id of parental sequence
   size_t start : start location within the parental sequence (GenBank notation)
   size_t end   : end location within the parental sequence (GenBank notation),
                  if end < start the sequence is located on minus strand
*/
Seq::Seq(string seqid, string title="", string seq="", string srcid="", size_t start=0, size_t end=0) {
    this->seqid = seqid;
    this->title = title;
    this->seq   = seq;
    this->srcid = srcid;
    this->start = start;
    this->end   = end;
}

/* Returns a stored sequecne data in FASTA format, by default save metadata in
   in header and set line width to 60 letters. Arguments:
   bool meta : whether to write sequence meta data to FASTA header, default true
   size_t lw : line length in characters
   Returns:
   string fasta : a string with the sequence written in FASTA format
*/
string Seq::fasta(bool meta=true, size_t lw=60) {
    stringstream fasta;
    fasta << ">" << this->seqid;
    if(this->title != "") fasta << " " << this->title;
    if(meta) {
        if(this->srcid != "") fasta << " srcid=" << this->srcid;
        if(this->start != 0)  fasta << " start=" << this->start;
        if(this->end   != 0)  fasta << " end="   << this->end;
    }
    fasta << endl;
    for(size_t i=0; i<this->seq.length(); i+=lw) {
        fasta << this->seq.substr(i, lw) << endl;
    }
    return fasta.str();
}

/* Compares two Seq objects for sorting purposes in respect to their location
   (start, end) in a source sequecne. Arguments:
   Seq* one   : a sequence to compare
   Seq* other : annother sequence to compare
   Returns:
   bool lower : true if the first sequence preceeds the other, otherwise false
*/
bool Seq::cmpSeqs(Seq* one, Seq* other){
    size_t coords[2];
    Seq* seqs[] = {one, other};
    for(size_t i=0; i<2; i++) {
        Seq* seq = seqs[i];
        if(seq->start <= seq->end) coords[i] = seq->start;
        else coords[i] = seq->end;
    }
    bool lower = coords[0] <= coords[1];
    return lower;
}

/* Loads sequcens from a FASTA file and return them in a vector. If the file path
   is not provided, read sequences form stdin. Arguments:
   string fpath : a path to a FASTA file with sequences
   Returns:
   vector<Seq*>* seqs : vector<Seq*>* object pointer
*/
vector<Seq*>* Seq::fromFile(string fpath_in) {
    vector<Seq*>* seqs = new vector<Seq*>();
    Seq* cur_seq = NULL;
    string line;
    
    istream* stream = &cin;
    ifstream inFile;
    if(fpath_in != "") {
        inFile.open(fpath_in);
        stream = &inFile;
    }
    while(getline(*stream, line)) {
        line.erase(0, line.find_first_not_of(" \n\r\t"));
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        if(line == "") continue;
        if(line[0] == '>') {
            string seqid;
            string title;
            size_t pos = line.find(' ', 1);
            if(pos == string::npos) seqid = line.substr(1);
            else {
                seqid = line.substr(1, pos-1);
                title = line.substr(pos+1);
            }
            cur_seq = new Seq(seqid, title);
            seqs->push_back(cur_seq);
        } else {
            cur_seq->seq += line;
        }
    }
    if(fpath_in != "") inFile.close();
    
    return seqs;
}

//------------------------------------------------------------------------------
/* ProtSeq class derived form Seq used for storing and manipulating protein
   sequences.
*/
class ProtSeq: public Seq {
    public:
        using Seq::Seq;
};
//------------------------------------------------------------------------------
/* DNASeq class derived form Seq used for storing and manipulating nucleotide
   sequences.
*/
class DNASeq: public Seq {
    private:
        static char cmpl(Alphabet*, char);
    public:
        using Seq::Seq;
        DNASeq* revCmpl(Alphabet*);
        ProtSeq* trans(CodonTable*, short, bool);
        vector<DNASeq*>* findORFs(size_t, Alphabet*, CodonTable*);
};

/* Returns a complementary letter to a given one. Arguments:
   Alphabet* alph : nucleotide sequence alphabet to be used
   char letter    : letter a complementary one is to be returned
   Returns:
   char cmpl_letter : a complemenatary letter to the given one
*/
char DNASeq::cmpl(Alphabet* alph, char letter) {
    char cmpl_letter;
    if(alph->mapping->find(letter) != alph->mapping->end()) {
        cmpl_letter = (*alph->mapping)[letter];
    } else letter = 'X';
    return cmpl_letter;
}
/* Returns a reverse and complementary sequence. Arguments:
   Alphabet* alph : nucleotide sequence alphabet to be used
   Returns:
   DNASeq* seq : DNASeq* object pointer
*/
DNASeq* DNASeq::revCmpl(Alphabet* alph) {
    string raw_seq = this->seq;
    reverse(raw_seq.begin(), raw_seq.end());
    transform(raw_seq.begin(), raw_seq.end(), raw_seq.begin(),
              [alph](char letter){return DNASeq::cmpl(alph, letter);});
    DNASeq* seq = new DNASeq(this->seqid, this->title, raw_seq, this->srcid,
                             this->start, this->end);
    return seq;
}

/* Translates a DNA sequence into amino acid one. Arguments:
   CodonTable* tab : codon table to be used for translation
   short start     : open reading frame shift (e.g. 0, 1 or 3), default 0
   bool tostop     : stop translation when stop codon is encountered, dafault true
   Returns:
   ProtSeq* seq : ProtSeq* object pointer
*/
ProtSeq* DNASeq::trans(CodonTable* tab, short start=0, bool tostop=true) {
    stringstream trans;
    //string_view raw_seq = this->seq;
    //size_t len = (raw_seq.length() - start) / 3 * 3;
    size_t len = (this->seq.length() - start) / 3 * 3;
    
    for(size_t i=start; i<len; i+=3) {
        string codon = (string)this->seq.substr(i, 3);
        char residue = 'X';
        if(tab->mapping->find(codon) != tab->mapping->end())
            residue = (*tab->mapping)[codon];
        if(tostop && residue == '*') break;
        trans << residue;
    }
    
    ProtSeq* seq = new ProtSeq(this->seqid, this->title, trans.str(), this->srcid,
                               this->start, this->end);
    return seq;
}

/* Searches for open reading frames (ORFs) in the DNA sequence. Arguments:
   size_t minlen   : minimal length of an ORF
   Alphabet* alph  : nucleotide sequence alphabet to be used
   CodonTable* tab : codon table to be used for translation (includes information
                     on start and stop codons)
   Returns:
   vector<DNASeq*>* orfs : vector<DNASeq*>* object pointer
*/
vector<DNASeq*>* DNASeq::findORFs(size_t minlen, Alphabet* alph, CodonTable* tab) {
    vector<DNASeq*>* orfs = new vector<DNASeq*>();
    DNASeq* revcmpl = this->revCmpl(alph);
    size_t count = 0;
    /* Iterate over the sequence and its reverse complement.
    */
    for(DNASeq* seqobj : {this, revcmpl}) {
        short strand = seqobj == this ? 1 : -1;
        string_view seq = seqobj->seq;
        /* Iterate over all 3 possible shifts of open reading frames.
        */
        for(unsigned short frame=0; frame<3; frame++)
        {
            size_t len = (seqobj->seq.length() - frame) / 3 * 3;
            size_t start = string::npos;
            /* Iterate over subsequent codon triplets.
            */
            for(size_t i=frame; i<len; i+=3) {
                string codon = (string)seq.substr(i, 3);
                /* If we already have encountered a start codon.
                */
                if(start != string::npos) {
                    /* If codon is a stop codon, check for an ORF.
                    */
                    if(tab->stops->find(codon) != tab->stops->end()) {
                        size_t end = i+3;
                        size_t orflen = end - start;
                        /* Check whether the ORF meets the length requirement.
                           If it does, generate a new sequence id from
                           the sequecne id and a subsequent number, convert start
                           and stop coordinates to GenBank format, and finally
                           create a new DNASeq object and push a pointer to it
                           to the final vector.
                        */
                        if(orflen >= minlen) {
                            string orfseq = (string)seq.substr(start, orflen);
                            char seqid[100];
                            sprintf(seqid, "%s_%06d", seqobj->seqid.c_str(), count);
                            start ++;
                            if(strand == -1) {
                                start = seqobj->seq.length() - start + 1;
                                end = seqobj->seq.length() - end + 1;
                            }
                            DNASeq* orf = new DNASeq(seqid, "", orfseq,
                                                     seqobj->seqid, start, end);
                            orfs->push_back(orf);
                            count ++;
                        }
                        start = string::npos;
                    }
                /* If we have not and the current codon is a start codon, save
                   its position.
                */
                } else if(tab->starts->find(codon) != tab->starts->end()) {
                    start = i;
                }
            }
        }
    }
    return orfs;
}
//------------------------------------------------------------------------------
/* A helper class for stroing command line arguments.
*/
class Args {
    public:
        // required
        string alph;
        string tab;
        string asmacc;
        string minlen;
        // optional
        string in;
        string seqs;
        string trans;
        static Args* parseArgs(int, char**);
};

/* Parses command line arguments. Arguments:
   int argc    : main entry function argument
   char** argv : main entry function argument
   Returns:
   Args* args : Args* object pointer
*/
Args* Args::parseArgs(int argc, char** argv) {
    Args* args = new Args();
    string* cur_arg = NULL;
    
    for(int i=1; i<argc; i++) {
        string_view arg = argv[i];
        if(arg.substr(0, 2) == "--") {
            if     (arg == "--alph")   cur_arg = &args->alph;
            else if(arg == "--tab")    cur_arg = &args->tab;
            else if(arg == "--asmacc") cur_arg = &args->asmacc;
            else if(arg == "--minlen") cur_arg = &args->minlen;
            else if(arg == "--in")     cur_arg = &args->in;
            else if(arg == "--seqs")   cur_arg = &args->seqs;
            else if(arg == "--trans")  cur_arg = &args->trans;
        } else if(cur_arg != NULL) {
            *cur_arg = arg;
            cur_arg = NULL;
        }
    }
    
    for(string* ptr=&args->alph; ptr <= &args->minlen; ptr++) {
        if(*ptr == "") {
            args = NULL;
            break;
        }
    }
    
    return args;
}
//------------------------------------------------------------------------------
/* Entry function that performs searching for all possible open reading frames
   (ORFs) in input FASTA format data in the context of given alphabet and
   translation table. Saves ORFs and their translations to FASTA files or prints
   them to stdout. Prints diagnostic messages to stderr.
*/
int main(int argc, char *argv[]) {
    /* Parse arguments.
    */
    Args *args = Args::parseArgs(argc, argv);
    if(args == NULL) {
        cerr << USAGE << endl;
        return 1;
    }
    
    /* Load the alphabet, the codon table from files, and finally the sequences
       from indicated files.
    */
    Alphabet* alph = Alphabet::fromFile(args->alph);
    CodonTable*  tab  = CodonTable::fromFile(args->tab);
    vector<DNASeq*>* seqs = (vector<DNASeq*>*)Seq::fromFile(args->in);
    
    /* Create and empty vector for ORFs (DNASeq* object pointers). Iterate over
       sequences and search for ORFs in each.
    */
    vector<DNASeq*> orfs;
    vector<DNASeq*>::iterator ptr;
    int minlen;
    sscanf(args->minlen.data(), "%d", &minlen);
    for(ptr=seqs->begin(); ptr!=seqs->end(); ptr++) {
        DNASeq* seq = *ptr;
        cerr << "Processing sequence: " << seq->seqid << endl;
        vector<DNASeq*>* orfs_batch = seq->findORFs(minlen, alph, tab);
        sort(orfs_batch->begin(), orfs_batch->end(), Seq::cmpSeqs);
        cerr << "In sequence " << seq->seqid << " of length "
            << seq->seq.length() << " found " << orfs_batch->size() << " orfs"
            << endl;
        orfs.insert(orfs.end(), orfs_batch->begin(), orfs_batch->end());
        delete orfs_batch;
    }
    cerr << "Find total " << orfs.size() << " orfs" << endl;
    
    /* Print in FASTA format ORF sequences and their translations to stdout
       or save to files.
    */
    ostream* seqstream = &cout;
    ofstream outFileSeqs;
    if(args->seqs != "") {
        outFileSeqs.open(args->seqs);
        seqstream = &outFileSeqs;
    }
    ostream* transtream = &cout;
    ofstream outFileTrans;
    if(args->trans != "") {
        outFileTrans.open(args->trans);
        transtream = &outFileTrans;
    }
    for(ptr=orfs.begin(); ptr!=orfs.end(); ptr++) {
        DNASeq* seq = *ptr;
        seq->title = "asmacc=" + args->asmacc + seq->title;
        *seqstream  << seq->fasta()             << endl;
        *transtream << seq->trans(tab)->fasta() << endl;
    }
    if(args->seqs != "") outFileSeqs.close();
    if(args->trans != "") outFileTrans.close();
    
    /* Clean up.
    */
    delete args, alph, tab, seqs;
    cerr << "Done" << endl;
}

