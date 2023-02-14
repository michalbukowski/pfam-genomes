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

const char USAGE[] = "USAGE: extractorfs "
                     "--alph SEQ_ALPH_FILE "
                     "--tab TRANS_TAB_FILE "
                     "--asmacc ASSEMBLY_ACCESSION "
                     "--minlen ORF_MINLEN "
                     "[--in INPUT_FNA_FILE] "
                     "[--seqs OUTPUT_FNA_FILE] "
                     "[--trans OUTPUT_FAA_FILE] "
                     "[--quiet]";

//------------------------------------------------------------------------------
class Alphabet {
    public:
        map<char, char>* mapping;
        Alphabet();
        ~Alphabet();
        static Alphabet* fromFile(string);
};

Alphabet::Alphabet() {
    this->mapping = new map<char, char>();
}

Alphabet::~Alphabet() {
    delete this->mapping;
}

Alphabet* Alphabet::fromFile(string fpath) {
    Alphabet* alph = new Alphabet();
    vector<string> lines;
    string line;
    
    ifstream inFile(fpath);
    while(getline(inFile, line)) {
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        lines.push_back(line);
    }
    inFile.close();
    
    if(lines.size() != 2) return NULL;
    
    for(size_t i=0; i<lines[0].size(); i++) {
        char letter_one = lines[0][i];
        char letter_two = lines[1][i];
        (*alph->mapping)[letter_one] = letter_two;
    }
    
    return alph;
}

//------------------------------------------------------------------------------
class CodonTable {
    public:
        set<string>* starts;
        set<string>* stops;
        map<string, char>* mapping;
        CodonTable();
        ~CodonTable();
        static CodonTable* fromFile(string);
};

CodonTable::CodonTable() {
    this->starts = new set<string>();
    this->stops  = new set<string>();
    this->mapping = new map<string, char>();
}

CodonTable::~CodonTable() {
    delete this->starts;
    delete this->stops;
    delete this->mapping;
}

CodonTable* CodonTable::fromFile(string fpath) {
    CodonTable* tab = new CodonTable();
    vector<string> lines;
    string line;
    
    ifstream inFile(fpath);
    while(getline(inFile, line)) {
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        lines.push_back(line);
    }
    inFile.close();
    
    if(lines.size() != 5) return NULL;
    
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
class Seq {
    public:
        string seqid, title, seq, srcid;
        size_t start, end;
        Seq(string, string, string, string, size_t, size_t);
        string fasta(bool, size_t);
        static bool cmpSeqs(Seq*, Seq*);
        static vector<Seq*>* fromFile(string);
};

Seq::Seq(string seqid, string title="", string seq="", string srcid="", size_t start=0, size_t end=0) {
    this->seqid = seqid;
    this->title = title;
    this->seq   = seq;
    this->srcid = srcid;
    this->start = start;
    this->end   = end;
}

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

bool Seq::cmpSeqs(Seq* one, Seq* other){
    size_t coords[2];
    Seq* seqs[] = {one, other};
    for(size_t i=0; i<2; i++) {
        Seq* seq = seqs[i];
        if(seq->start <= seq->end) coords[i] = seq->start;
        else coords[i] = seq->end;
    }
    return coords[0] <= coords[1];
}

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
class ProtSeq: public Seq {
    private:
    public:
        using Seq::Seq;
};
//------------------------------------------------------------------------------
class DNASeq: public Seq {
    private:
        static char cmpl(Alphabet*, char);
    public:
        using Seq::Seq;
        DNASeq* revCmpl(Alphabet*);
        ProtSeq* trans(CodonTable*, short, bool);
        vector<DNASeq*>* findORFs(size_t, Alphabet*, CodonTable*);
};

char DNASeq::cmpl(Alphabet* alph, char letter) {
    if(alph->mapping->find(letter) != alph->mapping->end()) {
        return (*alph->mapping)[letter];
    }
    return 'X';
}

DNASeq* DNASeq::revCmpl(Alphabet* alph) {
    string seq = this->seq;
    reverse(seq.begin(), seq.end());
    transform(seq.begin(), seq.end(), seq.begin(),
              [alph](char letter){return DNASeq::cmpl(alph, letter);});
    return new DNASeq(this->seqid, this->title, seq, this->srcid,
                      this->start, this->end);
}

ProtSeq* DNASeq::trans(CodonTable* tab, short start=0, bool tostop=true) {
    stringstream trans;
    string_view seq = this->seq;
    size_t len = (seq.length() - start) / 3 * 3;
    
    for(size_t i=start; i<len; i+=3) {
        string codon = (string)this->seq.substr(i, 3);
        char residue = 'X';
        if(tab->mapping->find(codon) != tab->mapping->end())
            residue = (*tab->mapping)[codon];
        if(tostop && residue == '*') break;
        trans << residue;
    }
    
    return new ProtSeq(this->seqid, this->title, trans.str(), this->srcid,
                       this->start, this->end);
}

vector<DNASeq*>* DNASeq::findORFs(size_t minlen, Alphabet* alph, CodonTable* tab) {
    vector<DNASeq*>* orfs = new vector<DNASeq*>();
    DNASeq* revcmpl = this->revCmpl(alph);
    size_t count = 0;
    for(DNASeq* seqobj : {this, revcmpl}) {
        short strand = seqobj == this ? 1 : -1;
        string_view seq = seqobj->seq;
        for(unsigned short frame=0; frame<3; frame++)
        {
            size_t len = (seqobj->seq.length() - frame) / 3 * 3;
            size_t start = string::npos;
            for(size_t i=frame; i<len; i+=3) {
                string codon = (string)seq.substr(i, 3);
                if(start != string::npos) {
                    if(tab->stops->find(codon) != tab->stops->end()) {
                        size_t end = i+3;
                        size_t orflen = end - start;
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
                } else if(tab->starts->find(codon) != tab->starts->end()) {
                    start = i;
                }
            }
        }
    }
    return orfs;
}
//------------------------------------------------------------------------------
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
int main(int argc, char *argv[]) {
    Args *args = Args::parseArgs(argc, argv);
    if(args == NULL) {
        cerr << USAGE << endl;
        return 1;
    }
    
    Alphabet* alph = Alphabet::fromFile(args->alph);
    CodonTable*  tab  = CodonTable::fromFile(args->tab);
    vector<DNASeq*>* seqs = (vector<DNASeq*>*)Seq::fromFile(args->in);
    
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
    
    delete args, alph, tab, seqs;
    cerr << "Done" << endl;
}

