#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <map>

using namespace std;

string read_reference(const string& filename) {
    fstream file(filename, fstream::in);

    string ref;
    string tmp;
    while (!file.eof()) {
        getline(file, tmp);
        if (tmp[0] == '>') {
            // cout << tmp << endl;
        } else {
            ref += tmp;
        }
    }
    file.close();
    return ref;
}

int*** init_3d_array(uint x, uint y, uint z) {
    // cout << "Initializing array of size (" << x << ", " << y << ", " << z << ")." << endl;
    int*** res = new int**[x];
    for (int i=0; i<x; i++) {
        res[i] = new int*[y];
        for (int j=0; j<y; j++) {
            res[i][j] = new int[z];
            for (int k=0; k<z; k++) {
                res[i][j][k] = 0;
            }
        }
    }
    return res;
}

void print_3d_array(int*** A, uint x, uint y, uint z) {
    for (int i=0; i<x; i++) {
        for (int j=0; j<y; j++) {
            for (int k=0; k<z; k++) {
                cout << A[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl << endl;
    }
}

void delete_3d_array(int ***A, uint x, uint y, uint z) {
    // cout << "Deleting array of size (" << x << ", " << y << ", " << z << ")." << endl;
    for (int i=0; i<x; i++) {
        for (int j=0; j<y; j++) {
            delete[] A[i][j];
        }
        delete[] A[i];
    }
    delete[] A;
}

void read_mutations(
    const string& filename, vector<string>& variants, map<string, string>& snps, map<string, uint>& thresholds
) {
    fstream file(filename, fstream::in);
    string tmp;
    string variant;
    uint thr;
    while (!file.eof()) {
        getline(file, tmp);
        stringstream ss(tmp);
        ss >> variant >> thr;
        variants.push_back(variant);
        thresholds.insert(pair(variant, thr));
        string snp;
        while (!ss.eof()) {
            ss >> snp;
            snps.insert(pair(snp, variant));
        }
    }
}

class sam_record {
public:
    string qname, rname, cigar, rnext, seq, qual;
    int flag, pos, mapq, pnext, tlen;

    sam_record() {
        flag = pos = mapq = pnext = tlen = 0;
        qname = rname = cigar = rnext = seq = qual = "";
    }

    explicit sam_record(const string& line) : sam_record(){
        stringstream ss(line);
        ss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;
    }
};

class record {
public:
    string qname;   // query name, 1st column in the SAM file
    int pos;        // position of the alignment, 4th column in the SAM file
    string seq;     // aligned sequence, 10th column, in the SAM file

    record() {qname=""; pos=0; seq="";}
    explicit record(uint length) : record() {
        seq = string(length, 'N');
    }
};

string read_next_seq(fstream& aln, uint n){
    static sam_record previous;
    record record(n);

    string line;
    while (!aln.eof()) {
        getline(aln, line);
        if (line[0] == '@') { continue; }
        if (line.length() == 0) { continue; }
        sam_record current = sam_record(line);
        // cout << current.cigar << " ";
        for (int i=0; i<current.seq.length(); i++) {
            record.seq[current.pos - 1 + i] = current.seq[i];
        }

        if (previous.qname != current.qname) {
            // cout << endl;
            previous = current;
            return record.seq;
        }
    }
    return "";
}

tuple<char, uint, char> parse_snp(string s) {
    uint n = s.length();
    tuple<char, uint, char> result = tuple(s[0], stoi(s.substr(1, n-1)), s[n-1]);
    return result;
}

string get_seq_variant(const string& seq, const map<string, string>& snp) {
    for (const auto& s : snp) {
        tuple<char, uint, char> mut = parse_snp(s.first);

        if (seq[get<1>(mut)-1] == get<2>(mut)) {
            cout << get<0>(mut) << " " << get<1>(mut) << " " << get<2>(mut) << " " << s.second << endl;
        }
    }
    return "";
}

int main() {
    // inputs:
    string ref_filename = "/home/andy/projects/cpp_freqtable/data/covid_ref.fa";
    string aln_filename = "/home/andy/projects/cpp_freqtable/data/covid_aln.sam";
    string mut_filename = "/home/andy/projects/cpp_freqtable/data/mut.txt";
    // outputs:
    string tsv_filename = "/home/andy/projects/cpp_freqtable/data/full_out.tsv";

    // read reference
    string ref = read_reference(ref_filename);
    cout << ref << endl;

    // read mutations
    vector<string> variants;
    map<string, string> snp;
    map<string, uint> threshold;
    read_mutations(mut_filename, variants, snp, threshold);

    // initialize alphabet
    vector<char> alphabet = {'A', 'C', 'G', 'T', 'N', '-'};

    int*** A = init_3d_array(ref.length(), variants.size(), alphabet.size());


    // for each virus
    //      read sequence
    //      assign variant
    //      add counts to fullout

    fstream aln(aln_filename, fstream::in);
    while (!aln.eof()) {
        string seq = read_next_seq(aln, ref.length());
        string variant = get_seq_variant(seq, snp);
        cout << seq << variant << endl;

    }

    // store fullout
    delete_3d_array(A, ref.length(), variants.size(), alphabet.size());
    return 0;
}
