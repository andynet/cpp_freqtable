#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <map>
#include <cstdlib>

#define OTHER "other"

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

void store_3d_array(int ***A, uint length, const vector<string>& variants, const vector<char>& alphabet, ostream& out) {
    for (int i=0; i<variants.size(); i++) {
        for (int j=0; j<length; j++) {
            for (int k=0; k<alphabet.size(); k++) {
                out << variants.at(i) << "\t" << j << "\t" << alphabet.at(k) << "\t" << A[j][i][k] << endl;
            }
        }
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

void read_mutations(const string& filename, vector<string>& variants, map<string, string>& snps, map<string, uint>& thresholds) {
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
    // qname flag rname pos mapq cigar rnext pnext tlen seq qual tags
    string qname, cigar, seq;
    int pos = 0;

    explicit sam_record(const string& line){
        stringstream ss(line);
        string del;
        ss >> qname >> del >> del >> pos >> del >> cigar >> del >> del >> del >> seq >> del;
    }
};

class record {
public:
    string qname;   // query name, 1st column in the SAM file
    string seq;     // aligned sequence, 10th column, in the SAM file
    string variant; // predicted using mut file

    record() {qname=""; seq=""; variant="";}
    explicit record(uint length) : record() {
        seq = string(length, 'N');
    }
};

int ctoi(char c) {
    return (int)c - 48;
}

void parse_cigar(const string& cigar, vector<uint>& op_size, vector<char>& op_type) {
    uint size = 0;
    for (const char c : cigar) {
        if (isdigit(c)) {
            size = (size * 10) + ctoi(c);
        } else {
            op_size.push_back(size);
            op_type.push_back(c);
            size = 0;
        }
    }
}

void match(const string& source, uint& rpos, uint& qpos, uint size, string& destination) {
    for (int i=0; i<size; i++) {
        destination[rpos + i] = source[qpos + i];
    }
    rpos += size;
    qpos += size;
}

void add_seq(const sam_record& sam_rec, record& rec) {
    vector<uint> op_size;
    vector<char> op_type;
    parse_cigar(sam_rec.cigar, op_size, op_type);

    uint rpos = sam_rec.pos - 1;
    uint qpos = 0;
    for (int j = 0; j < op_type.size(); j++) {
        if (op_type[j] == 'M') {
            match(sam_rec.seq, rpos, qpos, op_size[j], rec.seq);
        } else if (op_type[j] == 'D') {
            rpos += op_size[j];
        } else if (op_type[j] == 'I') {
            qpos += op_size[j];
        } else if (op_type[j] == 'H' || op_type[j] == 'S' ) { // intentionally do nothing
        } else {
            cout << "Unexpected operation type " << op_type[j] << ". Allowed values are {M, D, I, H, S}.";
        }
    }
}

class sam_reader {
private:
    fstream sam;
    string current_line;
    string current_qname;

public:
    explicit sam_reader(const string& aln_file) {
        sam = fstream(aln_file, fstream::in);
    }

    bool eof() {
        return sam.eof();
    }

    void skip_header() {
        getline(sam, current_line);
        while (current_line[0] == '@') {
            getline(sam, current_line);
        }
        stringstream ss(current_line);
        ss >> current_qname;
    }

    record get_next_record(uint n) {
        record rec(n);
        string new_line;
        string new_qname;
        do {
            sam_record sam_rec(current_line);
            rec.qname = sam_rec.qname;
            add_seq(sam_rec, rec);
            getline(sam, new_line);
            stringstream ss(new_line);
            ss >> new_qname;
            current_line = new_line;
        } while (current_qname == new_qname);
        // current_line = new_line;
        current_qname = new_qname;
        return rec;
    }
};

string detect_variant(const string& seq, const map<string, string>& snps, map<string, uint>& thresholds) {
    map<string, uint> counter;
    for (const auto& snp : snps) {
        uint pos = stoi(snp.first.substr(1, snp.first.length()-1));
        char alt = snp.first[snp.first.length()-1];
        if (seq.at(pos) == alt) {
            if (counter.contains(snp.second)) {
                counter[snp.second] += 1;
            } else {
                counter[snp.second] = 1;
            }
        }
    }
    string variant = OTHER;
    for (const auto& c : counter) {
        string tmp = string(c.first);
        if (c.second >= thresholds[tmp]) {
            variant = tmp;
        }
    }
    return variant;
}

void add_counts(int ***A, const string& seq, const string& variant, const vector<string>& variants, const vector<char>& alphabet) {
    auto tmp = find(variants.begin(), variants.end(), variant);
    uint variant_num = distance(variants.begin(), tmp);
    cout << variant_num << endl;
    for (int i=0; i<seq.length(); i++) {
        auto tmp2 = find(alphabet.begin(), alphabet.end(), seq.at(i));
        uint char_num = distance(alphabet.begin(), tmp2);
        A[i][variant_num][char_num] += 1;
    }
}

int main() {
    // inputs:
    string ref_filename = "/home/andy/projects/cpp_freqtable/data/covid_ref.fa";
    string aln_filename = "/home/andy/projects/cpp_freqtable/data/covid_aln.sam";
    string mut_filename = "/home/andy/projects/cpp_freqtable/data/mut.txt";
    // outputs:
    string tsv_filename = "/home/andy/projects/cpp_freqtable/data/full_out.tsv";

    string ref = read_reference(ref_filename);
    cout << ref << endl;

    vector<string> variants;
    map<string, string> snp;        // snp: variant
    map<string, uint> threshold;    // variant: threshold
    read_mutations(mut_filename, variants, snp, threshold);
    variants.emplace_back(OTHER);

    vector<char> alphabet = {'A', 'C', 'G', 'T', 'N', '-'};

    int*** A = init_3d_array(ref.length(), variants.size(), alphabet.size());

    sam_reader sr(aln_filename);
    sr.skip_header();
    while (!sr.eof()) {
        record rec = sr.get_next_record(ref.length());
        rec.variant = detect_variant(rec.seq, snp, threshold);
        add_counts(A, rec.seq, rec.variant, variants, alphabet);
    }
    record rec = sr.get_next_record(ref.length());
    rec.variant = detect_variant(rec.seq, snp, threshold);
    add_counts(A, rec.seq, rec.variant, variants, alphabet);

    fstream out = fstream(tsv_filename, fstream::out | fstream::trunc);
    // store_3d_array(A, ref.length(), variants, alphabet, out);
    store_3d_array(A, ref.length(), variants, alphabet, cout);
    delete_3d_array(A, ref.length(), variants.size(), alphabet.size());
    return 0;
}
