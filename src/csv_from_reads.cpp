#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <map>
#include <cstdlib>
#include <cstring>

#include <getopt.h>

#define OTHER "other"

using namespace std;

class sam_record {
public:
    string qname, cigar, seq;
    int pos = 0;

    explicit sam_record(const string& line){
        stringstream ss(line);
        string del;
        //    qname    flag   rname  pos    mapq   cigar    rnext  pnext  tlen   seq    qual tags
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

void add_seq(const sam_record& sam_rec, record& rec);

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
        if (current_line.empty()) {return record();}
        string new_line;
        record rec(n);
        string new_qname;
        do {
            sam_record sam_rec(current_line);
            rec.qname = sam_rec.qname;
            add_seq(sam_rec, rec);
            getline(sam, new_line);
            if (new_line.empty()) {
                new_qname = "";
            } else {
                stringstream ss(new_line);
                ss >> new_qname;
            }
            current_line = new_line;
        } while (current_qname == new_qname);
        current_qname = new_qname;
        return rec;
    }
};

string read_reference(const string& filename) {
    fstream file(filename, fstream::in);
    if (!file.is_open()) {
        cout << "File " << filename << " could not be open. Does it exist?\n";
        exit(EXIT_FAILURE);
    }

    string ref;
    string tmp;
    while (!file.eof()) {
        getline(file, tmp);
        if (tmp[0] == '>') {
            // skip
        } else {
            ref += tmp;
        }
    }
    file.close();
    return ref;
}

int*** init_3d_array(uint x, uint y, uint z) {
    int*** res = new int**[x];
    for (uint i=0; i<x; i++) {
        res[i] = new int*[y];
        for (uint j=0; j<y; j++) {
            res[i][j] = new int[z];
            for (uint k=0; k<z; k++) {
                res[i][j][k] = 0;
            }
        }
    }
    return res;
}

void store_3d_array(int ***A, uint length, const vector<string>& variants, const vector<char>& alphabet, const string& out_filename) {
    fstream out = fstream(out_filename, fstream::out | fstream::trunc);
    cout << "Storing table to " << out_filename << endl;
    for (uint i=0; i<variants.size(); i++) {
        for (uint j=0; j<length; j++) {
            for (uint k=0; k<alphabet.size(); k++) {
                out << variants.at(i) << '\t' << j << '\t' << alphabet.at(k) << '\t' << A[j][i][k] << endl;
            }
        }
    }
    cout << "Table stored successfully." << endl;
}

void delete_3d_array(int ***A, uint x, uint y) {
    for (uint i=0; i<x; i++) {
        for (uint j=0; j<y; j++) {
            delete[] A[i][j];
        }
        delete[] A[i];
    }
    delete[] A;
}

string trim(const string& line) {
    const char* white_space = " \t\v\r\n";
    size_t start = line.find_first_not_of(white_space);
    size_t end = line.find_last_not_of(white_space);
    return start == end ? string() : line.substr(start, end - start + 1);
}

void read_mutations(const string& filename, vector<string>& variants, multimap<string, string>& snps, map<string, uint>& thresholds) {
    fstream file(filename, fstream::in);
    string tmp;
    string variant;
    uint thr;
    while (!file.eof()) {
        getline(file, tmp);
        tmp = trim(tmp);
        stringstream ss(tmp);
        ss >> variant >> thr;
        variants.push_back(variant);
        thresholds.insert(pair<string, uint>(variant, thr));
        string snp;
        while (!ss.eof()) {
            ss >> snp;
            snps.insert(pair<string, string>(snp, variant));
        }
    }
}

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
    for (uint i=0; i<size; i++) {
        destination[rpos + i] = source[qpos + i];
    }
    rpos += size;
    qpos += size;
}

void deletion(uint& rpos, uint size, string& destination) {
    for (uint i=0; i<size; i++) {
        destination[rpos + i] = '-';
    }
    rpos += size;
}

void add_seq(const sam_record& sam_rec, record& rec) {
    vector<uint> op_size;
    vector<char> op_type;
    parse_cigar(sam_rec.cigar, op_size, op_type);

    uint rpos = sam_rec.pos - 1;
    uint qpos = 0;
    for (uint j = 0; j < op_type.size(); j++) {
        if (op_type[j] == 'M') {
            match(sam_rec.seq, rpos, qpos, op_size[j], rec.seq);
        } else if (op_type[j] == 'D') {
            deletion(rpos, op_size[j], rec.seq);
        } else if (op_type[j] == 'I' || op_type[j] == 'S') {
            qpos += op_size[j];
        } else if (op_type[j] == 'H' ) {
            // intentionally do nothing
        } else {
            // cout << "Unexpected operation type " << op_type[j] << ". Allowed values are {M, D, I, H, S}." << endl;
            rec.seq = "";
            return;
        }
    }
}

uint get_pos(const string& s) {
    return stoi(s.substr(1, s.length()-1)) - 1;
}

void add_counts(int ***A, const string& seq, const string& variant, const vector<string>& variants, const vector<char>& alphabet) {
    auto tmp = find(variants.begin(), variants.end(), variant);
    uint variant_num = distance(variants.begin(), tmp);
    for (uint i=0; i<seq.length(); i++) {
        auto tmp2 = find(alphabet.begin(), alphabet.end(), seq.at(i));
        uint char_num;
        char_num = distance(alphabet.begin(), tmp2);
        A[i][variant_num][char_num] += 1;
    }
}

void check_sanity(const string& reference, const multimap<string, string>& snp) {
    for (const auto& s : snp) {
        char ref = s.first[0];
        uint pos = get_pos(s.first);
        if (reference[pos] != ref) {
            cout << "reference[" << pos << "] is not " << ref << endl;
        }
    }
}

void correct_alphabet(string& seq, const vector<char>& alphabet, const char c) {
    for (char & i : seq) {
        i = (char)toupper(i);
        int equal = 0;
        for (char j : alphabet) {
            if (i == j) {
                equal++;
            }
        }
        if (equal == 0) {
            i = c;
        }
    }
}

typedef enum {REF, MUT, ALN, VAR, OUT} arg;
void parse_args(int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, "r:m:a:v:o:")) != -1) {
        switch (opt) {
            case 'r': argv[REF] = optarg; break;
            case 'm': argv[MUT] = optarg; break;
            case 'a': argv[ALN] = optarg; break;
            case 'v': argv[VAR] = optarg; break;
            case 'o': argv[OUT] = optarg; break;
            default:
                cerr << "Usage: " << argv[0] << " -r ref.fasta -m mut.txt -a aln.sam -o out.tsv" << endl;
                exit(EXIT_FAILURE);
        }
    }
}

int main(int argc, char **argv) {
    // run with:
    // ./cmake-build-debug/cpp_freqtable
    // -r data/covid_ref.fa -m data/mut.txt
    // -a data/msa-minimap.sam -v B.1.1.7 -o data/full_out.tsv

    parse_args(argc, argv);
    string ref_filename = argv[REF];    // can be multiline fasta
    string aln_filename = argv[ALN];    // sam has to end with empty line
    string mut_filename = argv[MUT];    // can end with whitespace, has to end with empty line
    string variant = argv[VAR];
    string out_filename = argv[OUT];

    string ref = read_reference(ref_filename);

    vector<string> variants;
    multimap<string, string> snp;   // snp: variant
    map<string, uint> threshold;    // variant: threshold
    read_mutations(mut_filename, variants, snp, threshold);
    variants.emplace_back(OTHER);

    check_sanity(ref, snp);

    vector<char> alphabet = {'A', 'C', 'G', 'T', 'N', '-'};

    int*** A = init_3d_array(ref.length(), variants.size(), alphabet.size());

    sam_reader sr(aln_filename);
    sr.skip_header();
    for (uint i=0; !sr.eof(); i++) {
        record rec = sr.get_next_record(ref.length());
        correct_alphabet(rec.seq, alphabet, 'N');
        add_counts(A, rec.seq, variant, variants, alphabet);
        if (i % 1000 == 0){ cout << i << " records processed." << endl;}
    }
    cout << "All records processed." << endl;

    store_3d_array(A, ref.length(), variants, alphabet, out_filename);
    delete_3d_array(A, ref.length(), variants.size());
    return 0;
}
