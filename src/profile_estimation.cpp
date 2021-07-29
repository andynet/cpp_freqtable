#include <string>
#include <htslib/sam.h>
#include <iostream>
#include <map>

#include "functions.c"

using namespace std;

typedef struct {} record;

int main(int argc, char** argv) {
    // parse_arguments
    const char *bamfile = "../data/test/subset.bam";
    string metadata = "../data/test/metadata.tsv";
    string variantsfile = "../data/test/variants.txt";
    string lineages = "../data/test/lineages.yaml";
    string output = "../data/test/output.tsv";

    char **variants = load_variants(variantsfile);

    // create_lineage_map
    map<string, string> lineage_map = create_lineage_map(lineages);

    // create_map_id2pangolin(meta)
    map<string, string> id2pangolin = create_id2pangolin(metadata);

    int ***table = init_table();

    // open(bam)
    samFile *fp = sam_open(bamfile, "r");
    sam_hdr_t *h = sam_hdr_read(fp);
    bam1_t *b = bam_init1();
    int res = sam_read1(fp, h, b);

    // while not the end:
    while (res >= 0) {
        record record = read_record(fp, h, b);
        record.variant = get_variant(record, id2pangolin, lineage_map);
        add_counts(table, record);
    }

    store(table, output);
    free(table);

//        << bam_get_qname(b) << '\n'
//        << bam_get_cigar(b) << '\n'
//        << bam_get_seq(b)   << '\n'
//        << bam_get_qual(b)  << '\n'
//        << bam_get_aux(b)   << '\n'
//        ;
        // bam_get_qname, bam_get_cigar, bam_get_seq, bam_get_qual and bam_get_aux
}
