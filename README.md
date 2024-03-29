# profile estimation

    | Op | Description                                                                                                                                                                                                                                                                                                               |   |   |   |
    |----|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---|---|---|
    | M  | Match (alignment column containing two letters). This could contain two different letters (mismatch) or two identical letters. USEARCH generates CIGAR strings containing Ms rather than X's and ='s (see below).                                                                                                         |   |   |   |
    | D  | Deletion (gap in the target sequence).                                                                                                                                                                                                                                                                                    |   |   |   |
    | I  | Insertion (gap in the query sequence).                                                                                                                                                                                                                                                                                    |   |   |   |
    | S  | Segment of the query sequence that does not appear in the alignment. This is used with soft clipping, where the full-length query sequence is given (field 10 in the SAM record). In this case, S operations specify segments at the start and/or end of the query that do not appear in a local alignment.               |   |   |   |
    | H  | Segment of the query sequence that does not appear in the alignment. This is used with hard clipping, where only the aligned segment of the query sequences is given (field 10 in the SAM record). In this case, H operations specify segments at the start and/or end of the query that do not appear in the SAM record. |   |   |   |
    | =  | Alignment column containing two identical letters. USEARCH can read CIGAR strings using this operation, but does not generate them.                                                                                                                                                                                       |   |   |   |
    | X  | Alignment column containing a mismatch, i.e. two different letters. USEARCH can read CIGAR strings using this operation, but does not generate them.

creating input bamfile
```
minimap2 -a -x asm5 -t 8 ref.fa sequences.fasta
```

installing htslib
```
sudo apt-get update  # Ensure the package list is up to date
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
tar -xf htslib-1.13.tar.bz2
cd htslib-1.13/
autoreconf -i
./configure
make
make install
```
