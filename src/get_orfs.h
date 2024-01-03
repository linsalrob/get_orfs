//
// Created by edwa0468 on 29/12/2023.
//

#ifndef CODONS_GET_ORFS
#define CODONS_GET_ORFS


struct options {
    char *fasta_file;
    int translation_table;
    bool debug;
    int minlen;
    int num_threads;
};



#define  __version__ 0.1

#endif //CODONS_GET_ORFS
