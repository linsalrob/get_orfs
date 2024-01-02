
/*
 * Translate a DNA sequence in all 6 frames using a translation table of your choice
 * This is a refactored version of translate_sequence.c. That code works, so I don't want to break it :)
 */

#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>
#include <getopt.h>

#include "kseq.h"

#include "colours.h"
#include "seqs_to_ints.h"
#include "get_orfs.h"
#include "translate.h"
#include "error.h"

KSEQ_INIT(gzFile, gzread);


void help() {
    printf("USAGE: get_orfs\n");
    printf("-t --translation_table <translation table (default=11)\n");
    printf("-f --fasta fasta file\n");
    printf("-l --length minimum length (default=1 amino acid)\n");
    printf("-d --debug print debugging help");
    printf("-v --version print the version and exit\n");
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        help();
        exit(0);
    }

    if (argc == 2 && ((strcmp(argv[1], "-v") == 0) || (strcmp(argv[1], "--version") == 0))) {
        printf("%s version: %f\n", argv[0], __version__);
        exit(0);
    }

    struct options *opt;
    opt = malloc(sizeof(struct options));
    if (!opt) {
        fprintf(stderr, "Unable to allocate memory for the options\n");
        exit(2);
    }

    opt->fasta_file = NULL;
    opt->translation_table = 11;
    opt->debug = false;
    opt->minlen = 1;

    static struct option long_options[] = {
            {"fasta", required_argument, 0, 'f'},
            {"translation_table", required_argument, 0, 't'},
            {"length", required_argument, 0, 'l'},
            {"debug", no_argument, 0, 'd'},
            { 0, 0, 0, 0 }
    };
    int option_index = 0;
    int gopt;
    while ((gopt = getopt_long(argc, argv, "f:t:l:d", long_options, &option_index)) != -1) {
        switch (gopt) {
            case 'f' :
                opt->fasta_file = strdup(optarg);
                break;
            case 't':
                opt->translation_table = atoi(optarg);
                break;
            case 'l':
                opt->minlen = atoi(optarg);
                break;
            case 'd':
                opt->debug = true;
                break;
            default:
                help();
                exit(EXIT_FAILURE);
        }
    }

    if (opt->fasta_file == NULL) {
        fprintf(stderr, "%sPlease provide a fasta file%s\n", RED, ENDC);
        help();
        exit(EXIT_FAILURE);
    }

    // invalid codon tables: 0, 7, 8, 17, 18, 19, 20,  > 31
    if (opt->translation_table < 0 || opt->translation_table > 31 ||
        opt->translation_table == 7 || opt->translation_table == 8 ||
        (opt->translation_table > 16 && opt->translation_table < 21)) {
        fprintf(stderr, "%sTranslation table %d is not valid. It can not be <0, 7, 8, 17-20, or >31%s\n",
                RED, opt->translation_table, ENDC);
        help();
        exit(EXIT_FAILURE);
    }


    gzFile fp = gzopen(opt->fasta_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->fasta_file, ENDC);
        exit(3);
    }
    kseq_t *seq = kseq_init(fp);
    int l;
    while ((l = kseq_read(seq)) >= 0) {

        translate_t *sequence = malloc(sizeof (translate_t));
        if (!sequence)
            error_and_exit("Unable to allocate memory for the sequence\n");

        sequence->dnaseq = malloc((seq->seq.l + 1) * sizeof (char));
        if (!sequence->dnaseq)
            error_and_exit("Unable to allocate memory for the DNA sequence\n");

        strcpy(sequence->dnaseq, seq->seq.s);
        sequence->name = malloc((strlen(seq->name.s) + 1) * sizeof (char));
        if (!sequence->name)
            error_and_exit("Unable to allocate memory for the sequence name\n");

        strcpy(sequence->name, seq->name.s);
        sequence->len = seq->seq.l;
        sequence->translation_table = opt->translation_table;
        sequence->verbose = opt->debug;
        sequence->num_orfs = 0;
        sequence->orfs = malloc(seq->seq.l * 3 * sizeof (char)); // translation is x2 for fwd/rev, x3 for each frame, but /3 for amino acidds
        if (!sequence->orfs)
            error_and_exit("Unable to allocate memory for the sequence ORFs\n");


        // Here I tried to dynamically figure out how much memory we might need, but I gave up and just set it to a million
        // THere are definitely better ways to do this.
        // the sequence name is "seq_name frame +x start stop"
        // assuming sequence length less than 1e9 (10 characters each for start/stop), we need 20 + seq_name + 7 for " frame " + 2 for [+/-]x + 2 spaces
        // * the number of ORFs we find. Assuming most orfs are 1 per kb of the sequence should be about seqlen /300
        //int seqnamesize = ((31 + strlen(seq->name.s)) * (seq->seq.l/10)) * sizeof (char);
        size_t seqnamesize = 1000000;

        if (opt->debug)
            fprintf(stderr, "%sWARN: Estimated we might need %ld bytes for ORF names%s\n", YELLOW, seqnamesize, ENDC);
        sequence->orf_names = malloc(seqnamesize); // somewhat randomly allocate the same size of the memory as above

        if (!sequence->orf_names)
            error_and_exit("Unable to allocate memory for the sequence ORF names\n");


        parallel_translate(sequence);

        int orfcount = 0;
        for (int orf=0; orf<sequence->num_orfs; orf++)
            if (strlen(sequence->orfs[orf]) > opt->minlen)
                printf(">orf%d [%s]\n%s\n", ++orfcount, sequence->orf_names[orf], sequence->orfs[orf]);

        for (int orf=0; orf<sequence->num_orfs; orf++) {
            free(sequence->orfs[orf]);
            free(sequence->orf_names[orf]);
        }
        free(sequence->orfs);
        free(sequence->orf_names);
        free(sequence->dnaseq);
        free(sequence->name);
/*
        for (int i=0; i < sequence->num_orfs; i++) {
            free(sequence->orf_names[i]);
            free(sequence->orfs[i]); // valgrind is complaining about this free??
        }
*/
        free(sequence);
    }
    free(opt->fasta_file);
    free(opt);
    kseq_destroy(seq);
    gzclose(fp);
}


