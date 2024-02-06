
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
    printf("-t --translation_table translation table (default=11)\n");
    printf("-f --fasta fasta file\n");
    printf("-l --length minimum length (default=1 amino acid)\n");
    printf("-j --jobs number of parallel threads to use. We use 6 for the translation regardless of -j (default=8)\n");
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
    opt->num_threads = 8;

    static struct option long_options[] = {
            {"fasta", required_argument, 0, 'f'},
            {"translation_table", required_argument, 0, 't'},
            {"length", required_argument, 0, 'l'},
            {"jobs", required_argument, 0, 'j'},
            {"debug", no_argument, 0, 'd'},
            { 0, 0, 0, 0 }
    };
    int option_index = 0;
    int gopt;
    while ((gopt = getopt_long(argc, argv, "f:t:l:j:d", long_options, &option_index)) != -1) {
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
            case 'j':
                opt->num_threads = atoi(optarg);
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
        sequence->num_threads = opt->num_threads;

        // we initiate the ORFs to be 3x the sequence. They should actually be 1/3 the sequence (for codons)
        // x 6x the sequence for reading frames, but this allows a bit more.
        // We remember the size so we know when to realloc
        sequence->orf_sz = seq->seq.l * 6;
        sequence->orfs = malloc(sequence->orf_sz); // translation is x2 for fwd/rev, x3 for each frame, but /3 for amino acidds
        if (!sequence->orfs)
            error_and_exit("Unable to allocate memory for the sequence ORFs\n");
        memset(sequence->orfs, 0, sequence->orf_sz);


        // We allocate a MB for the ORF names, but then if we need more we realloc it dynamically, later.
        // We remember the size so we know when to realloc
        sequence->orf_name_sz = 1000000;
        sequence->orf_names = malloc(sequence->orf_name_sz); // somewhat randomly allocate the same size of the memory as above
        memset(sequence->orf_names, 0, sequence->orf_name_sz);

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


