//
// Created by edwa0468 on 3/01/2024.
//
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
#include "translate.h"
#include "error.h"
#include "translate_seq.h"

KSEQ_INIT(gzFile, gzread);


void help() {
    printf("USAGE: translate_seq\n");
    printf("Translate one sequence in the frame given");
    printf("-t --translation_table translation table (default=11)\n");
    printf("-s --seq the sequence\n");
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

    opt->seq = NULL;
    opt->translation_table = 11;
    opt->debug = false;
    opt->num_threads = 8;

    static struct option long_options[] = {
            {"seq",               required_argument, 0, 's'},
            {"translation_table", required_argument, 0, 't'},
            {"debug",             no_argument,       0, 'd'},
            {0, 0,                                   0, 0}
    };
    int option_index = 0;
    int gopt;
    while ((gopt = getopt_long(argc, argv, "s:t:d", long_options, &option_index)) != -1) {
        switch (gopt) {
            case 's' :
                opt->seq = strdup(optarg);
                break;
            case 't':
                opt->translation_table = atoi(optarg);
                break;
            case 'd':
                opt->debug = true;
                break;
            default:
                help();
                exit(EXIT_FAILURE);
        }
    }

    if (opt->seq == NULL) {
        fprintf(stderr, "%sPlease provide a sequence to translate%s\n", RED, ENDC);
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

    size_t seqlen = strlen(opt->seq);
    printf("Sequence length is %ld\n", seqlen);
    unsigned char *enc = malloc(seqlen * sizeof(unsigned char));
    memset(enc, 0, seqlen * sizeof(unsigned char));
    encode_sequence(opt->seq, enc);
    int protsize = (seqlen / 3) + 1;
    char *protein = malloc(protsize * sizeof *protein);
    memset(protein, 0, protsize * sizeof *protein);
    translate_one_sequence(enc, seqlen, opt->translation_table, protein);
    printf("%s\n%s\n", opt->seq, protein);
}