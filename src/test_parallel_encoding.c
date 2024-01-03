//
// A simple test to make sure that I encode the sequence properly!
//

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "seqs_to_ints.h"
#include "test_parallel_encoding.h"


int main(int argc, char* argv[]) {
    char *seq = "ATCGATCGTCAGCATGCATCGCATCGAGCTCGTACGATCGACTAGCTACGCTACGTACGACTACGCTAGCATCGATCAGCATCACTATCGCTAGCTACGATCTAT";
    unsigned char *parallel_enc = malloc(strlen(seq) * sizeof(unsigned char));
    memset(parallel_enc, 0, strlen(seq) * sizeof(unsigned char));
    unsigned char *single_enc = malloc(strlen(seq) * sizeof(unsigned char));
    memset(single_enc, 0, strlen(seq) * sizeof(unsigned char));
    int num_threads = 8;
    parallel_encode_sequence(seq, parallel_enc, num_threads);
    encode_sequence(seq, single_enc);

    for (int i = 0; i <= strlen(seq); i++)
        printf("%c", seq[i]);
    printf("\n");

    for (int i = 0; i < strlen(seq); i++)
        printf("%d", single_enc[i]);
    printf("\n");

    for (int i = 0; i < strlen(seq); i++)
        printf("%d", parallel_enc[i]);
    printf("\n");

    for (int i = 0; i < strlen(seq); i++)
        if (parallel_enc[i] != single_enc[i])
            printf("There is a difference at position %d (parallel: %d single %d)\n", i, parallel_enc[i], single_enc[i]);
}