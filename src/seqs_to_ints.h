
#ifndef SEQS_TO_INT_H
#define SEQS_TO_INT_H

#include <stdint.h>

/*
 * convert a kmer of A,T,G,C into a uint64_t
 */
uint64_t kmer_encoding(char*, int, int);

/*
 * decode a uint64_t into a kmer of A,T,G,C. We need to know how long k is otherwise it will be left filled with A's!
 */

char* kmer_decoding(uint64_t, int);

/*
 * find the next primer for a sequence
 */
uint64_t next_kmer_encoding(char*, int, int, uint64_t);

int encode_base(int);
int decode_base(int);
void encode_sequence(char *, unsigned char *);
void parallel_encode_sequence(char *, unsigned char *, int);
void * threaded_encode_sequence(void *);
void reverse_complement_sequence(unsigned char*, unsigned char *, int);

typedef struct thtranslate_struct {
    char *seq; // the DNA sequence
    unsigned char *enc; // where we are going to save the encoding
    int from; // where we are going to start from (inclusive)
    int to; // where we are going to stop (exclusive)
} thtranslate_t;


#endif
