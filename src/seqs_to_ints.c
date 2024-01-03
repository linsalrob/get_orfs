#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <stdint.h>
#include "colours.h"
#include "seqs_to_ints.h"
#include "error.h"

// ALternate DNA encoding lookup table
//                         A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z     
static int dnaEncodeTable [26] = {0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0};
static int dnaDecodeTable [26] = {65, 67, 71, 84};



int encode_base(int base) {
	/*
	 * Convert a base (A, G, C, T) to a number.
	 * Note, you can use seq[i]
	 *
	 * Encoding:
	 * 	A : 0 : 00
	 * 	C : 1 : 01
	 * 	G : 2 : 10 
	 * 	T : 3 : 11
	 *
	 */
    if ((base >= (int)'A') && (base <= (int)'Z')) {
        return dnaEncodeTable[base - (int)'A'];
    }
    if ((base >= (int)'a') && (base <= (int)'z')) {
		return dnaEncodeTable[base - (int)'a'];
	}

	fprintf(stderr, "We can't encode a base that is not [a-z][A-Z]. We have |%c|\n", (char) base);
	return 0;
}

int decode_base(int val) {

	if (val < 0 || val > 3) {
		fprintf(stderr, "We can't decode a value of %d\n", val);
		return 78;
	}
	return dnaDecodeTable[val];
}



uint64_t kmer_encoding(char * seq, int start_position, int k) {
	/*
	 * Given a sequence, seq, start at start_position (0 indexed), and read k characters. 
	 * Convert that to a 64-bit int using 2-bit encoding
	 *
	 * The maximum k-mer length (k-start_position is 32)
	 */

	if ((k - start_position) > 32) {
		fprintf(stderr, "%s We can only encode k<=32 strings in 64 bits. Please reduce k %s\n", PINK, ENDC);
		exit(-1);
	}

	uint64_t enc = 0;


	for (int i=start_position; i < start_position+k; i++)
		enc = (enc << 2) + encode_base(seq[i]);

	return enc;
}

char* kmer_decoding(uint64_t enc, int k) {
	/* 
	 * convert an encoded string back to a base
	 */

	char* seq = malloc(sizeof(char *) * k);
	if (!seq) {
		fprintf(stderr, "Unable to allocate memory for the sequence\n");
		exit(2);
	}
	int posn = k;
	while (posn >= 0) {
		int l = (enc & 1);
		enc >>= 1;
		int h = (enc & 1);
		enc >>= 1;
		h = (h << 1) + l;
		seq[posn--] = (char) decode_base(h);
	}
	return seq;
}



uint64_t next_kmer_encoding(char* seq, int start_position, int k, uint64_t enc) {
	/*
	 * Given a sequence, a start position, k-mer size and a previous encoding
	 * calculate the encoding that starts at start_position but remove the base
	 * at seq[start_position-1]
	 */

	if ((start_position + k - 1) > strlen(seq)) {
		fprintf(stderr, "%s Can't calculate a k-mer beyond the end of the sequence. Start: %d, k-mer %d, sequence length %ld %s\n", RED, start_position, k, strlen(seq), ENDC);
		exit(2);
	}

	enc = enc - ((uint64_t) encode_base(seq[start_position - 1]) << (2 * (k-1)));
	enc = (enc << 2) + encode_base(seq[start_position + k - 1]);
	return enc;
}

void encode_sequence(char *seq, unsigned char *enc) {
    /*
     * encode a whole sequence provided in seq and store it in enc.
     *
     * Note that this encodes the sequence to be read from left to right, unlike the above
     * kmer encoding.
     *
     * We assume that you have malloc'd and initialised enc
     *
     * unsigned char * enc = malloc(strlen(seq) * sizeof (unsigned char));
     * memset(enc, 0, strlen(seq) * sizeof (unsigned char));
     */
    for (int i =0; i<strlen(seq); i++)
        enc[i] = encode_base(seq[i]);
}

void * threaded_encode_sequence(void *thtranslate) {
    /*
     * A threaded function that writes to the encoded sequence.
     *
     * From is inclusive, so 0 to ...
     * To is exclusive so the last position needs to be strlen(seq) + 1
     */
    thtranslate_t *translate_args = (thtranslate_t *) thtranslate;
    for (int i=translate_args->from; i < translate_args->to; i++)
        translate_args->enc[i] = encode_base(translate_args->seq[i]);
    return NULL;
}

void parallel_encode_sequence(char *seq, unsigned char *enc, int num_threads) {
    /*
     * Use threading to encode the sequence.
     *
     * We need to check we encode every base :)
     */
    pthread_t threads[num_threads];
    int starts[num_threads];
    memset(starts, 0, num_threads * sizeof (*starts));
    thtranslate_t **thread_args = malloc(num_threads * sizeof(thtranslate_t *));
    if (!thread_args)
        error_and_exit("Unable to allocate memory for the threading\n");

    int fragments = strlen(seq)/num_threads;
    int from = 0;
    int to = fragments;
    int thread_number = 0;
    while (to < strlen(seq) - fragments) {
        thread_args[thread_number] = malloc(sizeof(thtranslate_t));
        thread_args[thread_number]->from = from;
        thread_args[thread_number]->to = to;
        thread_args[thread_number]->seq = seq;
        thread_args[thread_number]->enc = enc;
        starts[thread_number] = pthread_create(&threads[thread_number], NULL, &threaded_encode_sequence, (void *) thread_args[thread_number]);
        from = to;
        to += fragments;
        thread_number++;
    }
    to = strlen(seq);
    thread_args[thread_number] = malloc(sizeof(thtranslate_t));
    thread_args[thread_number]->seq = seq;
    thread_args[thread_number]->enc = enc;
    thread_args[thread_number]->from = from;
    thread_args[thread_number]->to = to;
    starts[thread_number] = pthread_create(&threads[thread_number], NULL, &threaded_encode_sequence, (void *) thread_args[thread_number]);

    // join the threads and wait for them to finish before we can return anything
    for (int thread_n = 0; thread_n <= thread_number; thread_n++) { // RAE Clean up
        pthread_join(threads[thread_n], NULL);
        free(thread_args[thread_n]);
    }

    free(thread_args);
}

void reverse_complement_sequence(unsigned char *enc, unsigned char * rc_enc, int seqlen) {
    /*
     * Reverse complement a sequence encoded above
     */

    int rclookup[4] = {3, 2,1, 0};
    for (int i = 0; i < seqlen; i++)
        rc_enc[seqlen - i - 1] = rclookup[enc[i]];
}



