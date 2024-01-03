set(LIBRARY_FILES  colours.h error.c error.h kseq.h seqs_to_ints.c seqs_to_ints.h translate.c translate.h)

target_sources(get_orfs PUBLIC ${LIBRARY_FILES} get_orfs.c get_orfs.h)
target_sources(test_parallel_encoding PRIVATE ${LIBRARY_FILES} test_parallel_encoding.c test_parallel_encoding.h)