cmake_minimum_required(VERSION 3.16)
project(get_orfs VERSION 0.2)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES get_orfs.c get_orfs.h colours.h error.c error.h kseq.h seqs_to_ints.c seqs_to_ints.h translate.c translate.h)

add_executable(get_orfs ${SOURCE_FILES})
target_link_libraries(get_orfs z pthread)
install (TARGETS get_orfs DESTINATION bin)

