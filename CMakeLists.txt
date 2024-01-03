cmake_minimum_required(VERSION 3.16)
project(get_orfs VERSION 0.2 LANGUAGES C)

set(CMAKE_C_STANDARD 99)

add_executable(test_parallel_encoding "")
add_executable(get_orfs "")
add_subdirectory(src)
target_link_libraries(get_orfs z pthread)
install (TARGETS get_orfs DESTINATION bin)


target_link_libraries(test_parallel_encoding z pthread)
install (TARGETS test_parallel_encoding DESTINATION bin)
