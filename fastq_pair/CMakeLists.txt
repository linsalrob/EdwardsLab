cmake_minimum_required(VERSION 3.8)
project(fastq_pair)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES main.c robstr.c index_fastq.c)
add_executable(fastq_pair ${SOURCE_FILES})