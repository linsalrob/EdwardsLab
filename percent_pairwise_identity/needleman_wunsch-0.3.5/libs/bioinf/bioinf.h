/*
 bioinf.h
 project: bioinf
 author: Isaac Turner <isaac.turner@dtc.ox.ac.uk>
 Copyright (C) 1-Nov-2011
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BIOINF_HEADER_SEEN
#define BIOINF_HEADER_SEEN

#include <stdio.h>
#include <zlib.h>

typedef struct SEQ_FILE SEQ_FILE;

enum SEQ_FILE_TYPE { SEQ_UNKNOWN,SEQ_FASTA,SEQ_FASTQ,SEQ_PLAIN };

SEQ_FILE* seq_file_gzopen(char* path, char* mode);
SEQ_FILE* seq_file_open(char* path, char* mode);
SEQ_FILE* seq_file_init(FILE* file);
SEQ_FILE* seq_file_gzinit(gzFile file);
SEQ_FILE* seq_open_cmd_arg(char* cmd);
void seq_close_cmd_arg(SEQ_FILE* seq, char* cmd);
void seq_file_close(SEQ_FILE* seq);
void seq_file_free(SEQ_FILE* file);
void seq_file_read(SEQ_FILE* file, STRING_BUFFER* title, STRING_BUFFER* sequence);
enum SEQ_FILE_TYPE seq_file_get_type(SEQ_FILE* file);

// Pass either FILE* or gzFile with the other ptr NULL
char read_fasta_entry(STRING_BUFFER* title, STRING_BUFFER* sequence,
                      FILE *file, gzFile gz_file);

// Pass either FILE* or gzFile with the other ptr NULL
// If quality is NULL it is not stored
char read_fastq_entry(STRING_BUFFER* header, STRING_BUFFER* sequence,
                      STRING_BUFFER* quality, FILE *file, gzFile gz_file);

// Other stuff
long load_fasta_file(const char* ref_genome_file,
                     gzFile ref_genome_file_handle,
                     STRING_BUFFER*** chrom_names_ptr,
                     STRING_BUFFER*** chrom_seqs_ptr);

void sort_chroms(STRING_BUFFER*** chrom_names,
                 STRING_BUFFER*** chrom_seqs,
                 const int num_of_chroms);


// SAM file waiting here
typedef struct SAM_ENTRY
{
  char *qname;
  int flag;
  char *rname;
  int pos, mapq;
  char *cigar, *rnext;
  int pnext, tlen;
  char *seq, *qual;
} SAM_ENTRY;

SAM_ENTRY* sam_alloc();

void sam_entry_alloc(SAM_ENTRY* sam_entry);
void sam_entry_free(SAM_ENTRY* sam_entry);

char sam_read(FILE* sam_file, SAM_ENTRY* sam_entry);

#endif
