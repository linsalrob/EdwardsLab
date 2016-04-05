/*
 bioinf.c
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "string_buffer.h"

#include "bioinf.h"


#define INIT_NUM_OF_CHROMS 50
#define INIT_CHROM_NAME_LEN 80
// Ensure capacity for 260Mbp (largest human/chimp chr is ~250Mbp)
#define INIT_CHROM_SEQ_LEN 260000000
#define is_base_char(x) ((x) == 'a' || (x) == 'A' || \
                         (x) == 'c' || (x) == 'C' || \
                         (x) == 'g' || (x) == 'G' || \
                         (x) == 't' || (x) == 'T' || \
                         (x) == 'n' || (x) == 'N')

struct SEQ_FILE
{
  gzFile gz_file;
  FILE* file;
  enum SEQ_FILE_TYPE file_type;
};

void _set_seq_filetype(SEQ_FILE* file)
{
  int first_char;

  do
  {
    first_char = (file->file != NULL) ? getc(file->file)
                                      : gzgetc(file->gz_file);
  } while (first_char != -1 && (first_char == '\n' || first_char == '\r'));

  if(first_char == -1)
  {
    //fprintf(stderr, "Warning: empty sequence file\n");
    return;
  }
  else if(first_char == '>')
  {
    // Reading FASTA
    file->file_type = SEQ_FASTA;
  }
  else if(first_char == '@')
  {
    // Reading FASTQ
    file->file_type = SEQ_FASTQ;
  }
  else {
    file->file_type = SEQ_PLAIN;
  }

  // Put char back
  if(file->file != NULL)
  {
    ungetc(first_char, file->file);
  }
  else
  {
    gzungetc(first_char, file->gz_file);
  }
}

SEQ_FILE* seq_file_gzinit(gzFile gz_file)
{
  SEQ_FILE* seq_file = (SEQ_FILE*) malloc(sizeof(SEQ_FILE));
  seq_file->gz_file = gz_file;
  seq_file->file = NULL;
  seq_file->file_type = SEQ_UNKNOWN;
  
  _set_seq_filetype(seq_file);
  
  return seq_file;
}

SEQ_FILE* seq_file_init(FILE* file)
{
  SEQ_FILE* seq_file = (SEQ_FILE*) malloc(sizeof(SEQ_FILE));
  seq_file->gz_file = NULL;
  seq_file->file = file;
  seq_file->file_type = SEQ_UNKNOWN;

  _set_seq_filetype(seq_file);

  return seq_file;
}

SEQ_FILE* seq_file_open(char* path, char* mode)
{
  FILE* file = fopen(path, mode);

  if(file == NULL)
  {
    fprintf(stderr, "Error: Cannot open file '%s' with mode '%s'\n", path, mode);
    exit(EXIT_FAILURE);
  }

  return seq_file_init(file);
}

SEQ_FILE* seq_file_gzopen(char* path, char* mode)
{
  gzFile gz_file = gzopen(path, mode);

  if(gz_file == NULL)
  {
    fprintf(stderr, "Error: Cannot open file '%s' with mode '%s'\n", path, mode);
    exit(EXIT_FAILURE);
  }

  return seq_file_gzinit(gz_file);
}

void seq_file_close(SEQ_FILE* seq)
{
  if(seq->gz_file != NULL)
  {
    gzclose(seq->gz_file);
  }
  if(seq->file != NULL)
  {
    fclose(seq->file);
  }
}

// Open using a command line argument
SEQ_FILE* seq_open_cmd_arg(char* cmd_arg)
{
  if(cmd_arg == NULL)
  {
    return NULL;
  }

  if(strcmp(cmd_arg, "-") == 0)
  {
    // Read from STDIN
    return seq_file_init(stdin);
  }
  else
  {
    return seq_file_gzopen(cmd_arg, "r");
  }
}

// cmd is the command used to open
void seq_close_cmd_arg(SEQ_FILE* seq, char* cmd_arg)
{
  if(seq == NULL || cmd_arg == NULL)
  {
    return;
  }

  if(strcmp(cmd_arg, "-") != 0)
  {
    // If we were not reading from stdin
    // (don't want to try to close stdin)
    seq_file_close(seq);
  }

  seq_file_free(seq);
}

// still need to free file->file / file->gz_file after this
void seq_file_free(SEQ_FILE* seq)
{
  free(seq);
}

enum SEQ_FILE_TYPE seq_file_get_type(SEQ_FILE* file)
{
  return file->file_type;
}

void seq_file_read(SEQ_FILE* file, STRING_BUFFER* title, STRING_BUFFER* sequence)
{
  string_buff_reset(title);
  string_buff_reset(sequence);

  if(file->file_type == SEQ_FASTA)
  {
    read_fasta_entry(title, sequence, file->file, file->gz_file);
  }
  else if(file->file_type == SEQ_FASTQ)
  {
    read_fastq_entry(title, sequence, NULL, file->file, file->gz_file);
  }
  else if(file->file_type == SEQ_PLAIN)
  {
    t_buf_pos chars_read;

    do
    {
      chars_read
        = (file->file != NULL) ? string_buff_readline(sequence, file->file)
                               : string_buff_gzreadline(sequence, file->gz_file);

      string_buff_chomp(sequence);
    } while(chars_read > 0 && string_buff_strlen(sequence) == 0);
  }
}

// Returns 1 if a header is read, 0 otherwise
// Prints errors to stderr if there are syntax issues
char read_fastq_entry(STRING_BUFFER* header, STRING_BUFFER* sequence,
                      STRING_BUFFER* quality, FILE *file, gzFile gz_file)
{
  string_buff_reset(header);
  string_buff_reset(sequence);

  // Read until we have a header line
  int read_len;

  do
  {
    read_len = (file != NULL) ? string_buff_readline(header, file)
                              : string_buff_gzreadline(header, gz_file);

    string_buff_chomp(header);
  }
  while (read_len > 0 && string_buff_strlen(header) == 0);

  if(read_len == 0)
  {
    return 0;
  }
  else if(header->buff[0] != '@')
  {
    fprintf(stderr, "Warning: FASTQ header does not begin with 'Q' ('%s')\n",
            header->buff);
  }

  // Read sequence
  read_len = (file != NULL) ? string_buff_readline(sequence, file)
                            : string_buff_gzreadline(sequence, gz_file);

  string_buff_chomp(sequence);

  if(read_len == 0)
  {
    fprintf(stderr, "Warning: FASTQ missing sequence line (header: '%s')\n",
            header->buff);
  }

  // Read skip line ('+')
  if(quality != NULL)
  {
    string_buff_reset(quality);

    // Use quality buffer as temp buff to read in skip line
    read_len = (file != NULL) ? string_buff_readline(quality, file)
                              : string_buff_gzreadline(quality, gz_file);

    if(read_len == 0 || string_buff_get_char(quality, 0) != '+')
    {
      string_buff_chomp(quality);
      fprintf(stderr, "Warning: FASTQ skip line does not begin with '+' "
                      "(header: '%s', skip line: '%s')\n",
              header->buff, quality->buff);
    }

    // Now read quality line
    read_len = (file != NULL) ? string_buff_readline(quality, file)
                              : string_buff_gzreadline(quality, gz_file);

    string_buff_chomp(quality);

    if(read_len == 0)
    {
      fprintf(stderr, "Warning: FASTQ is missing a quality line\n");
    }
    else if(string_buff_strlen(quality) != string_buff_strlen(sequence))
    {
      fprintf(stderr, "Warning: FASTQ sequence and quality lines are not "
                      "the same length (header: '%s')\n", header->buff);
      fprintf(stderr, "Sequence: '%s'\n", sequence->buff);
      fprintf(stderr, "Quality : '%s'\n", quality->buff);
    }
  }
  else
  {
    // Just skip two lines
    if(file != NULL)
    {
      string_buff_skip_line(file);
      string_buff_skip_line(file);
    }
    else
    {
      string_buff_gzskip_line(gz_file);
      string_buff_gzskip_line(gz_file);
    }
  }

  return 1;
}

// Read an entry from a FASTA file
char read_fasta_entry(STRING_BUFFER* header, STRING_BUFFER* sequence,
                      FILE *file, gzFile gz_file)
{
  string_buff_reset(header);
  string_buff_reset(sequence);

  // Read until we have a header line
  int read_len;

  do
  {
    read_len = (file != NULL) ? string_buff_readline(header, file)
                              : string_buff_gzreadline(header, gz_file);

    string_buff_chomp(header);
  }
  while (read_len > 0 && string_buff_strlen(header) == 0);

  if(read_len == 0)
  {
    return 0;
  }
  else if(header->buff[0] != '>')
  {
    fprintf(stderr, "Warning: FASTA header does not begin with '>' ('%s')\n",
            header->buff);
  }

  char success = 0;

  while(1)
  {
    // Check line doesn't begin with '>'
    int first_char = (file != NULL) ? getc(file) : gzgetc(gz_file);
  
    if(first_char == -1)
    {
      break;
    }
    else if(first_char == '>')
    {
      // Push char back onto buffer
      if(file != NULL)
      {
        ungetc(first_char, file);
      }
      else
      {
        gzungetc(first_char, gz_file);
      }

      // Done
      break;
    }
    else
    {
      // Push char onto string
      string_buff_append_char(sequence, first_char);
      string_buff_chomp(sequence);
      success = 1;
    
      if(first_char != '\n' && first_char != '\r')
      {
        // Read the rest of the line
        t_buf_pos chars_read;

        chars_read = (file != NULL) ? string_buff_readline(sequence, file)
                                    : string_buff_gzreadline(sequence, gz_file);

        if(chars_read == 0)
        {
          break;
        }

        string_buff_chomp(sequence);
      }
    }
  }

  return success;
}
















// Load each chromosome's name and sequence into respective arrays, passed as
//  parameters
// returns number of sequences loaded
//         or -1 on error
long load_fasta_file(const char* ref_genome_file,
                     gzFile ref_genome_file_handle,
                     STRING_BUFFER*** chrom_names_ptr,
                     STRING_BUFFER*** chrom_seqs_ptr)
{
  int array_mem_size = INIT_NUM_OF_CHROMS*sizeof(STRING_BUFFER*);
  STRING_BUFFER** chrom_names = (STRING_BUFFER**) malloc(array_mem_size);
  STRING_BUFFER** chrom_seqs = (STRING_BUFFER**) malloc(array_mem_size);

  int chrom_array_capacity = INIT_NUM_OF_CHROMS;
  int chrom_index = 0;
  
  STRING_BUFFER* chrom_name = string_buff_init(INIT_CHROM_NAME_LEN);
  STRING_BUFFER* chrom_seq = string_buff_init(INIT_CHROM_SEQ_LEN);
  
  //fgets(line_init_buf, REF_GEONME_BUF_LEN, ref_genome_file) != NULL
  //int first_char;
  // first char and \0 only
  char* first_char_str = (char*) malloc(sizeof(char)*2);
  
  while(gzread(ref_genome_file_handle, first_char_str, 1) > 0)
  {
    char first_char = first_char_str[0];
  
    if(first_char == '>')
    {
      // Store existing chrom_name / chrom_seq
      if(chrom_name->len > 0 || chrom_seq->len > 0) {
        if(chrom_name->len == 0) {
          fprintf(stderr, "Error: FASTA seqence without name in file '%s'\n",
                  ref_genome_file);
          return -1;
        }
        
        if(chrom_seq->len == 0) {
          fprintf(stderr, "Error: FASTA name without sequence in file '%s': %s\n",
                  ref_genome_file, chrom_name->buff);
          return -1;
        }

        // Convert sequence to all upper case to simplify code later
        string_buff_to_uppercase(chrom_seq);

        chrom_names[chrom_index] = chrom_name;
        chrom_seqs[chrom_index] = chrom_seq;

        chrom_index++;
        
        if(chrom_index == chrom_array_capacity)
        {
          // extend arrays
          chrom_array_capacity *= 2;
          
          // number of bytes required to store each of the arrays
          array_mem_size = chrom_array_capacity*sizeof(STRING_BUFFER*);
          
          chrom_names = realloc(chrom_names, array_mem_size);
          chrom_seqs = realloc(chrom_seqs, array_mem_size);
        }
        
        // re-initialise the STRING_BUFFER objects
        // (using final size of previous buffer to chose the size of the new one)
        chrom_name = string_buff_init(chrom_name->size);
        chrom_seq = string_buff_init(chrom_seq->size);
      }
      
      // load rest of line into chrom_name
      if(string_buff_gzreadline(chrom_name, ref_genome_file_handle) == 0)
      {
        fprintf(stderr, "Error: Unexpected ref fasta End Of File '%s' "
                        "(in title line)\n",
                ref_genome_file);
        return -1;
      }

      // remove newline characters
      string_buff_chomp(chrom_name);
      
      // printf("title '%s'\n", chrom_name->buff);
    }
    else if(is_base_char(first_char))
    {
      // append char to chrom_seq
      if(chrom_seq->len+1 >= chrom_seq->size)
      {
        // Increase buffer size
        string_buff_resize(chrom_seq, 2*chrom_seq->size);
      }

      // copy first character to string buffer
      string_buff_append_char(chrom_seq, first_char);

      // load rest of line into chrom_seq
      //  (might not be any more characters to read)
      string_buff_gzreadline(chrom_seq, ref_genome_file_handle);

      // remove newline characters
      string_buff_chomp(chrom_seq);
    
      // printf("seq '%s'\n", chrom_seq->buff);
    }
    else if(first_char != '\n') {
      fprintf(stderr, "Error: Unexpected line start character('%c') "
                      "in FASTA file '%s'\n",
             first_char, ref_genome_file);
      return -1;
    }
  }
  
  free(first_char_str);
  first_char_str = NULL;
  
  // Store last chrom
  if(chrom_name->len > 0 || chrom_seq->len > 0)
  {
    if(chrom_name->len == 0)
    {
      fprintf(stderr, "Error: FASTA seqence without name in file '%s'\n",
              ref_genome_file);
      return -1;
    }
        
    if(chrom_seq->len == 0) {
      fprintf(stderr, "Error: FASTA name without sequence in file '%s': %s\n",
              ref_genome_file, chrom_name->buff);
      return -1;
    }
    
    // Convert sequence to all upper case to simplify code later
    string_buff_to_uppercase(chrom_seq);

    chrom_names[chrom_index] = chrom_name;
    chrom_seqs[chrom_index] = chrom_seq;

    chrom_index++;
  }
  else
  {
    // Discard remaining (unused) string buffers
    string_buff_free(chrom_name);
    string_buff_free(chrom_seq);
  }

  gzclose(ref_genome_file_handle);

  // Return arrays passed as parameters
  *chrom_names_ptr = chrom_names;
  *chrom_seqs_ptr = chrom_seqs;

  return chrom_index;
}


// Sort chromosome list by chromosome names
void sort_chroms(STRING_BUFFER*** chrom_names,
                 STRING_BUFFER*** chrom_seqs,
                 const int num_of_chroms)
{
  // Shrink and sort chrom_names, chrom_seqs
  int chroms_order[num_of_chroms];
  
  // Create an initial ordering
  int i;
  for(i = 0; i < num_of_chroms; i++)
  {
    chroms_order[i] = i;
  }
  
  // define sort compare function
  //  - sort by string ordering of chromosome names
  int cmp_chroms(const void *a, const void *b)
  {
    const int *ia = (const int *) a;
    const int *ib = (const int *) b;
    
    return strcmp(((*chrom_names)[*ia])->buff, ((*chrom_names)[*ib])->buff);
  }
  
  qsort(chroms_order, num_of_chroms, sizeof(int), cmp_chroms);
  
  /* // Debug: print chromosome's order
  printf("%i", chroms_order[0]);
  for(i = 1; i < num_of_chroms; i++)
  {
    printf(",%i", chroms_order[i]);
  }
  printf("\n");
  */

  int shrunk_array_mem_size = num_of_chroms * sizeof(STRING_BUFFER*);
  STRING_BUFFER** temp_chrom_names = (STRING_BUFFER**) malloc(shrunk_array_mem_size);
  STRING_BUFFER** temp_chrom_seqs = (STRING_BUFFER**) malloc(shrunk_array_mem_size);

  for(i = 0; i < num_of_chroms; i++)
  {
    temp_chrom_names[i] = (*chrom_names)[chroms_order[i]];
    temp_chrom_seqs[i] = (*chrom_seqs)[chroms_order[i]];
  }

  // free old arrays
  free(*chrom_names);
  free(*chrom_seqs);
  
  // re-assign
  *chrom_names = temp_chrom_names;
  *chrom_seqs = temp_chrom_seqs;
}
