/*
 string_buffer.h
 project: string_buffer
 url: https://github.com/noporpoise/StringBuffer
 author: Isaac Turner <turner.isaac@gmail.com>

 Copyright (c) 2011, Isaac Turner
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef STRING_BUFFER_FILE_SEEN
#define STRING_BUFFER_FILE_SEEN

#include <zlib.h> // needed for gzFile
#include <stdarg.h> // required for va_list

typedef unsigned long t_buf_pos;
typedef struct STRING_BUFFER STRING_BUFFER;

struct STRING_BUFFER
{
  char *buff;
  t_buf_pos len; // length of the string
  t_buf_pos size; // buffer size - includes '\0' (size >= len+1)
};

//
// Creation, reset, free and memory expansion
//

// Constructors
STRING_BUFFER* string_buff_init(const t_buf_pos size);
STRING_BUFFER* string_buff_create(const char* str);

// Destructors
void string_buff_free(STRING_BUFFER* sbuf);
// Free sbuf struct, but retain and return the char array
char* string_buff_free_get_str(STRING_BUFFER* sbuf);

// Clone a buffer (including content)
STRING_BUFFER* string_buff_clone(const STRING_BUFFER* sbuf);

// Get a copy of this STRING_BUFFER as a char array
// Returns NULL if not enough memory
char* string_buff_as_str(const STRING_BUFFER* sbuf);

// Clear the content of an existing STRING_BUFFER (sets size to 0)
void string_buff_reset(STRING_BUFFER* sbuf);

// Get size
inline t_buf_pos string_buff_strlen(const STRING_BUFFER* sbuf);
inline t_buf_pos string_buff_size(const STRING_BUFFER* sbuf);

//
// Resizing
//

// Ensure capacity for len characters plus '\0' character - exits on FAILURE
void string_buff_ensure_capacity(STRING_BUFFER *sbuf, const t_buf_pos len);

// reallocs to exact memory specified - return 1 on success 0 on failure
char string_buff_resize(STRING_BUFFER *sbuf, const t_buf_pos new_size);

// convenience function: prints error and exits with EXIT_FAILURE if it fails
void string_buff_resize_vital(STRING_BUFFER *sbuf, const t_buf_pos new_size);

// Shorten string without reallocating memory
void string_buff_shrink(STRING_BUFFER *sbuf, const t_buf_pos new_len);

//
// Useful String functions
//

// get/set chars
inline char string_buff_get_char(const STRING_BUFFER *sbuf,
                                 const t_buf_pos index);
inline void string_buff_set_char(STRING_BUFFER *sbuf, const t_buf_pos index,
                                 const char c);

// Add a character to the end of this STRING_BUFFER
void string_buff_append_char(STRING_BUFFER* sbuf, const char txt);
// Copy a STRING_BUFFER to the end of this STRING_BUFFER
void string_buff_append_buff(STRING_BUFFER* dst, STRING_BUFFER* src);
// Copy a character array to the end of this STRING_BUFFER
void string_buff_append_str(STRING_BUFFER* sbuf, const char* txt);
// Copy N characters from a character array to the end of this STRING_BUFFER
void string_buff_append_strn(STRING_BUFFER* sbuf, const char* txt,
                             const t_buf_pos len);

// Remove \r and \n characters from the end of this STRING_BUFFER
void string_buff_chomp(STRING_BUFFER *sbuf);

// Get a substring as a new null terminated char array
// (remember to free the returned char* after you're done with it!)
char* string_buff_substr(STRING_BUFFER *sbuf, const t_buf_pos start,
                         const t_buf_pos len);

// Change to upper or lower case
void string_buff_to_uppercase(STRING_BUFFER *sbuf);
void string_buff_to_lowercase(STRING_BUFFER *sbuf);

// Copy a string to this STRING_BUFFER, overwriting any existing characters
void string_buff_copy(STRING_BUFFER* dest, const t_buf_pos dest_pos,
                      const STRING_BUFFER* src, const t_buf_pos src_pos,
                      const t_buf_pos len);

void string_buff_str_copy(STRING_BUFFER* dst, const t_buf_pos dst_pos,
                          const char* src, const t_buf_pos len);

// Copy to a STRING_BUFFER, shifting any existing characters along
void string_buff_insert(STRING_BUFFER* dest, const t_buf_pos dest_pos,
                        const STRING_BUFFER* src, const t_buf_pos src_pos,
                        const t_buf_pos len);

void string_buff_str_insert(STRING_BUFFER* dst, const t_buf_pos dst_pos,
                            const char* src, const t_buf_pos len);

//
// sprintf
//

// sprintf to a STRING_BUFFER (adds string terminator after sprint)
void string_buff_sprintf(STRING_BUFFER *sbuf, const char* fmt, ...);
void string_buff_sprintf_at(STRING_BUFFER *sbuf, const t_buf_pos pos,
                            const char* fmt, ...);
void string_buff_vsprintf(STRING_BUFFER *sbuf, const t_buf_pos pos,
                          const char* fmt, va_list argptr);

// sprintf without terminating character
// Does not prematurely end the string if you sprintf within the string
// (terminates string if sprintf to the end)
void string_buff_sprintf_noterm(STRING_BUFFER *sbuf, const t_buf_pos pos,
                                const char* fmt, ...);

//
// Reading files
//

// Reading a FILE
t_buf_pos string_buff_reset_readline(STRING_BUFFER *sbuf, FILE *file);
t_buf_pos string_buff_readline(STRING_BUFFER *sbuf, FILE *gz_file);

// Reading a gzFile
t_buf_pos string_buff_reset_gzreadline(STRING_BUFFER *sbuf, gzFile gz_file);
t_buf_pos string_buff_gzreadline(STRING_BUFFER *sbuf, gzFile gz_file);

// Skip a line and return how many characters were skipped
t_buf_pos string_buff_skip_line(FILE *file);
t_buf_pos string_buff_gzskip_line(gzFile gz_file);


/**************************/
/* Other String functions */
/**************************/

char string_is_all_whitespace(const char* s);
char* string_next_nonwhitespace(const char* s);
char* string_trim(char* str);
size_t string_chomp(char* str);
size_t string_count_char(const char* str, const int c);
long string_split(const char* split, const char* txt, char*** result);

#endif
