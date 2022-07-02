#ifndef TOOLS_H
#define TOOLS_H

#include <stdarg.h>
#include <stddef.h>

void* xmalloc(size_t size);
void* xcalloc(size_t num, size_t size);
void* xrealloc(void* ptr, size_t size);
void* memdup(const void* src, size_t n);
char* xstrdup(const char* str);
void  xstrcat(char** str1, const char* str2);

void printff(const char* format, ...);
void err_exit(const char* format, ...);
void syserr_exit(const char* format, ...);

int xsprintf(char** str, const char* format, ...);

void print_bits(size_t const size, void const * const ptr);

#define GETNUM_NOFLAGS 0
#define GETNUM_NOEXIT  1

typedef enum {
	GETNUM_ESUCCESS =  0,
	GETNUM_ENULLSTR = -1,
	GETNUM_EEMPTSTR = -2,
	GETNUM_ENONNUM  = -3,
	GETNUM_ESTRTOL  = -4,
	GETNUM_ESTRTOD  = -5
} Getnum_err;

int   getint  (const char* numstr, int flags, ...);
float getfloat(const char* numstr, int flags, ...);

int binary_search(const void* base, const void* key, size_t nmemb, size_t size,
                  int comp(const void*, const void*, void*), void* comp_data);
#endif
