#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "tools.h"

typedef enum {
	GETNUM_INT = 0,
	GETNUM_FLOAT
} Number_type;

typedef union {
	long i;
	double f;
} Number;

static
Number getnum(const char* numstr, Number_type ntype, int flags, va_list args);

static
int binary_search_(const void* base, const void* key, size_t nmemb,
                   int comp(const void*, const void*, void*), void* comp_data,
                   int left, int right);

void* xmalloc(size_t size)
{
	void* ptr = malloc(size);
	if (!ptr)
		abort();

	return ptr;
};

void* xcalloc(size_t num, size_t size)
{
	void* ptr = calloc(num, size);
	if (!ptr)
		abort();

	return ptr;
}

void* xrealloc(void* ptr, size_t size)
{
	void* newptr = realloc(ptr, size);
	if (!newptr)
		abort();

	return newptr;
}


void* memdup(const void* src, size_t n)
{
	char* dest = xmalloc(n);

	return memcpy(dest, src, n);
}

char* xstrdup(const char* str)
{
	char*  str2 = NULL;
	size_t len;

	if (str) {
		len  = strlen(str);
		str2 = xmalloc(len +1);
		memcpy(str2, str, len +1);
	}

	return str2;
}

void xstrcat(char** str1, const char* str2)
{
	size_t len1;
	size_t len2;
	char*  cat;

	len1 = (*str1 == NULL) ? 0 : strlen(*str1);
	len2 = strlen(str2);

	cat = xmalloc(len1 +len2 +1);

	if (*str1)
		memcpy(cat, *str1, len1);
	memcpy(cat +len1, str2, len2);
	cat[len1 +len2] = '\0';

	free(*str1);
	*str1 = cat;
}

void printff(const char* format, ...)
{
	va_list args;

	va_start(args, format);
	vprintf(format, args);
	va_end(args);
	fflush(stdout);
}

void err_exit(const char* format, ...)
{
	va_list args;
	char nformat[1024];

	strncpy(nformat, format, 1023);
	strcat (nformat, "\n");

	va_start(args, format);
	vfprintf(stderr, nformat, args);
	va_end(args);

	fflush(stdout);
	exit(EXIT_FAILURE);
}

void syserr_exit(const char* format, ...)
{
	va_list args;
	char usermsg[1024];
	int errno_copy = errno;

	va_start(args, format);
	vsnprintf(usermsg, 1024, format, args);
	fprintf(stderr, "%s: %s\n", usermsg, strerror(errno_copy));
	va_end(args);

	fflush(stdout);
	exit(EXIT_FAILURE);
}

int xsprintf(char** str, const char* format, ...)
{
	va_list args;
	int len;
	int wlen;

	va_start(args, format);
	len = vsnprintf(NULL, 0, format, args);
	va_end(args);

	*str = xmalloc(len +1);

	va_start(args, format);
	wlen = vsnprintf(*str, len +1, format, args);
	va_end(args);

	return wlen;
}

// Assumes little endian
void print_bits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;

    for (i = size -1; i >= 0; i--) {
        for (j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}

int getint(const char* numstr, int flags, ...)
{
	va_list args;
	Number num;

	va_start(args, flags);
	num = getnum(numstr, GETNUM_INT, flags, args);
	va_end(args);

	return num.i;
}

float getfloat(const char* numstr, int flags, ...)
{
	va_list args;
	Number num;

	va_start(args, flags);
	num = getnum(numstr, GETNUM_FLOAT, flags, args);
	va_end(args);

	return num.f;
}

/* Converts the integral number represented by the string "numstr" to type
 * integer.
 */
static
Number getnum(const char* numstr, Number_type ntype, int flags, va_list args)
{
	Number num;
	char* endptr;
	Getnum_err  err  = 0;
	Getnum_err* perr = NULL;

	if (flags & GETNUM_NOEXIT)
		perr = va_arg(args, Getnum_err*);

	if (numstr) {
		errno = 0;

		switch (ntype) {
			case GETNUM_FLOAT:
				num.f = strtod(numstr, &endptr);
				if (errno)
					err = GETNUM_ESTRTOD;

				break;

			default:
				num.i = strtol(numstr, &endptr, 10);
				if (errno)
					err = GETNUM_ESTRTOL;
		}

		if (!err) {
			if (endptr == numstr)
				err = GETNUM_EEMPTSTR;

			else if (*endptr != '\0')
				err = GETNUM_ENONNUM;
		}
	}
	else
		err = GETNUM_ENULLSTR;

	if (flags ^ GETNUM_NOEXIT) {
		switch (err) {
		case GETNUM_ENULLSTR:
			err_exit("getnum: Null input string");
		case GETNUM_EEMPTSTR:
			err_exit("getnum: Empty input string");
		case GETNUM_ENONNUM:
			err_exit("getnum: Non-numeric characters: %s", numstr);
		case GETNUM_ESTRTOL:
			syserr_exit("getnum: strtol failed");
		case GETNUM_ESTRTOD:
			syserr_exit("getnum: strtod failed");
		default:
			break;
		}
	}

	if (perr)
		*perr = err;

	return num;
}

/* Performs a binary search in the given array (pre-sorted). If the key
 * searched is not found, the index of the key with the next greatest value will
 * be returned */
int binary_search(const void* base, const void* key, size_t nmemb, size_t size,
                  int comp(const void*, const void*, void*), void* comp_data)
{
	return binary_search_(base, key, nmemb, comp, comp_data, 0, size -1);
}

static
int binary_search_(const void* base_, const void* key, size_t nmemb,
                   int comp(const void*, const void*, void*), void* comp_data,
                   int left, int right)
{
	int mid;
	double ret;

	const char* base = (const char*)base_;

	if (right -left < 2)
		return right;

	mid = (left +right)/2;

	ret = comp(key, base +mid*nmemb, comp_data);
	if (ret > 0)
		return binary_search_(base, key, nmemb, comp, comp_data, mid, right);

	if (ret < 0)
		return binary_search_(base, key, nmemb, comp, comp_data, left, mid);

	return mid;
}
