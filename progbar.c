#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "progbar.h"

struct Progbar {
	char* str;
	size_t size;
	unsigned int num;
	unsigned int spindex;
};

const char g_spinner[] = "|/-\\";

Progbar* progress_bar_init(size_t size, unsigned int num)
{
	Progbar* b = malloc(sizeof(*b));
	if (!b)
		return NULL;

	b->str = malloc(size +2);
	if (!b->str) {
		free(b);
		return NULL;
	}

	memset(b->str, '=', size);
	b->size = size;
	b->num = num;
	b->spindex = 0;

	return b;
}

void progress_bar_free(Progbar* b)
{
	free(b->str);
	free(b);
}

void progress_bar_print(Progbar* b, unsigned int i)
{
	unsigned int progress = (int)((i/(float)b->num)*b->size);

	b->str[progress] = '>';
	b->str[progress +1] = '\0';

	printf("\r[ %-*s ] %c %.1f%%", (int)b->size, b->str,
	       g_spinner[b->spindex % 4], (i/(float)b->num)*100) ;
	fflush(stdout);

	b->spindex++;

	b->str[progress] = '=';
	b->str[progress +1] = '=';
}
