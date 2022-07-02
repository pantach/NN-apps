#ifndef PROGBAR_H
#define PROGBAR_H

typedef struct Progbar Progbar;

Progbar* progress_bar_init(size_t size, unsigned int num);
void progress_bar_free(Progbar* b);
void progress_bar_print(Progbar* b, unsigned int i);

#endif
