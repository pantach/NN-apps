#ifndef BITMAP_H
#define BITMAP_H

typedef struct Bitmap Bitmap;
typedef char Bit;

// Bitmap errors
typedef enum {
	BITMAP_SUCCESS = 0,   // No errors
	BITMAP_ENOMEM,        // Memory failure
	BITMAP_ENOKEY,        // No such key
	BITMAP_EDUPKEY        // Duplicate key
} Bitmap_errcode;

Bitmap*  bitmap_init(size_t size);
void  bitmap_free(Bitmap* map);
int   bitmap_insert(Bitmap* map, int key, Bit bit);
Bit   bitmap_get(Bitmap* map, int key);
char* bitmap_error(Bitmap_errcode error);
size_t bitmap_n(Bitmap* map);

#endif
