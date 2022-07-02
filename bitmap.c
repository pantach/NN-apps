#include <stdlib.h>
#include <stdarg.h>
#include "bitmap.h"
#include "math.h"

typedef struct Bucket Bucket;

struct Bucket {
	int key;
	Bit bit;
	Bucket* next;
};

struct Bitmap {
	Bucket** table; // Hash table index
	size_t size;    // Index size
	size_t n;       // Number of currently stored entries
};

static size_t  bitmap_hash(Bitmap* map, int key);

static Bucket* bucket_init(void);
static void    bucket_free(Bucket* b);
static int     bucket_insert(Bucket** b, int key, Bit bit);
static Bit     bucket_find(Bucket* b, int key);


/* Implementation */

char* bitmap_error(Bitmap_errcode errcode)
{
	struct Bitmap_err {
		Bitmap_errcode errcode;
		char* errmsg;
	} err[] = {
		{ BITMAP_SUCCESS, "Success" },
		{ BITMAP_ENOMEM,  "Out of memory" },
		{ BITMAP_ENOKEY,  "The specified key does not exist" },
		{ BITMAP_EDUPKEY, "The specified key already exists" }
	};

	return err[errcode].errmsg;
}

static size_t bitmap_hash(Bitmap* map, int key)
{
	return mod(key, map->size);
}

Bitmap* bitmap_init(size_t size)
{
	Bitmap* map = malloc(sizeof(*map));

	if (map) {
		map->table = calloc(size, sizeof(*map->table));
		if (!map->table) {
			free(map);
			return NULL;
		}

		map->size = size;
		map->n = 0;
	}

	return map;
}

void bitmap_free(Bitmap* map)
{
	for (int i = 0; i < map->size; i++)
		bucket_free(map->table[i]);

	free(map->table);
	free(map);
}

size_t bitmap_n(Bitmap* map)
{
	return map->n;
}

int bitmap_insert(Bitmap* map, int key, Bit bit)
{
	int hash;
	int ret;

	hash = bitmap_hash(map, key);

	ret = bucket_insert(&map->table[hash], key, bit);
	if (ret == BITMAP_SUCCESS)
		map->n++;

	return ret;
}

Bit bitmap_get(Bitmap* map, int key)
{
	int hash;

	hash = bitmap_hash(map, key);

	if (map->table[hash] == NULL)
		return BITMAP_ENOKEY;

	return bucket_find(map->table[hash], key);
}

static Bucket* bucket_init(void)
{
	Bucket* b = calloc(1, sizeof(*b));

	return b;
}

static void bucket_free(Bucket* b)
{
	Bucket* temp;

	if (b) {
		do {
			temp = b->next;
			free(b);
		} while ((b = temp));
	}
}

static int bucket_insert(Bucket** b, int key, Bit bit)
{
	while (*b) {
		if ((*b)->key == key)
			return BITMAP_EDUPKEY;

		b = &(*b)->next;
	}

	*b = bucket_init();
	if (!*b)
		return BITMAP_ENOMEM;

	(*b)->key = key;
	(*b)->bit = bit;

	return BITMAP_SUCCESS;
}

static Bit bucket_find(Bucket* b, int key)
{
	while (b) {
		if (b->key == key)
			return b->bit;

		b = b->next;
	}

	return BITMAP_ENOKEY;
}
