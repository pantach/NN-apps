#ifndef HT_H
#define HT_H

#include "math.h"
#include "vector.h"

typedef struct Hashtable Hashtable;

// Hashtable errors
typedef enum {
	HT_SUCCESS = 0,  // No errors
	HT_ENOMEM,       // Memory failure
	HT_ENOHASH       // No hashing function provided
} Hashtable_errcode;

// Hashtable modes
typedef enum {
	HT_MODE_LSH = 0, // LSH mode
	HT_MODE_CUBE     // Hypecube mode
} Hashtable_mode;

// Hashtable set option operations
typedef enum {
	HT_HASHFUNC = 0, // Hashing function
	HT_METRIC,       // Comparator function
	HT_FREEFUNC,     // Freeing function
	HT_HASHDATA,     // Callback data for the hashing function
	HT_CUBE_DIM,     // Dimension of the hypercube
	HT_CUBE_PROBES,  // Max hamming distance
	HT_CUBE_M        // Max points to check
} Hashtable_op;

typedef long   (*Hashtable_hashfunc)(const Math_vector*, void*);
typedef double (*Hashtable_metric)(const Math_vector*, const Math_vector*);
typedef void   (*Hashtable_freefunc)(const Math_vector*);

char* hashtable_error(Hashtable_errcode);

Hashtable* hashtable_init(Hashtable_mode mode, size_t size);
void hashtable_free(Hashtable* ht);
void hashtable_setopt(Hashtable* ht, Hashtable_op op, ...);
int  hashtable_insert(Hashtable* ht, Math_vector* v);
size_t hashtable_nentries(Hashtable* ht);

// LSH search
Vector* hashtable_search(Hashtable* ht, const Math_vector* q);

// Range search
Vector* hashtable_rsearch(Hashtable* ht, const Math_vector* q, float radius);

// Exhaustive search
Math_vector* hashtable_esearch(Hashtable* ht, const Math_vector* q);

#endif
