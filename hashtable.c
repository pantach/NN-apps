#include <stdlib.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>
#include "hashtable.h"

#define EPSILON 0.1

typedef struct Bucket Bucket;

typedef struct {
	long ID;            // ID of the bucket
	Math_vector* v;     // Vector stored
} Bucket_data;

struct Bucket {
	Bucket_data data;
	Bucket* next;   // Next Bucket
};

struct Cube_data {
	int dim;
	int probes;
	int max_points;
};

struct Hashtable {
	Hashtable_mode mode;  // LSH or CUBE
	Bucket** table;       // Hash table index
	size_t size;          // Index size
	size_t n;             // Number of currently stored entries
	float radius;         // Radius in range search
	void* hashdata;       // Callback data for the hashing function

	// Hypercube-specific data
	struct Cube_data cube;

	/* A pointer to the hashing function provided by the user. */
	Hashtable_hashfunc hashfunc;

	/* A pointer to the comparator function provided by the user. This function
	   should return 0 when key1 and key2 are the same */
	Hashtable_metric metric;

	/* A pointer to the freeing function provided by the user. */
	Hashtable_freefunc freefunc;
};

static Vector* hashtable_lsh_search(Hashtable* ht, const Math_vector* q);
static Vector* hashtable_cube_search(Hashtable* ht, const Math_vector* q);
static Vector* hashtable_search_(Hashtable* ht, const Math_vector* q);

static Bucket* bucket_init(void);
static void    bucket_free(Bucket* b, Hashtable_freefunc freefunc);
static int     bucket_insert(Bucket** b, long ID, Math_vector* v);
static Vector* bucket_search(Bucket* b, long ID);
static Vector* bucket_rsearch(Bucket* b, long ID, const Math_vector* q,
                              float radius, Hashtable_metric metric);
static Math_vector* bucket_nearest_search(Bucket* b, const Math_vector* q,
                                          Hashtable_metric metric);

static inline double relative_error(long ID, long qID)
{
	if (ID == 0 && (ID -qID) == 0)
		return 0;

	if (ID == 0)
		return INFINITY;

	return (double) (ID -qID)/ID;
}

/* Hypercube */
static long invert_bits(long bitstr, const Math_vector* bitpos);


/* Implementation */

static long invert_bits(long bitstr, const Math_vector* bitpos)
{
	for (int i = 0; i < math_vector_dim(bitpos); ++i)
		bitstr ^= (1 << (int) math_vector_val(bitpos, i));

	return bitstr;
}

char* hashtable_error(Hashtable_errcode errcode)
{
	struct Hashtable_err {
		Hashtable_errcode errcode;
		char* errmsg;
	} err[] = {
		{ HT_SUCCESS, "Success" },
		{ HT_ENOMEM,  "Out of memory" },
		{ HT_ENOHASH, "No hashing function provided" }
	};

	return err[errcode].errmsg;
}

Hashtable* hashtable_init(Hashtable_mode mode, size_t size)
{
	Hashtable* ht = calloc(1, sizeof(*ht));

	if (ht) {
		ht->table = calloc(size, sizeof(*ht->table));
		if (!ht->table) {
			free(ht);
			return NULL;
		}

		ht->mode  = mode;
		ht->size  = size;
	}

	return ht;
}

void hashtable_free(Hashtable* ht)
{
	for (int i = 0; i < ht->size; i++)
		bucket_free(ht->table[i], ht->freefunc);

	free(ht->table);
	free(ht);
}

void hashtable_setopt(Hashtable* ht, Hashtable_op op, ...)
{
	va_list args;

	va_start(args, op);
	switch (op) {
	case HT_HASHFUNC:
		ht->hashfunc = va_arg(args, Hashtable_hashfunc);
		break;
	case HT_METRIC:
		ht->metric = va_arg(args, Hashtable_metric);
		break;
	case HT_FREEFUNC:
		ht->freefunc = va_arg(args, Hashtable_freefunc);
		break;
	case HT_HASHDATA:
		ht->hashdata = va_arg(args, void*);
		break;
	case HT_CUBE_DIM:
		ht->cube.dim = va_arg(args, int);
		break;
	case HT_CUBE_PROBES:
		ht->cube.probes = va_arg(args, int);
		break;
	case HT_CUBE_M:
		ht->cube.max_points = va_arg(args, int);
		break;
	}

	va_end(args);
}

size_t hashtable_nentries(Hashtable* ht)
{
	return ht->n;
}

int hashtable_insert(Hashtable* ht, Math_vector* v)
{
	long ID;
	int hash;
	int ret;

	if (!ht->hashfunc)
		return HT_ENOHASH;

	ID = ht->hashfunc(v, ht->hashdata);
	hash = mod(ID, ht->size);

	ret = bucket_insert(&ht->table[hash], ID, v);
	if (ret == HT_SUCCESS)
		ht->n++;

	return ret;
}

Vector* hashtable_search(Hashtable* ht, const Math_vector* q)
{
	ht->radius = -1;

	return hashtable_search_(ht, q);
}

Vector* hashtable_rsearch(Hashtable* ht, const Math_vector* q, float radius)
{
	ht->radius = radius;

	return hashtable_search_(ht, q);
}

static Vector* hashtable_search_(Hashtable* ht, const Math_vector* q)
{
	switch (ht->mode) {
	case HT_MODE_LSH:
		return hashtable_lsh_search(ht, q);
	case HT_MODE_CUBE:
		return hashtable_cube_search(ht, q);
	}

	return NULL;
}

static Vector* hashtable_lsh_search(Hashtable* ht, const Math_vector* q)
{
	Vector* results_aggr = vector_init(NULL);
	Vector* results;
	long ID, i;
	int hash;

	if (!ht->hashfunc)
		return NULL;

	ID = ht->hashfunc(q, ht->hashdata);

	for (i = floor(ID -abs(ID*EPSILON)); i <= ceil(ID +abs(ID*EPSILON)); i++) {
		hash = mod(i, ht->size);

		if (ht->radius == -1)
			results = bucket_search(ht->table[hash], i);
		else
			results = bucket_rsearch(ht->table[hash], i, q, ht->radius,
									 ht->metric);

		vector_concat(results_aggr, results);
		vector_free(results);
	}

	return results_aggr;
}

static Vector* hashtable_cube_search(Hashtable* ht, const Math_vector* q)
{
	Vector* results_aggr = vector_init(NULL);
	Vector* results;
	Vector* bitcombos;
	Math_vector* bitpos;
	long ID;
	long altID;
	int hash;

	if (!ht->hashfunc)
		return NULL;

	ID = ht->hashfunc(q, ht->hashdata);

	for (int hamdist = 0; hamdist <= ht->cube.probes; ++hamdist) {
		bitcombos = k_combinations(ht->cube.dim, hamdist);

		while ((bitpos = vector_pop(bitcombos))) {

			altID = invert_bits(ID, bitpos);
			hash = mod(altID, ht->size);

			math_vector_free(bitpos);

			if (ht->table[hash]) {
				if (ht->radius == -1)
					results = bucket_search(ht->table[hash], altID);
				else
					results = bucket_rsearch(ht->table[hash], altID, q,
					                         ht->radius, ht->metric);

				vector_concat(results_aggr, results);
				vector_free(results);
			}

			if (results_aggr->size >= ht->cube.max_points) {
				vector_free(bitcombos);
				goto finish;
			}

		}
		vector_free(bitcombos);
	}

finish:
	return results_aggr;
}

Math_vector* hashtable_esearch(Hashtable* ht, const Math_vector* q)
{
	Math_vector* nearest = NULL;
	Math_vector* near;
	double nearest_dist  = DBL_MAX;
	double dist;

	for (int i = 0; i < ht->size; ++i) {
		near = bucket_nearest_search(ht->table[i], q, ht->metric);
		if (near) {
			if ((dist = ht->metric(q, near)) < nearest_dist) {
				nearest = near;
				nearest_dist = dist;
			}
		}
	}

	return nearest;
}

static Bucket* bucket_init(void)
{
	Bucket* b = calloc(1, sizeof(*b));

	return b;
}

static void bucket_free(Bucket* b, Hashtable_freefunc freefunc)
{
	Bucket* btemp;

	if (b) {
		do {
			if (freefunc)
				freefunc(b->data.v);

			btemp = b->next;
			free(b);
		} while ((b = btemp));
	}
}

static int bucket_insert(Bucket** b, long ID, Math_vector* v)
{
	while (*b != NULL)
		b = &(*b)->next;

	*b = bucket_init();
	if (!*b)
		return HT_ENOMEM;

	(*b)->data.ID = ID;
	(*b)->data.v  = v;

	return HT_SUCCESS;
}

static Vector* bucket_search(Bucket* b, long ID)
{
	Vector* results = vector_init(NULL);

	while (b != NULL) {
		if (b->data.ID == ID)
			vector_push(results, b->data.v);

		b = b->next;
	}

	return results;
}

static Vector* bucket_rsearch(Bucket* b, long ID, const Math_vector* q,
                              float radius, Hashtable_metric metric)
{
	Vector* results = vector_init(NULL);

	while (b != NULL) {
		if (b->data.ID == ID)
			if (metric(b->data.v, q) <= radius)
				vector_push(results, b->data.v);

		b = b->next;
	}

	return results;
}

static Math_vector* bucket_nearest_search(Bucket* b, const Math_vector* q,
                                          Hashtable_metric metric)
{
	Math_vector* nearest = NULL;
	double nearest_dist  = DBL_MAX;
	double dist;

	while (b != NULL) {
		if ((dist = metric(b->data.v, q)) < nearest_dist) {
			nearest = b->data.v;
			nearest_dist = dist;
		}

		b = b->next;
	}

	return nearest;
}
