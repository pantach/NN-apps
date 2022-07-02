#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "cube.h"
#include "tools.h"
#include "io.h"
#include "hashtable.h"
#include "bitmap.h"

#define MAX_PATH 512
#define DATASET_PCT 0.01   // The percentage of the dataset to use for the
					  	   // calculation of w

/* A struct holding all the necessary variables for the hashing function */
struct Hashdata {
	Bitmap** map;    // An array of k maps keeping track of f's mappings
	int k;           // Number of f/h functions
	int w;           // Window size
	Math_vector** v; // Array of k random v vectors
	float* t;        // Array of k random t's
};

struct Cube {
	Hashtable* ht;
	struct Hashdata data;
	int M;
	double knn_time;
	double ex_time;
};

static struct Hashdata hashdata_init(Vector* dataset, int k, int w);
static void hashdata_free(struct Hashdata data);
static int calc_w(Vector* dataset);
static Math_vector* create_v(size_t dim);

static long F(Math_vector* p, void* hashdata);
static Bit f(Math_vector* p, struct Hashdata* data, int i);
static int h(Math_vector* p, Math_vector* v, float t, int w);

static int sort_by_euclid_dist(const void* v, const void* u, void* q);
static int sort_by_id(const void* pv, const void* pu, void* unused);
static int id_comp(const void* v, const void* u, void* unsused);


/* Calculate the hashing variables
 * Create and fill the hashtable
 * */
Cube* cube_init(Vector* dataset, size_t ht_size, size_t k)
{
	Cube* cube;
	int w;

	cube = xmalloc(sizeof(*cube));

	w = calc_w(dataset);

	cube->ht   = hashtable_init(HT_MODE_CUBE, ht_size);
	cube->data = hashdata_init(dataset, k, w);

	hashtable_setopt(cube->ht, HT_HASHFUNC, F);
	hashtable_setopt(cube->ht, HT_HASHDATA, &cube->data);
	hashtable_setopt(cube->ht, HT_CUBE_DIM, k);

	for (int j = 0; j < dataset->size; ++j) {
		int err = hashtable_insert(cube->ht, vector_get(dataset, j));
		if (err)
			err_exit("%s\n", hashtable_error(err));
	}

	return cube;
}

void cube_free(Cube* cube)
{
	if (cube) {
		hashtable_free(cube->ht);
		hashdata_free(cube->data);
		free(cube);
	}
}

void cube_setopt(Cube* cube, Cube_opt opt, ...)
{
	Cube_metric metric;
	int probes;
	int M;
	va_list args;

	va_start(args, opt);
	switch (opt) {
	case CUBE_METRIC:
		metric = va_arg(args, Cube_metric);
		hashtable_setopt(cube->ht, HT_METRIC, metric);
		break;
	case CUBE_PROBES:
		probes = va_arg(args, int);
		hashtable_setopt(cube->ht, HT_CUBE_PROBES, probes);
		break;
	case CUBE_M:
		M = va_arg(args, int);
		hashtable_setopt(cube->ht, HT_CUBE_M, M);
		cube->M = M;
		break;
	}
	va_end(args);
}

double cube_kNN_time(Cube* cube)
{
	return cube->knn_time;
}

double cube_exhaust_time(Cube* cube)
{
	return cube->ex_time;
}

Math_vector* cube_NN(Cube* cube, Math_vector* q)
{
	Vector* results;
	Math_vector* NN = NULL;
	Math_vector* v;
	double min_dist = DBL_MAX;
	double dist;
	clock_t start, end;

	start = clock();
	results = hashtable_search(cube->ht, q);

	vector_sort (results, sort_by_id, NULL);
	vector_rmdup(results, id_comp, NULL);

	for (int i = 0; i < results->size; i++) {
		v = vector_get(results, i);
		dist = euclidean_distance(q, v);
		if (dist < min_dist) {
			min_dist = dist;
			NN = v;
		}
	}
	end = clock();
	cube->knn_time = (double) (end -start)/CLOCKS_PER_SEC;

	vector_free(results);

	return NN;
}

Vector* cube_kNN(Cube* cube, Math_vector* q, int N)
{
	clock_t start, end;
	Vector* results;
	int max;

	start = clock();
	results = hashtable_search(cube->ht, q);
	end = clock();

	max = (N < cube->M) ? N : cube->M;

	vector_sort (results, sort_by_id, NULL);
	vector_rmdup(results, id_comp, NULL);
	vector_sort (results, sort_by_euclid_dist, q);
	vector_trunc(results, max);

	cube->knn_time = (double) (end -start)/CLOCKS_PER_SEC;

	return results;
}

Vector* cube_range_search(Cube* cube, Math_vector* q, float radius)
{
	Vector* results = hashtable_rsearch(cube->ht, q, radius);

	vector_sort (results, sort_by_id, NULL);
	vector_trunc(results, cube->M);

	return results;
}

Math_vector* cube_exhaust_search(Cube* cube, Math_vector* q)
{
	Math_vector* result;
	clock_t start, end;

	start = clock();
	result = hashtable_esearch(cube->ht, q);
	end = clock();

	cube->ex_time = (double) (end -start)/CLOCKS_PER_SEC;

	return result;
}

static int sort_by_id(const void* pv, const void* pu, void* unused)
{
	Math_vector* v = *(Math_vector**)pv;
	Math_vector* u = *(Math_vector**)pu;

	return math_vector_id(v) -math_vector_id(u);
}


static int id_comp(const void* v, const void* u, void* unsused)
{
	return math_vector_compid(v, u);
}

static int sort_by_euclid_dist(const void* pv, const void* pu, void* pq)
{
	Math_vector* v = *(Math_vector**)pv;
	Math_vector* u = *(Math_vector**)pu;
	Math_vector* q =  (Math_vector*) pq;
	double vq_dist = euclidean_distance(v, q);
	double uq_dist = euclidean_distance(u, q);

	if (vq_dist < uq_dist)
		return -1;
	else if (vq_dist > uq_dist)
		return 1;
	else
		return 0;
}

/* The function F concatenates all the bit strings associated with each f_i into
 * a k-length bit string, which acts as a hashing function for the hhypecube */
static long F(Math_vector* p, void* hashdata)
{
	struct Hashdata* data = (struct Hashdata*)hashdata;
	long res = 0;

	for (int i = 0; i < data->k; ++i)
		res |= f(p, data, i) << i;

	return res;
}

/* The function f_i(h_i) maps each output of h_i to a random bit value, stores
 * it and then returns it */
static Bit f(Math_vector* p, struct Hashdata* data, int i)
{
	int hval;
	Bit bit;
	int err;

	hval = h(p, data->v[i], data->t[i], data->w);
	bit  = bitmap_get(data->map[i], hval);
	if (bit == BITMAP_ENOKEY) {
		bit = randint(0, 1);
		err = bitmap_insert(data->map[i], hval, bit);
		if (err)
			err_exit("Cannot insert value into bitmap: %s", bitmap_error(err));
	}

	return bit;
}

/* The h_i hashing function */
static int h(Math_vector* p, Math_vector* v, float t, int w)
{
	return floor((dot(p, v) + t)/(float)w);
}

static struct Hashdata hashdata_init(Vector* dataset, int k, int w)
{
	struct Hashdata data;
	Math_vector* p = vector_get(dataset, 0);

	data.map = xmalloc(k * sizeof(*data.map));
	data.v   = xmalloc(k * sizeof(*data.v));
	data.t   = xmalloc(k * sizeof(*data.t));

	data.k = k;
	data.w = w;

	for (int i = 0; i < data.k; ++i) {
		data.map[i] = bitmap_init(1000);
		data.t[i]   = randfloat(0, data.w);
		data.v[i]   = create_v(math_vector_dim(p));
	}

	return data;
}

static void hashdata_free(struct Hashdata data)
{
	for (int i = 0; i < data.k; ++i) {
		math_vector_free(data.v[i]);
		bitmap_free(data.map[i]);
	}

	free(data.v);
	free(data.t);
	free(data.map);
}

static int double_sort(const void* p1, const void* p2)
{
	const double* d1 = p1;
	const double* d2 = p2;

	if (*d1 > *d2)
		return 1;

	if (*d1 < *d2)
		return -1;

	return 0;
}

static int calc_subset_size(Vector* dataset)
{
	if (dataset->size <= 10)
		return dataset->size;
	if (dataset->size <= 100)
		return 0.5 * dataset->size;
	if (dataset->size <= 1000)
		return 0.1 * dataset->size;
	else
		return 0.05 * dataset->size;
}

/* Calculates w by calculating the median distance between pairs of different
 * random vectors, of a sample subset of the dataset */
static int calc_w(Vector* dataset)
{
	int subset_size = calc_subset_size(dataset);
	double distance[subset_size];
	Math_vector* v;
	Math_vector* u;
	int r1;
	int r2;
	int median;

	if (dataset->size < 2)
		return 1;

	// Pick two different random vectors, calculate their distance and add it
	// to the array
	for (int i = 0; i < subset_size; i++) {
		v = vector_get(dataset, r1 = randint(0, dataset->size -1));
		do {
			u = vector_get(dataset, r2 = randint(0, dataset->size -1));
		} while (r1 == r2);

		distance[i] = euclidean_distance(v, u);
	}

	qsort(distance, subset_size, sizeof(double), double_sort);

	// Return the median distance divided by 10
	median = distance[subset_size/2]/10;
	if (median == 0)
		return 4;
	else
		return median;
}

/* Creates the d-vector v ~ N(0,1)^d */
static Math_vector* create_v(size_t dim)
{
	Math_vector* v = math_vector_init(VECTOR_MODE, 0);

	for (int i = 1; i < dim; ++i)
		math_vector_push(v, uniform_rand(0, 1));

	return v;
}
