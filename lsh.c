#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include "lsh.h"
#include "tools.h"
#include "hashtable.h"

#define DATASET_PCT 0.01    // The percentage of the dataset to use for the
					  	  	// calculation of w
#define RMULT 1.5           // Max multiplier for random r
#define EPSILON 0.01        // Curve filtering allowance

typedef struct {
	double delta;   // Grid's delta value
	float t;        // Grid's random t shift value
} Grid;

/* A struct holding all the necessary variables for the hashing function */
struct Hashdata {
	int k;           // Number of h functions
	int w;           // Window size
	Math_vector** v; // Array of k random v vectors
	int* r;          // Array of k random r's (linear comb coefficients)
	float* t;        // Array of k random t's

	/* LSH for curves only */
	Grid grid;
};

typedef struct {
	Hashtable* table;
	struct Hashdata data;
} Hash;

struct Sort_data {
	const Math_vector* query;
	const LSH_metric metric;
};

struct LSH {
	LSH_mode mode;
	LSH_metric metric;
	Hash* hash;
	size_t size;
	double knn_time;
	double ex_time;
};

static struct Hashdata hashdata_init(Vector* dataset, int k, int w, Grid* grid);
static void hashdata_free(struct Hashdata* data);

static Grid grid_create(double delta);
static size_t grid_vector_dim(const Math_vector* p, const Grid* grid);

static Math_vector* create_v(size_t dim);
static int calc_w(Vector* dataset, LSH_metric metric);
static int calc_subset_size(Vector* dataset);

static long g(const Math_vector* p, void* hashdata);
static int  h(const Math_vector* p, const Math_vector* v, float t, int w);

static long hash_curve(const Math_vector* p, void* hashdata);
static Math_vector* curve_snap_to_grid(const Math_vector* curve,
                                       const Grid* grid);

static inline struct Sort_data sort_data_init(const LSH* lsh,
                                              const Math_vector* q);
static int sort_by_metric(const void* v, const void* u, void* sort_data);
static int sort_by_id(const void* pv, const void* pu, void* unused);
static int id_comp(const void* v, const void* u, void* unsused);
static int double_sort(const void* p1, const void* p2);


/* Calculate the hashing variables
 * Create and fill the hashtables
 * */
LSH* LSH_init(LSH_mode mode, LSH_metric metric, Vector* dataset, size_t ht_size,
              size_t L, size_t k, ...)
{
	LSH* lsh;
	int w;
	double delta;
	va_list curve_arg;

	lsh = xmalloc(sizeof(*lsh));
	lsh->hash = xmalloc(L * sizeof(*lsh->hash));

	lsh->mode = mode;
	lsh->metric = metric;
	lsh->size = L;

	w = calc_w(dataset, metric);

	if (mode == LSH_MODE_CURVE) {
		va_start(curve_arg, k);
		if (mode == LSH_MODE_CURVE)
			delta = va_arg(curve_arg, double);
		va_end(curve_arg);
	}

	for (int i = 0; i < L; ++i) {
		Hash* const hash = &lsh->hash[i];
		Hashtable_hashfunc hashfunc;

		hash->table = hashtable_init(HT_MODE_LSH, ht_size);

		if (mode == LSH_MODE_CURVE) {
			Grid grid = grid_create(delta);
			hash->data = hashdata_init(dataset, k, w, &grid);
			hashfunc = hash_curve;
		}
		else {
			hash->data = hashdata_init(dataset, k, w, NULL);
			hashfunc = g;
		}

		hashtable_setopt(hash->table, HT_HASHFUNC, hashfunc);
		hashtable_setopt(hash->table, HT_HASHDATA, &hash->data);
		hashtable_setopt(hash->table, HT_METRIC,   lsh->metric);

		for (int j = 0; j < dataset->size; ++j) {
			Math_vector* v = vector_get(dataset, j);

			int err = hashtable_insert(hash->table, v);
			if (err)
				err_exit("%s\n", hashtable_error(err));
		}
	}

	return lsh;
}

void LSH_free(LSH* lsh)
{
	if (lsh) {
		for (int i = 0; i < lsh->size; ++i) {
			hashtable_free(lsh->hash[i].table);
			hashdata_free(&lsh->hash[i].data);
		}

		free(lsh->hash);
		free(lsh);
	}
}

double LSH_kNN_time(const LSH* lsh)
{
	return lsh->knn_time;
}

double LSH_exhaust_time(const LSH* lsh)
{
	return lsh->ex_time;
}

Math_vector* LSH_NN(LSH* lsh, const Math_vector* q)
{
	Vector* results_aggr = vector_init(NULL);
	Vector* results;
	Math_vector* NN = NULL;
	Math_vector* v;
	double min_dist = DBL_MAX;
	double dist;
	clock_t start, end;
	int i;

	start = clock();
	for (i = 0; i < lsh->size; ++i) {
		results = hashtable_search(lsh->hash[i].table, q);
		vector_concat(results_aggr, results);
		vector_free(results);
	}

	vector_sort (results_aggr, sort_by_id, NULL);
	vector_rmdup(results_aggr, id_comp, NULL);

	for (i = 0; i < results_aggr->size; i++) {
		v = vector_get(results_aggr, i);
		dist = lsh->metric(q, v);
		if (dist < min_dist) {
			min_dist = dist;
			NN = v;
		}
	}
	end = clock();
	lsh->knn_time = (double) (end -start)/CLOCKS_PER_SEC;

	vector_free(results_aggr);

	return NN;
}

Vector* LSH_kNN(LSH* lsh, const Math_vector* q, int N)
{
	struct Sort_data sort_data = sort_data_init(lsh, q);
	Vector* results_aggr = vector_init(NULL);
	Vector* results;
	clock_t start, end;

	start = clock();
	for (int i = 0; i < lsh->size; ++i) {
		results = hashtable_search(lsh->hash[i].table, q);
		vector_concat(results_aggr, results);
		vector_free(results);
	}
	end = clock();

	lsh->knn_time = (double) (end -start)/CLOCKS_PER_SEC;

	vector_sort (results_aggr, sort_by_id, NULL);
	vector_rmdup(results_aggr, id_comp, NULL);
	vector_sort (results_aggr, sort_by_metric, &sort_data);
	vector_trunc(results_aggr, N);

	return results_aggr;
}

Vector* LSH_range_search(LSH* lsh, const Math_vector* q, float radius)
{
	Vector* results_aggr = vector_init(NULL);
	Vector* results;

	for (int i = 0; i < lsh->size; ++i) {
		results = hashtable_rsearch(lsh->hash[i].table, q, radius);
		vector_concat(results_aggr, results);
		vector_free(results);
	}
	vector_sort (results_aggr, sort_by_id, NULL);
	vector_rmdup(results_aggr, id_comp, NULL);

	return results_aggr;
}

static inline struct Sort_data
sort_data_init(const LSH* lsh, const Math_vector* q)
{
	return (struct Sort_data) { .query = q, .metric = lsh->metric };
}

static int sort_by_metric(const void* pv, const void* pu, void* sort_data)
{
	struct Sort_data* data = sort_data;
	const Math_vector* q = data->query;
	LSH_metric metric = data->metric;
	Math_vector* v = *(Math_vector**)pv;
	Math_vector* u = *(Math_vector**)pu;
	double vq_dist = metric(v, q);
	double uq_dist = metric(u, q);

	if (vq_dist < uq_dist)
		return -1;
	else if (vq_dist > uq_dist)
		return 1;
	else
		return 0;
}

static int sort_by_id(const void* pv, const void* pu, void* unused)
{
	Math_vector* v = *(Math_vector**)pv;
	Math_vector* u = *(Math_vector**)pu;

	return id_comp(v, u, NULL);
}

static int id_comp(const void* v, const void* u, void* unused)
{
	return math_vector_compid(v, u);
}

/* The g (ID) hashing function
 *
 * This function performs the summation "Σr_i*h_i" minus the last step (mod
 * Tablesize) which is performed within the hashtable
 *
 * */
static long g(const Math_vector* p, void* hashdata)
{
	struct Hashdata* data = (struct Hashdata*)hashdata;
	long sum = 0;

	for (int i = 0; i < data->k; ++i)
		sum += (data->r[i] * h(p, data->v[i], data->t[i], data->w)) % LONG_MAX;

	//printf("ID: %ld\n", sum);
	return sum;
}

/* The h_i hashing function */
static int h(const Math_vector* p, const Math_vector* v, float t, int w)
{
	//printf("t: %f | dot: %f | w: %d\n", t, dot(p, v), w);
	return floor((dot(p, v) + t)/(float)w);
}

static long hash_curve(const Math_vector* p, void* hashdata)
{
	struct Hashdata* data = (struct Hashdata*) hashdata;
	Math_vector* grid_vector;
	long ID;

	grid_vector = curve_snap_to_grid(p, &data->grid);

	ID = g(grid_vector, hashdata);

	math_vector_free(grid_vector);

	return ID;
}

 Math_vector* curve_filter(const Math_vector* curve, double epsilon)
{
	Math_vector* filtered_curve = math_vector_init(CURVE_MODE, NULL);
	Point a, b, c;
	int i;

	math_vector_push(filtered_curve, math_vector_point(curve, 0));

	for (i = 1; i < math_vector_dim(curve) -1; i++) {
		a = math_vector_point(curve, i -1);
		b = math_vector_point(curve, i);
		c = math_vector_point(curve, i +1);

		if (fabs(a.x -b.x) > epsilon && fabs(b.x -c.x) > epsilon &&
		    fabs(a.y -b.y) > epsilon && fabs(b.y -c.y) > epsilon)
		{
			math_vector_push(filtered_curve, math_vector_point(curve, i));
		}
	}

	math_vector_push(filtered_curve, math_vector_point(curve, i));

	return filtered_curve;
}

static Math_vector* curve_snap_to_grid(const Math_vector* curve,
                                       const Grid* grid)
{
	Math_vector* grid_vector = math_vector_init(VECTOR_MODE, NULL);
	const double delta = grid->delta;
	const float t = grid->t;
	Point point;
	Point grid_point;
	Point prev_point;

	for (int i = 0; i < math_vector_dim(curve); i++) {
		point = math_vector_point(curve, i);

		grid_point.x = floor((point.x -t)/delta +0.5)*delta +t;
		grid_point.y = floor((point.y -t)/delta +0.5)*delta +t;

		// If the last value has been snapped to the same x coordinate, ignore
		if (i > 0) {
			if (prev_point.x == grid_point.x)
				continue;
		}
		prev_point = grid_point;

		math_vector_push(grid_vector, grid_point.y);
	}

	return grid_vector;
}

/* Computes a vector's dimesion after snapping to the specified grid */
static size_t grid_vector_dim(const Math_vector* p, const Grid* grid)
{
	Math_vector* filtered_curve;
	Math_vector* grid_vector;
	size_t dim;

	filtered_curve = curve_filter(p, EPSILON);
	grid_vector    = curve_snap_to_grid(filtered_curve, grid);

	dim = math_vector_dim(grid_vector);

	math_vector_free(grid_vector);
	math_vector_free(filtered_curve);

	return dim;
}

static struct Hashdata hashdata_init(Vector* dataset, int k, int w, Grid* grid)
{
	struct Hashdata data = {0};
	Math_vector* v;
	size_t avg_vector_dim = 0;
	size_t dataset_subset;

	dataset_subset = DATASET_PCT*dataset->size;
	dataset_subset = dataset_subset > 10 ? dataset_subset : 10;

	for (int i = 0; i < dataset_subset; ++i) {
		v = vector_get(dataset, randint(0, dataset->size -1));

		if (grid)
			avg_vector_dim += grid_vector_dim(v, grid);
		else
			avg_vector_dim += math_vector_dim(v);
	}
	avg_vector_dim = avg_vector_dim/dataset_subset;

	data.v = xmalloc(k * sizeof(*data.v));
	data.r = xmalloc(k * sizeof(*data.r));
	data.t = xmalloc(k * sizeof(*data.t));

	data.k = k;
	data.w = w;
	if (grid)
		data.grid = *grid;

	for (int i = 0; i < data.k; ++i) {
		data.t[i] = randfloat(0, data.w);
		data.r[i] = randint(1, ceil(RMULT*w));
		data.v[i] = create_v(avg_vector_dim);
	}

	return data;
}

static void hashdata_free(struct Hashdata* data)
{
	for (int i = 0; i < data->k; ++i)
		math_vector_free(data->v[i]);

	free(data->v);
	free(data->r);
	free(data->t);
}

static Grid grid_create(double delta)
{
	return (Grid) { delta, randfloat(0, delta) };
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
static int calc_w(Vector* dataset, LSH_metric metric)
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

		distance[i] = metric(v, u);
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
	Math_vector* v = math_vector_init(VECTOR_MODE, NULL);

	for (int i = 0; i < dim; ++i)
		math_vector_push(v, uniform_rand(0, 1));

	return v;
}
