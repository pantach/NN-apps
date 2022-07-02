#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include "math.h"
#include "tools.h"

#define VECTOR_GROWTH_FACTOR 2

#include <stdio.h>

struct Math_vector {
	Math_vector_mode mode;
	char* id;
	union {
		double* val;
		Point* point;
	};
	size_t dim;
	size_t capacity;
	void* data;
};

static void k_combinations_recurs(Vector* results, Vector* stack, int n, int k,
                                  int d, int depth);

static double** frechet(const Math_vector* v, const Math_vector* u);

static void math_vector_resize(Math_vector* v, size_t num);

static Math_vector* point_pair_create(const Math_vector* v, int i,
                                      const Math_vector* u, int j);

/* Returns 1 of "num" is an even number, 0 if it is odd */
int even(long num)
{
	return (num % 2 == 0);
}

/* Euclidean modulo */
int mod(int a, int b)
{
	int m = a % b;

	if (m < 0)
		m = (b < 0) ? m - b : m + b;

	return m;
}

double max(int nargs, ...)
{
	va_list args;
	double maxval = -DBL_MAX;
	double val;

	va_start(args, nargs);
	for (int i = 0; i < nargs; i++) {
		val = va_arg(args, double);
		if (val > maxval)
			maxval = val;
	}
	va_end(args);

	return maxval;
}

double min(int nargs, ...)
{
	va_list args;
	double minval = DBL_MAX;
	double val;

	va_start(args, nargs);
	for (int i = 0; i < nargs; i++) {
		val = va_arg(args, double);
		if (val < minval)
			minval = val;
	}
	va_end(args);

	return minval;
}

unsigned int min_index(int nargs, ...)
{
	va_list args;
	double minval = DBL_MAX;
	double val;
	int index = -1;

	va_start(args, nargs);
	for (int i = 0; i < nargs; i++) {
		val = va_arg(args, double);
		if (val < minval) {
			minval = val;
			index = i;
		}
	}
	va_end(args);

	return index;
}

/* Generates a random number following a uniform distribution */
double uniform_rand(double mu, double sigma)
{
	double U1, U2, W, mult;
	static double X1, X2;
	static int call = 0;

	if (call == 1) {
		call = !call;
		return (mu + sigma * (double) X2);
	}

	do {
		U1 = -1 + ((double) rand () / RAND_MAX) * 2;
		U2 = -1 + ((double) rand () / RAND_MAX) * 2;
		W = pow (U1, 2) + pow (U2, 2);
	} while (W >= 1 || W == 0);

	mult = sqrt ((-2 * log (W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;

	return (mu + sigma * (double) X1);
}

/* Generates a random float number between min and max. Max should be >= min */
float randfloat(float min, float max)
{
	float num;

	max -= min;
	num = (((double)rand()/(double)(RAND_MAX)) * max) +min;

	return num;
}

/* Generates a random integer number between min and max. Max should be >= min */
int randint(int min, int max)
{
	int num;

	max -= min;
	num = rand()%(max+1) + min;

	return num;
}

/* Generates a random integer */
int randZ(void)
{
	return randint(-RAND_MAX/2, +RAND_MAX/2);
}

/* Generates and returns all the k-combinations c(n, k) as an array of
 * coordinate vectors
 * */
Vector* k_combinations(int n, int k)
{
	Vector* results = NULL;

	if (n >= k && k >= 0) {
		results = vector_init(math_vector_free_generic);

		Vector* stack = vector_init(NULL);

		k_combinations_recurs(results, stack, n, k, 0, 0);
		vector_free(stack);
	}

	return results;
}

static void k_combinations_recurs(Vector* results, Vector* stack, int n, int k,
                                  int d, int depth)
{
	int i;

	if (depth == k) {
		Math_vector* v = math_vector_init(VECTOR_MODE, NULL);

		for (i = 0; i < vector_size(stack); i++)
			math_vector_push(v, (double)(int)(intptr_t) vector_get(stack, i));

		vector_push(results, v);
		return;
	}

	for (i = d; i < n; i++) {
		if (i >= depth) {
			vector_push(stack, (void*)(intptr_t)i);
			k_combinations_recurs(results, stack, n, k, i +1, depth +1);
			vector_pop(stack);
		}
	}
}

double cartesian_distance(Point p, Point q)
{
	return sqrt(pow(p.x -q.x, 2) +pow(p.y -q.y, 2));
}

Point midpoint(Point p, Point q)
{
	return (Point) { (p.x +q.x)/2, (p.y +q.y)/2 };
}

/* Checks whether two vectors are euqual */
int vequal(const Math_vector* v, const Math_vector* u)
{
	if (math_vector_dim(v) != math_vector_dim(u))
		return 0;

	for (int i = 0; i < v->dim; ++i)
		if (math_vector_val(v, i) != math_vector_val(u, i))
			return 0;

	return 1;
}

/* Calculates the dot product between two vectors */
double dot(const Math_vector* v, const Math_vector* u)
{
	size_t dim = (v->dim > u->dim) ? u->dim : v->dim;
	double sum = 0;

	for (int i = 0; i < dim; ++i)
		sum += v->val[i] * u->val[i];

	return sum;
}

/* Calculates the euclidean distance between two vectors */
double euclidean_distance(const Math_vector* v, const Math_vector* u)
{
	size_t dim = (v->dim > u->dim) ? u->dim : v->dim;
	double sum = 0;

	for (int i = 0; i < dim; ++i)
		sum += pow(v->val[i] - u->val[i], 2);

	return sqrt(sum);
}

/* Calculates the square of the euclidean distance between two vectors */
double euclidean_distance_sq(const Math_vector* v, const Math_vector* u)
{
	size_t dim = (v->dim > u->dim) ? u->dim : v->dim;
	double sum = 0;

	for (int i = 0; i < dim; ++i)
		sum += pow(v->val[i] - u->val[i], 2);

	return sum;
}

/* Calculates the norm of a vector */
double norm(const Math_vector* v)
{
	double sum = 0;

	for (int i = 0; i < v->dim; ++i)
		sum += pow(v->val[i], 2);

	return sqrt(sum);
}

void vsum(Math_vector* sum, const Math_vector* v)
{
	int i;

	for (i = 0; i < sum->dim; i++)
		sum->val[i] += v->val[i];

	for (i = sum->dim; i < v->dim; i++)
		math_vector_push(sum, v->val[i]);
}

Math_vector* vadd(const Math_vector* v, const Math_vector* u)
{
	Math_vector* sum = math_vector_init(VECTOR_MODE, 0);
	const Math_vector* higher;
	const Math_vector* lower;
	int i;

	if (v->dim >= u->dim) {
		higher = v;
		lower  = u;
	}
	else {
		higher = u;
		lower  = v;
	}

	for (i = 0; i < lower->dim; i++)
		math_vector_push(sum, lower->val[i] +higher->val[i]);

	for (i = lower->dim; i < higher->dim; i++)
		math_vector_push(sum, higher->val[i]);

	return sum;
}

Math_vector* vmult(const Math_vector* v, double scalar)
{
	Math_vector* product = math_vector_init(VECTOR_MODE, 0);

	for (int i = 0; i < v->dim; i++)
		math_vector_push(product, scalar * v->val[i]);

	return product;
}

Math_vector* math_vector_init(Math_vector_mode mode, char* id)
{
	Math_vector* v = xcalloc(1, sizeof(*v));

	v->mode = mode;
	v->id = xstrdup(id);

	return v;
}

void math_vector_free(Math_vector* v)
{
	if (v) {
		if (v->mode == VECTOR_MODE)
			free(v->val);
		else if (v->mode == CURVE_MODE)
			free(v->point);

		free(v->id);
		free(v);
	}
}

static void math_vector_resize(Math_vector* v, size_t num)
{
	if (v->mode == VECTOR_MODE)
		v->val = xrealloc(v->val, num*sizeof(*v->val));
	else if (v->mode == CURVE_MODE)
		v->point = xrealloc(v->point, num*sizeof(*v->point));

	v->capacity = num;
}

void math_vector_push(Math_vector* v, ...)
{
	va_list args;

	if (v->dim == v->capacity)
		math_vector_resize(v, (v->capacity +1)*VECTOR_GROWTH_FACTOR);

	va_start(args, v);
	if (v->mode == VECTOR_MODE) {
		double val = va_arg(args, double);
		v->val[v->dim] = val;
	}
	else if (v->mode == CURVE_MODE) {
		Point val = va_arg(args, Point);
		v->point[v->dim] = val;
	}
	va_end(args);

	v->dim++;
}

double math_vector_val(const Math_vector* v, int i)
{
	return v->val[i];
}

Point math_vector_point(const Math_vector* v, int i)
{
	return v->point[i];
}

size_t math_vector_dim(const Math_vector* v)
{
	return v->dim;
}

char* math_vector_id(const Math_vector* v)
{
	if (v)
		return v->id;

	return "null";
}

void* math_vector_data(const Math_vector* v)
{
	return v->data;
}

void math_vector_set_val(Math_vector* v, int i, ...)
{

	va_list args;

	va_start(args, i);
	if (v->mode == VECTOR_MODE) {
		double val = va_arg(args, double);
		v->val[i] = val;
	}
	else if (v->mode == CURVE_MODE) {
		Point val = va_arg(args, Point);
		v->point[i] = val;
	}
	va_end(args);
}

void math_vector_set_data(Math_vector* v, void* data)
{
	v->data = data;
}

int math_vector_compid(const Math_vector* v, const Math_vector* u)
{
	return strcmp(v->id, u->id);
}

Math_vector* math_vector_dup(const Math_vector* v)
{
	Math_vector* dup = math_vector_init(v->mode, v->id);

	for (int i = 0; i < v->dim; i++) {
		if (v->mode == VECTOR_MODE)
			math_vector_push(dup, v->val[i]);
		else if (v->mode == CURVE_MODE)
			math_vector_push(dup, v->point[i]);
	}

	dup->data = v->data;

	return dup;
}

void math_vector_free_generic(void* v)
{
	math_vector_free(v);
}

char* math_vector_printstr(const Math_vector* v)
{
	char* str = NULL;
	char* buf;
	int i;

	if (v->mode == VECTOR_MODE) {
		for (i = 0; i < v->dim; i++) {
			if (i == v->dim -1)
				xsprintf(&buf, "%.2f", v->val[i]);
			else
				xsprintf(&buf, "%.2f, ", v->val[i]);

			xstrcat(&str, buf);
			free(buf);
		}
	}
	else if (v->mode == CURVE_MODE) {
		for (i = 0; i < v->dim; i++) {
			if (i == v->dim -1)
				xsprintf(&buf, "(%.2f, %.2f)", v->point[i].x, v->point[i].y);
			else
				xsprintf(&buf, "(%.2f, %.2f), ", v->point[i].x, v->point[i].y);

			xstrcat(&str, buf);
			free(buf);
		}
	}

	return str;
}

static double** frechet(const Math_vector* v, const Math_vector* u)
{
	double** C;
	Point p, q;
	int i, j;

	// Allocate dynamic memory
	C = xmalloc(v->dim * sizeof(double));
	for (i = 0; i < v->dim; i++)
		C[i] = xmalloc(u->dim * sizeof(double));

	// Execute the DFD dynamic programming algorithm
	for (i = 0; i < v->dim; i++) {
		for (j = 0; j < u->dim; j++) {
			p = v->point[i];
			q = u->point[j];

			if (i == 0 && j == 0)
				C[0][0] = cartesian_distance(p, q);
			else if (i == 0)
				C[0][j] = MAX(C[0][j -1], cartesian_distance(p, q));
			else if (j == 0)
				C[i][0] = MAX(C[i -1][0], cartesian_distance(p, q));
			else
				C[i][j] = MAX(MIN(C[i -1][j], C[i -1][j -1], C[i][j -1]),
				              cartesian_distance(p, q));
		}
	}

	return C;
}

double frechet_distance(const Math_vector* v, const Math_vector* u)
{
	double** C = frechet(v, u);
	double dist;

	dist = C[v->dim -1][u->dim -1];

	// Free memory
	for (int i = 0; i < v->dim; i++)
		free(C[i]);
	free(C);

	return dist;
}

static Math_vector* point_pair_create(const Math_vector* v, int i,
                                      const Math_vector* u, int j)
{
	Math_vector* point_pair;
	Point p;
	Point q;

	point_pair = math_vector_init(CURVE_MODE, NULL);

	p = v->point[i];
	q = u->point[j];

	math_vector_push(point_pair, p);
	math_vector_push(point_pair, q);

	return point_pair;
}


Vector* frechet_traversal(const Math_vector* v, const Math_vector* u)
{
	double** C = frechet(v, u);
	Vector* pairs;
	Math_vector* point_pair;
	unsigned int minidx;
	int i, j;

	pairs = vector_init(math_vector_free_generic);

	i = v->dim -1;
	j = u->dim -1;

	point_pair = point_pair_create(v, i, u, j);
	vector_prepend(pairs, point_pair);

	while (i +j != 0) {

		if (i -1 == -1)
			minidx = 1;
		else if (j -1 == -1)
			minidx = 0;
		else
			minidx = MIN_INDEX(C[i -1][j], C[i][j -1], C[i -1][j -1]);

		if (minidx == 0)
			i--;
		else if (minidx == 1)
			j--;
		else {
			i--;
			j--;
		}
		//printf("minidx: %u\n", minidx);
		//printf("%d, %d\n", i, j);

		point_pair = point_pair_create(v, i, u, j);
		vector_prepend(pairs, point_pair);
	}

	// Free memory
	for (i = 0; i < v->dim; i++)
		free(C[i]);
	free(C);

	return pairs;
}

Math_vector* frechet_mean(const Math_vector* v, const Math_vector* u)
{
	Vector* pairs;
	Math_vector* mean_curve;
	Math_vector* point_pair;
	Point p;
	Point q;

	if (v == NULL && u == NULL)
		return NULL;
	if (v == NULL)
		return math_vector_dup(u);
	if (u == NULL)
		return math_vector_dup(v);

	pairs = frechet_traversal(v, u);

	mean_curve = math_vector_init(CURVE_MODE, NULL);

	for (int i = 0; i < pairs->size; i++) {
		point_pair = vector_get(pairs, i);

		p = point_pair->point[0];
		q = point_pair->point[1];

		math_vector_push(mean_curve, midpoint(p, q));
	}

	vector_free(pairs);

	return mean_curve;
}

/*
int main(void)
{
	Math_vector* v = math_vector_init(CURVE_MODE, NULL);
	char* str;

	for (int i = 0; i < 20; i++)
		math_vector_push(v, (Point) { randint(0, 100), randfloat(-20, 20) });

	str = math_vector_printstr(v);

	printf("%s\n", str);

	free(str);
	math_vector_free(v);

	return 0;
}
*/
