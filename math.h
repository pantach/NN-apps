#ifndef MATH_H
#define MATH_H

#include <stdarg.h>
#include <stdbool.h>
#include "vector.h"

#define NARGS(...) (sizeof((int[]){__VA_ARGS__}) / sizeof(int))
#define MAX(...) max(NARGS(__VA_ARGS__), __VA_ARGS__)
#define MIN(...) min(NARGS(__VA_ARGS__), __VA_ARGS__)
#define MIN_INDEX(...) min_index(NARGS(__VA_ARGS__), __VA_ARGS__)

typedef struct Math_vector Math_vector;
typedef struct Curve Curve;

typedef struct {
	double x;
	double y;
} Point;

typedef enum {
	VECTOR_MODE = 0,
	CURVE_MODE
} Math_vector_mode;

/* Generic */
int even(long num);
int mod(int a, int b);
double max(int nargs, ...);
double min(int nargs, ...);
unsigned int min_index(int nargs, ...);

/* Random numbers */
double uniform_rand(double mu, double sigma);
float  randfloat(float min, float max);
int    randint(int min, int max);
int    randZ(void);

/* Combinatorics */
Vector* k_combinations(int n, int k);

/* Point operations */
double cartesian_distance(Point p, Point q);
Point  midpoint(Point p, Point q);

/* Vector operations */
int    vequal(const Math_vector* v, const Math_vector* u);
double dot(const Math_vector* v, const Math_vector* u);
double euclidean_distance(const Math_vector* v, const Math_vector* u);
double euclidean_distance_sq(const Math_vector* v, const Math_vector* u);
double norm(const Math_vector* v);
void   vsum(Math_vector* sum, const Math_vector* v);
Math_vector* vadd(const Math_vector* v, const Math_vector* u);
Math_vector* vmult(const Math_vector* v, double scalar);

/* Curve operations */
double frechet_distance(const Math_vector* v, const Math_vector* u);
Math_vector* frechet_mean(const Math_vector* v, const Math_vector* u);

/* Vector creation and use */
Math_vector* math_vector_init(Math_vector_mode mode, char* id);
void   math_vector_free(Math_vector* v);
void   math_vector_free_generic(void* v);
void   math_vector_push(Math_vector* v, ...);
double math_vector_val(const Math_vector* v, int i);
Point  math_vector_point(const Math_vector* v, int i);
size_t math_vector_dim(const Math_vector* v);
char*  math_vector_id(const Math_vector* v);
void*  math_vector_data(const Math_vector* v);
void   math_vector_set_val(Math_vector* v, int i, ...);
void   math_vector_set_data(Math_vector* v, void* data);
char*  math_vector_printstr(const Math_vector* v);
int    math_vector_compid(const Math_vector* v, const Math_vector* u);
Math_vector* math_vector_dup(const Math_vector* v);

#endif
