#ifndef LSH_H
#define LSH_H

#include <stdarg.h>
#include "math.h"
#include "vector.h"

typedef struct LSH LSH;
typedef double (*LSH_metric)(const Math_vector*, const Math_vector*);

typedef enum {
	LSH_MODE_VECTOR = 0,
	LSH_MODE_CURVE
} LSH_mode;

LSH* LSH_init(LSH_mode mode, LSH_metric metric, Vector* dataset, size_t ht_size,
              size_t L, size_t k, ...);
void LSH_free(LSH* lsh);
Math_vector* LSH_NN(LSH* lsh, const Math_vector* q);
Vector* LSH_kNN(LSH* lsh, const Math_vector* q, int N);
Vector* LSH_range_search(LSH* lsh, const Math_vector* q, float radius);

double LSH_kNN_time(const LSH* lsh);
double LSH_exhaust_time(const LSH* lsh);

Math_vector* curve_filter(const Math_vector* curve, double epsilon);
#endif
