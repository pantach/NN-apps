#ifndef CUBE_H
#define CUBE_H

#include "math.h"
#include "vector.h"

typedef struct Cube Cube;
typedef double (*Cube_metric)(Math_vector*, Math_vector*);

typedef enum {
	CUBE_METRIC = 0,
	CUBE_PROBES,
	CUBE_M
} Cube_opt;

Cube* cube_init(Vector* dataset, size_t ht_size, size_t k);
void cube_setopt(Cube* cube, Cube_opt opt, ...);
Math_vector* cube_NN(Cube* cube, Math_vector* q);
Vector* cube_kNN(Cube* cube, Math_vector* q, int N);
Vector* cube_range_search(Cube* cube, Math_vector* q, float radius);
Math_vector* cube_exhaust_search(Cube* cube, Math_vector* q);
void cube_free(Cube* cube);

double cube_kNN_time(Cube* cube);
double cube_exhaust_time(Cube* cube);

#endif
