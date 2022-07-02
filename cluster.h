#ifndef CLUSTER_H
#define CLUSTER_H

#include <stddef.h>
#include <stdbool.h>
#include "math.h"
#include "vector.h"

typedef struct Clustering Clustering;
typedef double (*Clustering_metric)(const Math_vector*, const Math_vector*);

typedef enum {
	ASSIGN_METHOD_CLASSIC = 0,
	ASSIGN_METHOD_RLSH,
	ASSIGN_METHOD_RCUBE,
	ASSIGN_METHOD_FRECHET,
} Assign_method;

typedef enum {
	UPDATE_METHOD_VECTOR = 0,
	UPDATE_METHOD_CURVE
} Update_method;

Clustering* clustering_init(Assign_method amethod, Update_method umethod,
                            const char* datapath, const char* confpath,
                            const char* outpath);
void clustering_free(Clustering* c);

void clustering_init_pp(Clustering* c);
void clustering_exec(Clustering* c);

void clustering_print(Clustering* c, bool complete);

#endif
