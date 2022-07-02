#ifndef IO_H
#define IO_H

#include "vector.h"
#include "math.h"

typedef struct {
	int K;
	int L;
	int klsh;
	int M;
	int kcube;
	int probes;
} Cluster_conf;

Vector* read_dataset(const char* input_file, Math_vector_mode mode);
void free_dataset(Vector* dataset);
Cluster_conf read_cluster_conf(const char* filepath);

Vector* tokenize(const char* str, const char* delim);

#endif
