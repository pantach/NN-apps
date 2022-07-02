#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"
#include "tools.h"

static void free_dataset_entry(void* entry);

Vector* read_dataset(const char* input_file, Math_vector_mode mode)
{
	Vector* dataset;
	FILE*   fp;
	char*   line = NULL;
	char*   line_id;
	size_t  line_size = 0;
	Vector* line_field;
	Math_vector* v;

	dataset = vector_init(free_dataset_entry);

	fp = fopen(input_file, "rb");
	if (!fp)
		syserr_exit("Cannot open \"%s\"", input_file);

	for (int i = 0; getline(&line, &line_size, fp) != -1; ++i) {
		line_field = tokenize(line, " \t\r\n");
		line_id    = vector_get(line_field, 0);

		v = math_vector_init(mode, line_id);

		for (int j = 1; j < line_field->size; ++j) {
			double val = getfloat(vector_get(line_field, j), 0);

			if (mode == VECTOR_MODE)
				math_vector_push(v, val);

			else if (mode == CURVE_MODE)
				math_vector_push(v, (Point) { j-1, val });
		}

		vector_push(dataset, v);
		vector_free(line_field);
	}

	free(line);
	fclose(fp);

	return dataset;
}

void free_dataset(Vector* dataset)
{
	vector_free(dataset);
}

static void free_dataset_entry(void* entry)
{
	math_vector_free(entry);
}

Cluster_conf read_cluster_conf(const char* filepath)
{
	Cluster_conf params = { .K = -1, .L = 3, .klsh = 4, .M = 10, .kcube = 3,
	                        .probes = 2 };
	FILE*   fp;
	char*   line = NULL;
	size_t  line_size = 0;
	Vector* line_field;

	fp = fopen(filepath, "rb");
	if (!fp)
		syserr_exit("Cannot open \"%s\"", filepath);

	for (int i = 0; getline(&line, &line_size, fp) != -1; ++i) {
		line_field = tokenize(line, " :\r\n");

		const char* const param_name = vector_get(line_field, 0);
		const char* const param      = vector_get(line_field, 1);

		if (!strcmp(param_name, "number_of_clusters"))
			params.K = getint(param, 0);

		else if (!strcmp(param_name, "number_of_vector_hash_tables"))
			params.L = getint(param, 0);

		else if (!strcmp(param_name, "number_of_vector_hash_functions"))
			params.klsh = getint(param, 0);

		else if (!strcmp(param_name, "max_number_M_hypercube"))
			params.M = getint(param, 0);

		else if (!strcmp(param_name, "number_of_hypercube_dimensions"))
			params.kcube = getint(param, 0);

		else if (!strcmp(param_name, "number_of_probes"))
			params.probes = getint(param, 0);

		else
			err_exit("Error reading \"%s\"\nUnrecognized parameter name: %s\n",
			         param_name);

		vector_free(line_field);
	}

	if (params.K == -1)
		err_exit("Please specify the number of clusters\n");

	free(line);
	fclose(fp);

	return params;
}

Vector* tokenize(const char* str, const char* delim)
{
	Vector* v;
	char* strcp;
	char* pch;

	v = vector_init(free);
	if (!v)
		abort();

	strcp = xstrdup(str);

	pch = strtok(strcp, delim);
	while (pch) {
		vector_push(v, xstrdup(pch));
		pch = strtok(NULL, delim);
	}

	free(strcp);

	return v;
}
