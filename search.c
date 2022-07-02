#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "lsh.h"
#include "cube.h"
#include "math.h"
#include "tools.h"
#include "io.h"

#define MAX_PATH 512

typedef enum {
	PARSE_ALGORITHM = 0,
	PARSE_FRECHET_MODE
} Parse_opt;

typedef enum {
	ALGORITHM_LSH = 0,
	ALGORITHM_CUBE,
	ALGORITHM_FRECHET
} Algorithm;

typedef enum {
	FRECHET_DISCRETE = 0,
	FRECHET_CONTINUOUS
} Frechet_mode;

// Command line arguments
struct CLA {
	Algorithm alg;
	Frechet_mode fmode;
	Math_vector_mode vmode;
	char input_file[MAX_PATH];
	char output_file[MAX_PATH];
	char query_file[MAX_PATH];
	int k;
	int L;
	int M;
	int probes;
	double delta;
};

typedef double (*Metric)(const Math_vector*, const Math_vector*);

static void parse_cla(int argc, char** argv);
static void parse_opt(Parse_opt opt, const char* input);
static void check_cla(void);
static char* algname(Algorithm alg, Frechet_mode fmode);
static void env_init(void);
static void read_usr_input(void);
static int  prompt_rerun(void);
static Math_vector* dataset_NN(Vector* dataset, const Math_vector* q,
                               Metric metric, double* time);

// Globals
struct CLA g_cla = { .alg = -1, .fmode = -1, .k = -1, .L = 5, .M = 10,
                     .probes = 2, .delta = 1 };

static void print_usage(void)
{
	err_exit("./search –i <input file> –q <query file> –k <int> -L <int> "
	         " -M <int> -probes <int> -ο <output file> "
	         " -algorithm <LSH or Hypercube or Frechet> "
	         " -metric <discrete or continuous | only for –algorithm Frechet> "
	         " -delta <double>");
}

static void parse_cla(int argc, char** argv)
{
	for (int i = 1; i < argc; ++i) {
		if (!strcmp(argv[i], "-i"))
			strncpy(g_cla.input_file, argv[++i], MAX_PATH -1);

		else if (!strcmp(argv[i], "-o"))
			strncpy(g_cla.output_file, argv[++i], MAX_PATH -1);

		else if (!strcmp(argv[i], "-q"))
			strncpy(g_cla.query_file, argv[++i], MAX_PATH -1);

		else if (!strcmp(argv[i], "-k"))
			g_cla.k = getint(argv[++i], 0);

		else if (!strcmp(argv[i], "-L"))
			g_cla.L = getint(argv[++i], 0);

		else if (!strcmp(argv[i], "-M"))
			g_cla.M = getint(argv[++i], 0);

		else if (!strcmp(argv[i], "-probes"))
			g_cla.probes = getint(argv[++i], 0);

		else if (!strcmp(argv[i], "-delta"))
			g_cla.delta = getfloat(argv[++i], 0);

		else if (!strcmp(argv[i], "-algorithm"))
			parse_opt(PARSE_ALGORITHM, argv[++i]);

		else if (!strcmp(argv[i], "-metric"))
			parse_opt(PARSE_FRECHET_MODE, argv[++i]);

		else {
			fprintf(stderr, "Unknown argument: %s\n", argv[i]);
			print_usage();
		}
	}
}

static void check_cla(void)
{
	if (g_cla.k == -1) {
		if (g_cla.alg == ALGORITHM_LSH || g_cla.alg == ALGORITHM_FRECHET)
			g_cla.k = 4;

		else if (g_cla.alg == ALGORITHM_CUBE)
			g_cla.k = 14;
	}
}

static void parse_opt(Parse_opt opt, const char* input)
{
	switch (opt) {
		case PARSE_ALGORITHM:
			if (!strcasecmp(input, "LSH")) {
				g_cla.alg = ALGORITHM_LSH;
				g_cla.vmode = VECTOR_MODE;
			}

			else if (!strcasecmp(input, "Hypercube")) {
				g_cla.alg = ALGORITHM_CUBE;
				g_cla.vmode = VECTOR_MODE;
			}

			else if (!strcasecmp(input, "Frechet")) {
				g_cla.alg = ALGORITHM_FRECHET;
				g_cla.vmode = CURVE_MODE;
			}

			else {
				fprintf(stderr, "Unknown algorithm: %s\n", input);
				print_usage();
			}
			break;
		case PARSE_FRECHET_MODE:
			if (!strcasecmp(input, "discrete"))
				g_cla.fmode = FRECHET_DISCRETE;

			else if (!strcasecmp(input, "continuous"))
				g_cla.fmode = FRECHET_CONTINUOUS;

			else {
				fprintf(stderr, "Unknown metric: %s\n", input);
				print_usage();
			}
			break;
	}
}

static void read_usr_input(void)
{
	char input[MAX_PATH];

	if (!strcmp(g_cla.input_file, "")) {
		printf("\nPlease write the path to the input file\n");
		scanf("%s", g_cla.input_file);
	}

	if (!strcmp(g_cla.query_file, "")) {
		printf("\nPlease write the path to the query file\n");
		scanf("%s", g_cla.query_file);
	}

	if (!strcmp(g_cla.output_file, "")) {
		printf("\nPlease write the path to the output file\n");
		scanf("%s", g_cla.output_file);
	}

	if (g_cla.alg == -1) {
		fprintf(stderr, "Please specify the algotithm\n");
		scanf("%s", input);
		parse_opt(PARSE_ALGORITHM, input);
	}

	if (g_cla.alg == ALGORITHM_FRECHET && g_cla.fmode == -1) {
		fprintf(stderr, "Please specify the metric\n");
		scanf("%s", input);
		parse_opt(PARSE_FRECHET_MODE, input);
	}

}

static char* algname(Algorithm alg, Frechet_mode fmode)
{
	switch (alg) {
	case ALGORITHM_LSH:
		return "LSH_Vector";

	case ALGORITHM_CUBE:
		return "Hypercube";

	case ALGORITHM_FRECHET:
		switch (fmode) {
		case FRECHET_DISCRETE:
			return "LSH_Frechet_Discrete";

		case FRECHET_CONTINUOUS:
			return "LSH_Frechet_Continuous";
		}
	}

	return "<Name error>";
}

static int prompt_rerun(void)
{
	char answer = 'a';

	printf("\nWould you like a rerun? (y/n)\n");
	scanf(" %c", &answer);
	while (answer != 'y' && answer != 'n') {
		printf("Please type 'y' or 'n'\n");
		scanf(" %c", &answer);
	}
	if (answer == 'y') {
		g_cla.input_file[0]  = '\0';
		g_cla.output_file[0] = '\0';
		g_cla.query_file[0]  = '\0';

		return 1;
	}
	else
		return 0;
}

static void env_init(void)
{
	srand(time(NULL));
}

int main(int argc, char** argv)
{
	Vector* dataset;
	Vector* queries;
	FILE* output;
	LSH* lsh = NULL;
	Cube* cube = NULL;
	Metric metric;
	size_t htab_size;
	double true_search_time;
	double avg_true_search_time = 0;
	double avg_appr_search_time = 0;
	double AF, MAF = 0;
	int i;

	env_init();
	parse_cla(argc, argv);
	check_cla();

	do {
		read_usr_input();

		printff("\nParsing input files...");

		dataset = read_dataset(g_cla.input_file, g_cla.vmode);
		if (dataset->size == 0)
			err_exit("Empty dataset\n");

		queries = read_dataset(g_cla.query_file, g_cla.vmode);
		if (queries->size == 0)
			err_exit("Empty query file\n");

		printff("Done\n");

		output = fopen(g_cla.output_file, "wb");
		if (!output)
			syserr_exit("Could not open \"%s\"", g_cla.output_file);

		htab_size = (dataset->size/8 >= 8) ? dataset->size/8 : 8;

		printff("Initializing hashtables...");
		switch (g_cla.alg) {
		case ALGORITHM_LSH:
			metric = euclidean_distance;
			lsh = LSH_init(LSH_MODE_VECTOR, metric, dataset, htab_size, g_cla.L,
			               g_cla.k);
			break;

		case ALGORITHM_FRECHET:
			metric = frechet_distance;
			lsh = LSH_init(LSH_MODE_CURVE, metric, dataset, htab_size, g_cla.L,
			               g_cla.k, g_cla.delta);
			break;

		case ALGORITHM_CUBE:
			metric = euclidean_distance;
			cube = cube_init(dataset, htab_size, g_cla.k);
			cube_setopt(cube, CUBE_METRIC, metric);
			cube_setopt(cube, CUBE_PROBES, g_cla.probes);
			cube_setopt(cube, CUBE_M, g_cla.M);
			break;
		}
		printff("Done\n");

		printff("Querying hashtables...");
		for (i = 0; i < vector_size(queries); ++i) {

			Math_vector* q;
			Math_vector* appr_neighbor;
			Math_vector* true_neighbor;
			double appr_dist;
			double true_dist;

			q = vector_get(queries, i);

			printf("\r%60s", "");
			printff("\rQuerying hashtables...query: %s", math_vector_id(q));

			// NN search
			switch (g_cla.alg) {
			case ALGORITHM_LSH:
			case ALGORITHM_FRECHET:
				appr_neighbor = LSH_NN(lsh, q);
				avg_appr_search_time += LSH_kNN_time(lsh);
				break;

			case ALGORITHM_CUBE:
				appr_neighbor = cube_NN(cube, q);
				avg_appr_search_time += cube_kNN_time(cube);
				break;
			}

			// Exhaustive search
			true_neighbor = dataset_NN(dataset, q, metric, &true_search_time);
			avg_true_search_time += true_search_time;

			// Calculate respective distances
			if (appr_neighbor)
				appr_dist = metric(q, appr_neighbor);
			else
				appr_dist = INFINITY;

			true_dist = metric(q, true_neighbor);

			// Calculate approximation factor, avoid division by zero
			if (true_dist == 0) {
				if (appr_dist == 0)
					AF = 1;
				else
					AF = INFINITY;
			}
			else
				AF = appr_dist/true_dist;

			// Calculate maximum approximation factor
			if (AF > MAF) MAF = AF;

			fprintf(output,
			        "Query: %s\n"
			        "Algorithm: %s\n"
			        "Approximate Nearest neighbor: %s\n"
			        "True neighbor: %s\n"
			        "distanceApproximate: %f\n"
			        "distanceTrue: %f\n\n",
			        math_vector_id(q),
			        algname(g_cla.alg, g_cla.fmode),
			        math_vector_id(appr_neighbor),
			        math_vector_id(true_neighbor),
			        appr_dist,
			        true_dist);
		}

		// Calculate average times
		avg_appr_search_time = avg_appr_search_time/vector_size(queries);
		avg_true_search_time = avg_true_search_time/vector_size(queries);

		fprintf(output,
		        "\n"
		        "tApproximateAverage: %f\n"
		        "tTrueAverage: %f\n"
				"MAF: %f\n",
				avg_appr_search_time,
				avg_true_search_time,
				MAF);

		printff("\rQuerying hashtables...Done%16s\n", " ");

		if (output) fclose(output);

		LSH_free(lsh);
		cube_free(cube);

		free_dataset(queries);
		free_dataset(dataset);

	} while (prompt_rerun());

	printf("\nDone!\n");
	return 0;
}

static Math_vector* dataset_NN(Vector* dataset, const Math_vector* q,
							   Metric metric, double* time)
{
	Math_vector* NN = NULL;
	Math_vector* v;
	double min_dist = DBL_MAX;
	double dist;
	clock_t start, end;

	start = clock();
	for (int i = 0; i < dataset->size; i++) {
		v = vector_get(dataset, i);
		dist = metric(q, v);
		if (dist < min_dist) {
			min_dist = dist;
			NN = v;
		}
	}
	end = clock();
	*time = (double) (end -start)/CLOCKS_PER_SEC;

	return NN;
}
