#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include "cluster.h"
#include "lsh.h"
#include "cube.h"
#include "tools.h"
#include "io.h"
#include "progbar.h"

#define MAX_PATH 512
#define MAX_ITERATIONS 15
#define CLUSTERING_PRECISION_PCT 0.01
#define EPSILON 0.1        // Curve filtering allowance

struct Clustering {
	Assign_method amethod;
	Update_method umethod;
	Clustering_metric metric;
	Cluster_conf params;
	Vector* dataset;
	FILE* output;
	Vector* clusters;
	Vector* centroids;
	Vector* silval;
	double time;
};

// Command line arguments
struct CLA {
	Assign_method amethod;
	Update_method umethod;
	char input_file [MAX_PATH];
	char output_file[MAX_PATH];
	char config_file[MAX_PATH];
	bool silhouette;
	bool complete;
};

typedef enum {
	PARSE_ASSIGN_METHOD = 0,
	PARSE_UPDATE_METHOD
} Parse_opt;

typedef struct {
	Math_vector* centroid;
	Vector* points;
} Cluster;

typedef struct {
	Cluster* cluster;
	Cluster* mark;
} Point_data;

static void env_init(void);

static void parse_cla(int argc, char** argv);
static int  parse_opt(Parse_opt opt, const char* method);
static void check_cla(void);
static void read_usr_input(void);
static int  prompt_rerun(void);

static void     clustering_create_cluster(Clustering* c, Math_vector* centroid);
static void     clustering_clear_points(Clustering* c);
static void     clustering_clear_marks(Clustering* c);
static Cluster* clustering_assign_to_nearest(Clustering* c, Math_vector* point);
static Vector*  clustering_centroids(Clustering* c);
static void     clustering_update_centroids(Clustering* c);
static double   clustering_min_centroid_dist(Clustering* c);
static Cluster* clustering_second_best_cluster(Clustering* c,
				   const Cluster* cluster, const Math_vector* point);

static void     clustering_lloyds(Clustering* c);
static void     clustering_reverse_LSH(Clustering* c);
static void     clustering_reverse_cube(Clustering* c);
static void     clustering_silhouette(Clustering* c);

static Cluster* cluster_init(Math_vector* centroid);
static void     cluster_free(Cluster* cluster);
static void     cluster_free_generic(void* cluster);
static void     cluster_clear(Cluster* cluster);
static void     cluster_replace_centroid(Cluster* cluster, Math_vector* centroid);
static Math_vector* cluster_mean_vector(Cluster* cluster);
static Math_vector* cluster_mean_curve(Cluster* cluster);
static char*    cluster_printstr(Cluster* cluster, bool complete);

static Cluster* mark_of(Math_vector* point);
static Cluster* cluster_of(Math_vector* point);
static void     mark(Math_vector* point, Cluster* cluster);
static Cluster* cluster_assign(Cluster* cluster, Math_vector* point);

static double min_dist(const Math_vector* v, Vector* points,
                       Clustering_metric metric);
static int mark_points(const Vector* points, Cluster* cluster,
                       Clustering_metric metric);
static int centroid_comp(const void* x_, const void* val_, void* unused);

// Globals
struct CLA g_cla = { .amethod = -1, .umethod = -1, .complete = false,
                     .silhouette = false };

static void print_usage(const char* exename)
{
	err_exit("Usage:\n%s –i <input file> –c <configuration file> "
	         "-o <output file> -assignement <Classic OR LSH or Hypercube or "
	         "LSH_Frechet> -update <\"Mean Frechet\" or \"Mean Vector\">"
	         "-complete <optional> -silhouette <optional>\n"
	         ,exename);
}

static void parse_cla(int argc, char** argv)
{
	int err;

	for (int i = 1; i < argc; ++i) {
		if (!strcmp(argv[i], "-i"))
			strncpy(g_cla.input_file, argv[++i],  MAX_PATH -1);

		else if (!strcmp(argv[i], "-o"))
			strncpy(g_cla.output_file, argv[++i], MAX_PATH -1);

		else if (!strcmp(argv[i], "-c"))
			strncpy(g_cla.config_file, argv[++i], MAX_PATH -1);

		else if (!strcmp(argv[i], "-complete"))
			g_cla.complete = true;

		else if (!strcmp(argv[i], "-silhouette"))
			g_cla.silhouette = true;

		else if (!strcmp(argv[i], "-assignment")) {
			err = parse_opt(PARSE_ASSIGN_METHOD, argv[++i]);
			if (err)
				print_usage(argv[0]);
		}

		else if (!strcmp(argv[i], "-update")) {
			err = parse_opt(PARSE_UPDATE_METHOD, argv[++i]);
			if (err)
				print_usage(argv[0]);
		}

		else {
			fprintf(stderr, "Unknown argument: %s\n", argv[i]);
			print_usage(argv[0]);
		}
	}

	if (g_cla.amethod == -1) {
		fprintf(stderr, "Please specify the assignment method\n");
		print_usage(argv[0]);
	}

	if (g_cla.umethod == -1) {
		fprintf(stderr, "Please specify the update method\n");
		print_usage(argv[0]);
	}
}

static int parse_opt(Parse_opt opt, const char* method)
{
	switch (opt) {
		case PARSE_ASSIGN_METHOD:
			if (!strcasecmp(method, "classic"))
				g_cla.amethod = ASSIGN_METHOD_CLASSIC;

			else if (!strcasecmp(method, "LSH"))
				g_cla.amethod = ASSIGN_METHOD_RLSH;

			else if (!strcasecmp(method, "Hypercube"))
				g_cla.amethod = ASSIGN_METHOD_RCUBE;

			else if (!strcasecmp(method, "LSH_Frechet") ||
					(!strcasecmp(method, "Frechet")))
				g_cla.amethod = ASSIGN_METHOD_FRECHET;

			else {
				fprintf(stderr, "Unknown method: %s\n", method);
				return 1;
			}
			break;
		case PARSE_UPDATE_METHOD:
			if (!strcasecmp(method, "Mean vector"))
				g_cla.umethod = ASSIGN_METHOD_CLASSIC;

			else if (!strcasecmp(method, "Mean curve"))
				g_cla.umethod = ASSIGN_METHOD_RLSH;

			else {
				fprintf(stderr, "Unknown method: %s\n", method);
				return 1;
			}
			break;
	}

	return 0;
}

static void check_cla(void)
{
	char* alg = NULL;

	if (g_cla.umethod == UPDATE_METHOD_VECTOR &&
	    g_cla.amethod == ASSIGN_METHOD_FRECHET) {
		fprintf(stderr, "You cannot combine \"%s\" assignement with \"%s\"\n",
		        "Frechet", "mean vector");
		exit(EXIT_FAILURE);
	}
	if (g_cla.umethod == UPDATE_METHOD_CURVE) {
	    if (g_cla.amethod == ASSIGN_METHOD_RLSH)
		    alg = "LSH vector";
		else if (g_cla.amethod == ASSIGN_METHOD_RCUBE)
		    alg = "LSH hypercube";

		if (alg) {
			fprintf(stderr, "You cannot combine \"%s\" assignement with \"%s\"\n",
					"vector", alg);
			exit(EXIT_FAILURE);
		}
	}
}

static char* algname(Assign_method alg)
{
	struct {
		Assign_method method;
		char* name;
	} method[] = {
		{ ASSIGN_METHOD_CLASSIC, "Classic Lloyds'" },
		{ ASSIGN_METHOD_RLSH,    "Reverse LSH" },
		{ ASSIGN_METHOD_RCUBE,   "Reverse hypercube" },
		{ ASSIGN_METHOD_FRECHET, "Reverse Frechet" }};

	return method[alg].name;
}


static void read_usr_input(void)
{
	if (!strcmp(g_cla.input_file, "")) {
		printf("\nPlease write the path to the input file\n");
		scanf("%s", g_cla.input_file);
	}

	if (!strcmp(g_cla.config_file, "")) {
		printf("\nPlease write the path to the config file\n");
		scanf("%s", g_cla.config_file);
	}

	if (!strcmp(g_cla.output_file, "")) {
		printf("\nPlease write the path to the output file\n");
		scanf("%s", g_cla.output_file);
	}
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
		g_cla.input_file [0] = '\0';
		g_cla.output_file[0] = '\0';
		g_cla.config_file[0] = '\0';

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
	Clustering* c;

	env_init();
	parse_cla(argc, argv);
	check_cla();

	do {
		read_usr_input();

		printff("\nParsing input files...");
		c = clustering_init(g_cla.amethod, g_cla.umethod, g_cla.input_file,
		                    g_cla.config_file, g_cla.output_file);
		printff("Done\n");

		printff("Initializing algorithm...");
		clustering_init_pp(c);
		printff("Done\n");

		printff("Executing %s algorithm...", algname(g_cla.amethod));
		clustering_exec(c);
		printff(" Done\n");

		if (g_cla.silhouette) {
			printff("\nExecuting silhouette...\n");
			clustering_silhouette(c);
			printff(" Done\n");
		}

		clustering_print(c, g_cla.complete);

		clustering_free(c);

	} while (prompt_rerun());

	printff("Done!\n");
	return 0;
}

Clustering* clustering_init(Assign_method amethod, Update_method umethod,
                            const char* datapath, const char* confpath,
                            const char* outpath)
{
	Clustering* c = xmalloc(sizeof(*c));

	switch (umethod) {
		case UPDATE_METHOD_VECTOR:
			c->metric = euclidean_distance;
			c->dataset = read_dataset(datapath, VECTOR_MODE);
			break;
		case UPDATE_METHOD_CURVE:
			c->metric = frechet_distance;
			c->dataset = read_dataset(datapath, CURVE_MODE);
			break;
	}

	c->amethod = amethod;
	c->umethod = umethod;
	c->params = read_cluster_conf(confpath);
	c->output = fopen(outpath, "wb");
	if (!c->output)
		syserr_exit("Could not open \"%s\"", outpath);

	c->clusters = vector_init(cluster_free_generic);
	c->silval = vector_init(free);
	c->centroids = NULL;
	c->time = 0;

	// Initialize point data
	for (int i = 0; i < c->dataset->size; i++) {
		Math_vector* const point = vector_get(c->dataset, i);
		Point_data* data = xcalloc(1, sizeof(*data));
		math_vector_set_data(point, data);
	}

	return c;
}

void clustering_free(Clustering* c)
{
	// Free point data
	for (int i = 0; i < c->dataset->size; i++) {
		Math_vector* const point = vector_get(c->dataset, i);
		Point_data* const data   = math_vector_data(point);
		free(data);
	}

	if (c->output) fclose(c->output);
	free_dataset(c->dataset);
	vector_free(c->silval);
	vector_free(c->centroids);
	vector_free(c->clusters);
	free(c);
}

static void clustering_create_cluster(Clustering* c, Math_vector* centroid)
{
	vector_push(c->clusters, cluster_init(centroid));
}

static void clustering_clear_points(Clustering* c)
{
	for (int i = 0; i < c->clusters->size; i++)
		cluster_clear(vector_get(c->clusters, i));
}

static void clustering_clear_marks(Clustering* c)
{
	for (int i = 0; i < c->dataset->size; i++) {
		Math_vector* const point = vector_get(c->dataset, i);
		mark(point, NULL);
	}
}

static Cluster*
clustering_second_best_cluster(Clustering* c, const Cluster* cluster,
                               const Math_vector* point)
{
	Cluster* candidate_cluster;
	Cluster* nearest_cluster = NULL;
	double dist;
	double min_dist = DBL_MAX;

	for (int i = 0; i < c->clusters->size; i++) {
		candidate_cluster = vector_get(c->clusters, i);
		if (candidate_cluster != cluster) {
			dist = c->metric(point, candidate_cluster->centroid);
			if (dist < min_dist) {
				min_dist = dist;
				nearest_cluster = candidate_cluster;
			}
		}
	}

	return nearest_cluster;
}

/* Assigns a point to the nearest cluster. Returns the previously assigned
 * cluster. */
static Cluster* clustering_assign_to_nearest(Clustering* c, Math_vector* point)
{
	Cluster* nearest_cluster = NULL;
	Cluster* prev_cluster = NULL;
	Cluster* cluster;
	double min_dist = DBL_MAX;
	double dist;

	for (int i = 0; i < c->clusters->size; i++) {
		cluster = vector_get(c->clusters, i);

		dist = c->metric(point, cluster->centroid);
		if (dist < min_dist) {
			min_dist = dist;
			nearest_cluster = cluster;
		}
	}

	if (nearest_cluster) {
		prev_cluster = cluster_of(point);
		cluster_assign(nearest_cluster, point);
	}

	return prev_cluster;
}

static void clustering_update_centroids(Clustering* c)
{
	for (int i = 0; i < c->clusters->size; i++) {
		Cluster* const cluster = vector_get(c->clusters, i);

		if (c->umethod == UPDATE_METHOD_VECTOR)
			cluster_replace_centroid(cluster, cluster_mean_vector(cluster));

		else if(c->umethod == UPDATE_METHOD_CURVE)
			cluster_replace_centroid(cluster, cluster_mean_curve(cluster));
	}
}

static Vector* clustering_centroids(Clustering* c)
{
	vector_free(c->centroids);
	c->centroids = vector_init(NULL);

	for (int i = 0; i < c->clusters->size; i++) {
		Cluster* const cluster = vector_get(c->clusters, i);
		vector_push(c->centroids, cluster->centroid);
	}

	return c->centroids;
}

/* Calculates the min distance between the centroids of the given clusters */
static double clustering_min_centroid_dist(Clustering* c)
{
	double min_dist = DBL_MAX;
	double dist;
	Cluster* cluster1;
	Cluster* cluster2;

	for (int i = 0; i < c->clusters->size; i++) {
		for (int j = i +1; j < c->clusters->size; j++) {
			cluster1 = vector_get(c->clusters, i);
			cluster2 = vector_get(c->clusters, j);

			dist = c->metric(cluster1->centroid, cluster2->centroid);
			if (dist < min_dist)
				min_dist = dist;
		}
	}

	return min_dist;
}

/* Initializes the clustering process with the k-means++ algorithm. */
void clustering_init_pp(Clustering* c)
{
	Math_vector* centroid;
	Vector* dataset;
	double partial_sum[c->dataset->size];
	double D[c->dataset->size];
	double D_max;
	float x;
	int r;
	int i, j;

	dataset = vector_dup(c->dataset, NULL);

	// Pick a random vector as the centroid of the first cluster
	r = randint(0, dataset->size -1);
	centroid = vector_extract(dataset, r);
	clustering_create_cluster(c, math_vector_dup(centroid));

	// Create K clusters
	for (i = 1; i < c->params.K; ++i) {

		// Calculate the minimum distances and the normalization denominator
		// (max_i D(i))
		D_max = 0;
		for (j = 0; j < dataset->size; ++j) {
			const Math_vector* const point = vector_get(dataset, j);

			D[j] = min_dist(point, clustering_centroids(c), c->metric);
			if (D[j] > D_max)
				D_max = D[j];
		}

		// Calculate the partial sums. Use the normalized distances to avoid
		// overflows
		partial_sum[0] = pow(D[0]/D_max, 2);
		for (j = 1; j < dataset->size; ++j)
			partial_sum[j] = pow(D[j]/D_max, 2) +partial_sum[j -1];

		// Pick a uniformly distributed float x from [0, P(n -t)]
		x = randfloat(0, partial_sum[j -1]);

		// Perform a binary search to locate the corresponding centroids
		r = binary_search(partial_sum, &x, sizeof(double), dataset->size,
		                  centroid_comp, NULL);

		// Create a cluster with the chosen centroid
		centroid = vector_extract(dataset, r);
		clustering_create_cluster(c, math_vector_dup(centroid));
	}

	vector_free(dataset);
}

static int centroid_comp(const void* x_, const void* val_, void* unused)
{
	float x = *(float*) x_;
	double val = *(double*) val_;

	if (x > val)
		return 1;

	if (x < val)
		return -1;

	return 0;
}

void clustering_exec(Clustering* c)
{
	switch (c->amethod) {
		case ASSIGN_METHOD_CLASSIC:
			clustering_lloyds(c);
			break;
		case ASSIGN_METHOD_RLSH:
			clustering_reverse_LSH(c);
			break;
		case ASSIGN_METHOD_RCUBE:
			clustering_reverse_cube(c);
			break;
		case ASSIGN_METHOD_FRECHET:
			clustering_reverse_LSH(c);
			break;
	}
}

static void clustering_lloyds(Clustering* c)
{
	Math_vector* point;
	Cluster* prev_cluster;
	Progbar* bar;
	clock_t time_s, time_e;
	int point_reassigns;
	int iter = 0;
	int i;

	time_s = clock();

	bar = progress_bar_init(60, c->dataset->size);
	printf("\n\nConvergence:\n");
	do {
		point_reassigns = 0;

		// Clear clusters off existing points
		clustering_clear_points(c);

		// Assign points to nearest clusters
		for (i = 0; i < c->dataset->size; ++i) {
			point = vector_get(c->dataset, i);
			prev_cluster = clustering_assign_to_nearest(c, point);

			if (prev_cluster != cluster_of(point))
				point_reassigns++;
		}

		// Update each cluster with the new (mean) centroid
		clustering_update_centroids(c);

		progress_bar_print(bar, c->dataset->size -point_reassigns);

	} while (point_reassigns > CLUSTERING_PRECISION_PCT * c->dataset->size &&
	         iter++ < MAX_ITERATIONS);

	progress_bar_free(bar);

	time_e = clock();

	c->time = (double) (time_e -time_s)/CLOCKS_PER_SEC;
}

static void clustering_reverse_LSH(Clustering* c)
{
	LSH* lsh;
	Progbar* bar;
	clock_t time_s, time_e;
	int point_reassigns = 0;
	int iter = 0;
	size_t htab_size;
	const double delta = 1;

	htab_size = (c->dataset->size/32 >= 8) ? c->dataset->size/8 : 8;

	time_s = clock();

	// Construct an LSH hashtable with the points
	if (c->umethod == UPDATE_METHOD_CURVE)
		lsh = LSH_init(LSH_MODE_CURVE, c->metric, c->dataset, htab_size,
		               c->params.L, c->params.klsh, delta);
	else
		lsh = LSH_init(LSH_MODE_VECTOR, c->metric, c->dataset, htab_size,
		               c->params.L, c->params.klsh);

	bar = progress_bar_init(60, c->dataset->size);
	printf("\n\nConvergence:\n");

	do {
		Cluster* cluster;
		Vector* points;
		double radius;
		int mark_changes = 0;
		int mark_changes_prev = 0;
		int i;

		// Clear clusters off existing points. Clear marks off points
		clustering_clear_points(c);
		clustering_clear_marks(c);

		// Perform a range search and mark the points belonging to each cluster
		// until no more mark changes occur
		for (radius = clustering_min_centroid_dist(c)/2;
			 mark_changes > 0 || (mark_changes -mark_changes_prev) >= 0;
			 radius = 2 * radius)
		{
			mark_changes_prev = mark_changes;
			mark_changes = 0;
			for (i = 0; i < c->clusters->size; ++i) {
				cluster = vector_get(c->clusters, i);
				points  = LSH_range_search(lsh, cluster->centroid, radius);
				mark_changes += mark_points(points, cluster, c->metric);
				vector_free(points);
			}
		}

		// Assign marked points to their corresponing clusters. Assign unmarked
		// points to the nearest clusters
		point_reassigns = 0;
		for (i = 0; i < c->dataset->size; ++i) {
			Math_vector* const point = vector_get(c->dataset, i);
			Cluster* prev_cluster;

			if ((cluster = mark_of(point)) != NULL)
				prev_cluster = cluster_assign(cluster, point);
			else
				prev_cluster = clustering_assign_to_nearest(c, point);

			if (prev_cluster != cluster_of(point))
				point_reassigns++;
		}

		// Update each cluster with the new (mean) centroid
		clustering_update_centroids(c);

		progress_bar_print(bar, c->dataset->size -point_reassigns);

	} while (point_reassigns > CLUSTERING_PRECISION_PCT * c->dataset->size &&
	         iter++ < MAX_ITERATIONS);

	progress_bar_free(bar);
	LSH_free(lsh);

	time_e = clock();

	c->time = (double) (time_e -time_s)/CLOCKS_PER_SEC;
}

static void clustering_reverse_cube(Clustering* c)
{
	Cube* cube;
	Progbar* bar;
	int point_reassigns;
	clock_t time_s, time_e;
	size_t htab_size;
	int iter = 0;

	htab_size = (c->dataset->size/32 >= 8) ? c->dataset->size/8 : 8;

	time_s = clock();

	// Construct an hypercube hashtable with the points
	cube = cube_init(c->dataset, htab_size, c->params.kcube);
	cube_setopt(cube, CUBE_METRIC, c->metric);
	cube_setopt(cube, CUBE_PROBES, c->params.probes);
	cube_setopt(cube, CUBE_M, c->params.M);

	bar = progress_bar_init(60, c->dataset->size);
	printf("\n\nConvergence:\n");

	do {
		Cluster* cluster;
		Vector* points;
		double radius;
		int mark_changes = 0;
		int mark_changes_prev = 0;
		int i;

		// Clear clusters off existing points. Clear marks off points
		clustering_clear_points(c);
		clustering_clear_marks(c);

		// Perform a range search and mark the points belonging to each cluster
		// until no more mark changes occur
		mark_changes = 0;
		mark_changes_prev = 0;
		for (radius = clustering_min_centroid_dist(c)/2;
			 mark_changes > 0 || (mark_changes -mark_changes_prev) >= 0;
			 radius = 2 * radius)
		{
			mark_changes_prev = mark_changes;
			mark_changes = 0;
			for (i = 0; i < c->clusters->size; ++i) {
				cluster = vector_get(c->clusters, i);
				points  = cube_range_search(cube, cluster->centroid, radius);
				mark_changes += mark_points(points, cluster, c->metric);
				vector_free(points);
			}
		}

		// Assign marked points to their corresponing clusters. Assign unmarked
		// points to the nearest clusters
		point_reassigns = 0;

		for (i = 0; i < c->dataset->size; ++i) {
			Math_vector* const point = vector_get(c->dataset, i);
			Cluster* prev_cluster;

			if ((cluster = mark_of(point)) != NULL)
				prev_cluster = cluster_assign(cluster, point);
			else
				prev_cluster = clustering_assign_to_nearest(c, point);

			if (prev_cluster != cluster_of(point))
				point_reassigns++;
		}

		// Update each cluster with the new (mean) centroid
		clustering_update_centroids(c);

		progress_bar_print(bar, c->dataset->size -point_reassigns);

	} while (point_reassigns > CLUSTERING_PRECISION_PCT * c->dataset->size &&
	         iter++ < MAX_ITERATIONS);

	progress_bar_free(bar);
	cube_free(cube);

	time_e = clock();

	c->time = (double) (time_e -time_s)/CLOCKS_PER_SEC;
}

void clustering_silhouette(Clustering* c)
{
	Progbar* bar = progress_bar_init(60, c->dataset->size);
	Cluster* alt_cluster;
	double s_avg[c->clusters->size];
	double s_avg_total = 0;
	double s;
	double a;
	double b;
	double max_ab;
	int point_i = 0;
	int i, j, k;

	for (i = 0; i < c->clusters->size; i++) {
		Cluster* const cluster = vector_get(c->clusters, i);
		s_avg[i] = 0;

		for (j = 0; j < cluster->points->size; j++) {
			Math_vector* const point1 = vector_get(cluster->points, j);

			alt_cluster = clustering_second_best_cluster(c, cluster, point1);

			a = b = 0;
			for (k = 0; k < cluster->points->size; k++) {
				Math_vector* const point2 = vector_get(cluster->points, k);
				a += c->metric(point1, point2);
			}
			a = a/cluster->points->size;

			for (k = 0; k < alt_cluster->points->size; k++) {
				Math_vector* const point2 = vector_get(alt_cluster->points, k);
				b += c->metric(point1, point2);
			}
			b = b/alt_cluster->points->size;

			max_ab = (a > b) ? a : b;

			s = (b -a)/max_ab;
			s_avg[i] += s;

			if (j % 20 == 0)
				progress_bar_print(bar, point_i +j +1);
		}
		point_i += j;

		s_avg[i] = s_avg[i]/cluster->points->size;
		s_avg_total += s_avg[i];
	}
	progress_bar_print(bar, point_i);

	// Store silhouette values for printing later
	double* val;

	for (i = 0; i < c->clusters->size; i++) {
		val = xmalloc(sizeof(*val));
		*val = s_avg[i];
		vector_push(c->silval, val);
	}

	val = xmalloc(sizeof(*val));
	*val = s_avg_total/c->clusters->size;
	vector_push(c->silval, val);

	progress_bar_free(bar);
}

void clustering_print(Clustering* c, bool complete)
{
	int i;

	fprintf(c->output, "Algorithm: %s\n", algname(c->amethod));

	for (i = 0; i < c->clusters->size; i++) {
		Cluster* const cluster = vector_get(c->clusters, i);
		char* linestr = cluster_printstr(cluster, false);
		fprintf(c->output, "CLUSTER-%d %s\n", i +1, linestr);
		free(linestr);
	}

	fprintf(c->output, "clustering_time: %f\n", c->time);

	for (i = 0; i < c->silval->size; i++) {
		double s = *(double*) vector_get(c->silval, i);
		if (i == 0)
			fprintf(c->output, "[%f, ", s);
		else if (i == c->silval->size -1)
			fprintf(c->output, "%f]",  s);
		else
			fprintf(c->output, "%f, ", s);
	}

	if (complete) {
		fprintf(c->output, "\n\n");
		for (i = 0; i < c->clusters->size; i++) {
			Cluster* const cluster = vector_get(c->clusters, i);
			char* linestr = cluster_printstr(cluster, complete);
			fprintf(c->output, "CLUSTER-%d %s\n", i +1, linestr);
			free(linestr);
		}
	}
}

static Cluster* cluster_init(Math_vector* centroid)
{
	Cluster* cluster = xmalloc(sizeof(*cluster));

	cluster->centroid = centroid;
	cluster->points   = vector_init(NULL);

	return cluster;
}

static void cluster_free(Cluster* cluster)
{
	math_vector_free(cluster->centroid);
	vector_free(cluster->points);
	free(cluster);
}

static void cluster_free_generic(void* cluster)
{
	cluster_free(cluster);
}

static Cluster* mark_of(Math_vector* point)
{
	Point_data* data = math_vector_data(point);

	return data->mark;
}

static Cluster* cluster_of(Math_vector* point)
{
	Point_data* data = math_vector_data(point);

	return data->cluster;
}

static void mark(Math_vector* point, Cluster* cluster)
{
	Point_data* data = math_vector_data(point);

	data->mark = cluster;
}

static Cluster* cluster_assign(Cluster* cluster, Math_vector* point)
{
	Point_data* data = math_vector_data(point);
	Cluster* prev_cluster = data->cluster;

	data->cluster = cluster;
	vector_push(cluster->points, point);

	return prev_cluster;
}

static void cluster_replace_centroid(Cluster* cluster, Math_vector* centroid)
{
	if (cluster->centroid)
		math_vector_free(cluster->centroid);

	cluster->centroid = centroid;
}

/* Clears the cluster off its points */
static void cluster_clear(Cluster* cluster)
{
	vector_free(cluster->points);
	cluster->points = vector_init(NULL);
}

static Math_vector* cluster_mean_vector(Cluster* cluster)
{
	Math_vector* sum = math_vector_init(VECTOR_MODE, 0);
	Math_vector* mean;
	Math_vector* point;

	for (int i = 0; i < cluster->points->size; i++) {
		point = vector_get(cluster->points, i);
		vsum(sum, point);
	}

	mean = vmult(sum, (float) 1/cluster->points->size);

	math_vector_free(sum);
	return mean;
}

static Math_vector* cluster_mean_curve(Cluster* cluster)
{
	Vector* tree_upper_lvl;
	Vector* tree_lower_lvl;
	Vector* filtered_curves;
	Math_vector* c;
	Math_vector* d;
	Math_vector* mean;
	Math_vector* filtered;
	int i;

	if (cluster->points->size == 0)
		return math_vector_dup(cluster->centroid);

	tree_lower_lvl = vector_init(NULL);
	for (int i = 0; i < cluster->points->size; i++)
		vector_push(tree_lower_lvl, vector_get(cluster->points, i));

	while (tree_lower_lvl->size > 1) {

		tree_upper_lvl = vector_init(math_vector_free_generic);
		for (i = 0; i < tree_lower_lvl->size; i += 2) {
			c = vector_get(tree_lower_lvl, i);
			d = vector_get(tree_lower_lvl, i +1);
			mean = frechet_mean(c, d);
			vector_push(tree_upper_lvl, mean);
		}

		filtered_curves = vector_init(math_vector_free_generic);
		for (i = 0; i < tree_upper_lvl->size; i++) {
			c = vector_get(tree_upper_lvl, i);
			filtered = curve_filter(c, EPSILON);
			vector_push(filtered_curves, filtered);
		}
		vector_free(tree_lower_lvl);
		vector_free(tree_upper_lvl);

		tree_lower_lvl = filtered_curves;
	}

	mean = math_vector_dup(vector_get(tree_lower_lvl, 0));
	vector_free(tree_lower_lvl);

	return mean;
}

static char* cluster_printstr(Cluster* cluster, bool complete)
{
	char* centstr;
	char* linestr = NULL;
	char* buf;

	centstr = math_vector_printstr(cluster->centroid);

	if (complete == false)
		xsprintf(&linestr, "{size: %zu, centroid: %s}", cluster->points->size,
		         centstr);
	else {
		xsprintf(&linestr, "{");

		for (int i = 0; i < cluster->points->size; i++) {
			Math_vector* point = vector_get(cluster->points, i);

			if (i == cluster->points->size -1)
				xsprintf(&buf, "%s", math_vector_id(point));
			else
				xsprintf(&buf, "%s, ", math_vector_id(point));

			xstrcat(&linestr, buf);
			free(buf);
		}

		xstrcat(&linestr, "}");
	}

	free(centstr);

	return linestr;
}

/* Calculates the min distance between the given vector and the given set of
 * clusters */
static double min_dist(const Math_vector* v, Vector* points,
                       Clustering_metric metric)
{
	double mindist = DBL_MAX;
	double dist;

	for (int i = 0; i < points->size; i++) {
		dist = metric(v, vector_get(points, i));
		if (dist < mindist)
			mindist = dist;
	}

	return mindist;
}

static int mark_points(const Vector* points, Cluster* cluster,
                       Clustering_metric metric)
{
	int mark_changes = 0;

	for (int i = 0; i < points->size; i++) {
		Math_vector* const point = vector_get(points, i);
		Point_data* const data   = math_vector_data(point);

		if (data->mark == NULL || (cluster != data->mark &&
		    metric(point, cluster->centroid) <
		    metric(point, data->mark->centroid)))
		{
			mark(point, cluster);
			mark_changes++;
		}
	}

	return mark_changes;
}
