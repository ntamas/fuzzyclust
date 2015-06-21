/* vim: set ts=4 sw=4 sts=4 et */
/*
    Fuzzy community detection in complex networks.
    (C) 2007-2011 Tamas Nepusz <nepusz@hal.elte.hu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <string.h>

#define VERSION "0.11"

#define sqr(x) ((x) * (x))

typedef igraph_matrix_t membership_matrix_t;

#define OPTIMIZE_NO 0
#define OPTIMIZE_YES 1
#define OPTIMIZE_AUTO 2

typedef enum { SIMILARITY_ADJACENCY=0, SIMILARITY_JACCARD } similarity_type_t;

typedef struct {
    similarity_type_t similarity_type;
    igraph_bool_t weights_from_modularity;
    int no_of_clusters;
    igraph_bool_t adaptive_cluster_count;
    int verbosity;
    igraph_real_t dominance_threshold;
    igraph_real_t prune_threshold;
    long int prune_frequency;
    char* filename;
    char* weight_file;
	unsigned int seed;
	igraph_bool_t has_seed;
} parameters_t;

parameters_t params;              /* Global array for holding command line args */

#ifndef __MSVC__
#  define INFO(...) { if (params.verbosity>=1) fprintf(stderr, __VA_ARGS__); }
#  define WARNING(...) fprintf(stderr, __VA_ARGS__)
#  define FATAL(...) { fprintf(stderr, __VA_ARGS__); exit(1); }
#else
inline void INFO(const char* format, ...) {
    if (params.verbosity<1) return;
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}
inline void WARNING(const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}
inline void FATAL(const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(1);
}
#endif

typedef igraph_real_t goal_function_t(void*);
typedef int iteration_function_t(void*, igraph_bool_t);

typedef struct {
    long int src, dst;            /* Source, destination of an edge */
    igraph_real_t similarity;     /* Similarity of the two vertices */
    igraph_real_t weight;         /* Weight of this specification */
    igraph_real_t target;         /* Target similarity we want to achieve */
} edge_t;

typedef struct {
    long int n, k;                /* Size of membership matrix */
	igraph_bool_t is_directed;    /* Is it a directed state? */
    membership_matrix_t u;        /* Current membership matrix */
    membership_matrix_t v;        /* Additional membership matrix if directed */
    igraph_real_t *self_sims;     /* Self-similarity values */
    igraph_real_t goal;           /* Value of the goal function (if calculated) */
} state_t;

typedef struct {
    const igraph_t *graph;        /* Graph being clustered */
	igraph_bool_t is_directed;    /* Is it a directed fuzzy clustering? */
    edge_t *edges;                /* Edge similarity spec. list of the graph */
    long int edge_count;          /* Length of the edge similarity specifications */
    igraph_vector_t outdegrees;   /* Degree list of the graph */
    igraph_vector_t indegrees;    /* Degree list of the graph */
    long int k;                   /* Current number of clusters */
    state_t current;              /* Current state */
    state_t local_best;           /* Local best state, for adaptive cluster count */
    state_t best;                 /* Best state */
    long int successes;           /* Number of successes in a row */
    igraph_real_t step_size;      /* Step size */
    igraph_real_t dominance_threshold;
    igraph_real_t epsilon;
    goal_function_t *goal_function;
    iteration_function_t *iterate;
} fuzzy_clustering_t;

// Auxiliary function to create a lexicographically sorted edge list
int edge_sort_aux(const void* a, const void* b) {
    edge_t *e1 = (edge_t*)a;
    edge_t *e2 = (edge_t*)b;
    if (e1->src == e2->src) {
        return e1->dst - e2->dst;
    }
    return e1->src - e2->src;
}

// Sampling from a gamma distribution
igraph_real_t sample_gamma(long int alpha, igraph_real_t beta) {
    /* Algorithm taken from Wikipedia. Works only for integer alphas */
    long int i;
    igraph_real_t result, r;

    result=0.0;
    for (i=0; i<alpha; i++) {
        r=0.0;
        while (r<0.001) { r = (float)rand() / RAND_MAX; }
        result -= log(r);
    }
    
    return result/beta;
}

/******************* Membership matrix functions *******************/
void membership_matrix_normalize(membership_matrix_t *c);

// Initializes a membership matrix
int membership_matrix_init(membership_matrix_t *c,
        long int n, long int k) {
    long int i, j;
    IGRAPH_CHECK(igraph_matrix_init(c, n, k));
    /* Fill it with random numbers between 0 and 1 */
    for (i=0; i<n; i++)
        for (j=0; j<k; j++)
            MATRIX(*c,i,j) = sample_gamma(1, 1);

    membership_matrix_normalize(c);
    return 0;
}

// Normalizes the columns of a membership matrix
void membership_matrix_normalize(membership_matrix_t *c) {
    long int i, j, n, k;
    igraph_real_t sum;
    n = igraph_matrix_nrow(c);
    k = igraph_matrix_ncol(c);
    for (i=0; i<n; i++) {
        for (j=0, sum=0; j<k; j++) sum += MATRIX(*c, i, j);
        for (j=0; j<k; j++) MATRIX(*c, i, j) /= sum;
    }
}

// Hadamard power of a membership matrix -- not used at the moment
void membership_matrix_hadamard_power(membership_matrix_t *c,
    igraph_real_t power) {
    long int i, j, n, k;
    igraph_real_t sum;
    n = igraph_matrix_nrow(c);
    k = igraph_matrix_ncol(c);
    for (i=0; i<n; i++) {
        for (j=0, sum=0; j<k; j++) {
            MATRIX(*c, i, j) = pow(MATRIX(*c, i, j), power);
            sum += MATRIX(*c, i, j);
        }
        for (j=0; j<k; j++) MATRIX(*c, i, j) /= sum;
    }
}

// Prunes a membership matrix by setting elements below a given epsilon to zero
void membership_matrix_prune(membership_matrix_t *c, igraph_real_t eps) {
    long int i, j, n, k;
    igraph_bool_t changed;
    igraph_real_t sum;

    n = igraph_matrix_nrow(c);
    k = igraph_matrix_ncol(c);
    for (i=0; i<n; i++) {
        for (j=0, sum=0, changed=0; j<k; j++) {
            if (MATRIX(*c, i, j) < eps) {
                MATRIX(*c, i, j) = 0.0;
                changed = 1;
            } else sum += MATRIX(*c, i, j);
        }
        if (changed)
            for (j=0; j<k; j++) MATRIX(*c, i, j) /= sum;
    }
}

// Returns the similarity of two vertices using a membership matrix (two if directed)
igraph_real_t membership_matrix_get_similarity(const membership_matrix_t *u,
        const membership_matrix_t *v, long int v1, long int v2) {
    long int i, k;
    igraph_real_t result = 0.0;
    k = igraph_matrix_ncol(u);
    for (i=0; i<k; i++) result += MATRIX(*u, v1, i)*MATRIX(*v, v2, i);
    return result;
}

// Adds a new community (or more, depending on dk) to the membership matrix
int membership_matrix_add_community(membership_matrix_t *c, long int dk) {
    /* WARNING: this function directly manipulates the internal representation
     * of igraph_matrix_t. Changes in the internals of igraph_matrix_t in the
     * igraph library might break this function. This function should work with
     * igraph 0.5 and 0.6 */
    long int n, k;
    k = igraph_matrix_ncol(c);
    n = igraph_matrix_nrow(c);
    IGRAPH_CHECK(igraph_vector_resize(&c->data, n * (k+dk)));
    /* Due to the current implementation of igraph_vector_resize and the
     * fact that igraph matrices are stored in row major format, we must fill the
     * new elements of the vector with zeroes and we're ready */
    memset(&VECTOR(c->data)[n*k], 0, sizeof(igraph_real_t) * n * dk);
    c->ncol += dk;
    return 0;
}

// Retrieves the size of each community
int membership_matrix_get_community_sizes(const membership_matrix_t *m,
		igraph_vector_t *vector) {
	long int i, j, n, k;
	
    n = igraph_matrix_nrow(m);
    k = igraph_matrix_ncol(m);
	IGRAPH_CHECK(igraph_vector_resize(vector, k));

	for (i=0; i < k; i++) {
		double result = 0.0;
		for (j=0; j < n; j++) {
			result += MATRIX(*m, j, i);
		}
		VECTOR(*vector)[i] = result;
	}

	return 0;
}

// Removes empty community from the membership matrix
int membership_matrix_remove_communities(membership_matrix_t *c, int* remapping) {
    /* WARNING: this function directly manipulates the internal representation
     * of igraph_matrix_t. Changes in the internals of igraph_matrix_t in the
     * igraph library might break this function. This function should work with
     * igraph 0.5 and 0.6 */
    long int i, j, n, k;
	membership_matrix_t new_c;

    k = igraph_matrix_ncol(c);
    n = igraph_matrix_nrow(c);

	j = -1;
	for (i = 0; i < k; i++)
		if (remapping[i] > j)
			j = remapping[i];
	j += 1;

	IGRAPH_CHECK(igraph_matrix_init(&new_c, n, j));
	IGRAPH_FINALLY(igraph_matrix_destroy, &new_c);

	for (j = 0; j < k; j++) {
		int new_j = remapping[j];
		if (new_j == -1) continue;

		for (i = 0; i < n; i++) {
			MATRIX(new_c, i, new_j) = MATRIX(*c, i, j);
		}
	}

	membership_matrix_normalize(&new_c);
	igraph_matrix_update(c, &new_c);
	igraph_matrix_destroy(&new_c);
	IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

// Destroys a membership matrix
void membership_matrix_destroy(membership_matrix_t *c) {
    igraph_matrix_destroy(c);
}

/****************** State objects ************************/
void state_recalculate(state_t *s);

// Initializes a state object for the algorithm
int state_init(state_t *s, long int n, long int k, igraph_bool_t is_directed) {
    s->n = n;
    s->k = k;
	s->is_directed = is_directed;

    /* Create membership matrix */
    IGRAPH_CHECK(membership_matrix_init(&s->u, n, k));
    IGRAPH_FINALLY(membership_matrix_destroy, &s->u);
    if (is_directed) {
		IGRAPH_CHECK(membership_matrix_init(&s->v, n, k));
		IGRAPH_FINALLY(membership_matrix_destroy, &s->v);
	}

    /* Create storage space for self similarity values */
    s->self_sims = (igraph_real_t*)calloc(n, sizeof(igraph_real_t));
    if (s->self_sims == 0) {
        IGRAPH_ERROR("can't create state object", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(free, s->self_sims);

    state_recalculate(s);

    s->goal = -1; /* -1 = not calculated yet */

	if (is_directed) IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

// Recalculates the derived variables of the state object
void state_recalculate(state_t *s) {
    long int i, j;
   
	if (!s->is_directed) {
		/* Recalculate self-similarities */
		for (i=0; i < s->n; i++) {
			s->self_sims[i] = 0.0;
			for (j=0; j < s->k; j++)
				s->self_sims[i] += sqr(MATRIX(s->u, i, j));
		}
	} else {
		/* Recalculate self-similarities */
		for (i=0; i < s->n; i++) {
			s->self_sims[i] = 0.0;
			for (j=0; j < s->k; j++)
				s->self_sims[i] += MATRIX(s->u, i, j)*MATRIX(s->v, i, j);
		}
	}
}

// Copies a state object to another
void state_copy(state_t *dest, const state_t *src) {
    dest->n = src->n;
    dest->k = src->k;
    dest->goal = src->goal;
	dest->is_directed = src->is_directed;
    igraph_matrix_copy(&dest->u, &src->u);
	if (src->is_directed) {
		igraph_matrix_copy(&dest->v, &src->v);
	}
    memcpy(dest->self_sims, src->self_sims, sizeof(igraph_real_t)*dest->n);
}

// Adds a new community to the state. Used when adaptively determining the number of communities.
int state_add_community(state_t *s, long int n) {
    s->k += n;
    
    IGRAPH_CHECK(membership_matrix_add_community(&s->u, n));
	if (s->is_directed) {
		IGRAPH_CHECK(membership_matrix_add_community(&s->v, n));
	}

    return 0;
}

// Removes the empty communities from the state.
int state_remove_communities(state_t *s, int *remapping) {
	long int n = 1;

	s->k -= n;
	IGRAPH_CHECK(membership_matrix_remove_communities(&s->u, remapping));

	if (s->is_directed) {
		IGRAPH_CHECK(membership_matrix_remove_communities(&s->v, remapping));
	}

	return 0;
}

// Retrieves the sizes of the communities to an initialized igraph vector
int state_get_community_sizes(const state_t *s, igraph_vector_t *vector) {
	return membership_matrix_get_community_sizes(&s->u, vector);
}

// Takes the Hadamard power of the membership matrix in the current state. Not used ATM.
void state_hadamard_power(state_t *s, igraph_real_t power) {
    membership_matrix_hadamard_power(&s->u, power);
	if (s->is_directed) membership_matrix_hadamard_power(&s->v, power);
    state_recalculate(s);
}

// Prunes the membership matrix in the current state.
void state_prune(state_t *s, igraph_real_t eps) {
    membership_matrix_prune(&s->u, eps);
    if (s->is_directed) membership_matrix_prune(&s->v, eps);
    state_recalculate(s);
}

// Destroys the state object
void state_destroy(state_t *s) {
    if (s->is_directed) {
		igraph_matrix_destroy(&s->v);
	}
	igraph_matrix_destroy(&s->u);
    free(s->self_sims);
}

/****************** Fuzzy clustering structure functions ****************/
void fuzzy_clustering_recalculate(fuzzy_clustering_t *f);
igraph_real_t fuzzy_clustering_simple_goal_function(const fuzzy_clustering_t *f);
int fuzzy_clustering_simple_step(fuzzy_clustering_t *f, igraph_bool_t prune_now);

// Initializes the algorithm with a given graph, number of clusters and dominance threshold
int fuzzy_clustering_init(fuzzy_clustering_t *f, const igraph_t *g, int numcl,
    igraph_real_t dom_thr) {
    long int i, j, k, m, n;
    igraph_integer_t from, to;
    igraph_eit_t eit;

    f->graph = g;
    f->k = numcl;
    f->step_size = 0.25;
    f->dominance_threshold = dom_thr;
    f->epsilon = 1e-2;
    f->successes = 0;
	f->is_directed = igraph_is_directed(g);

    IGRAPH_CHECK(state_init(&f->current,igraph_vcount(g),f->k,f->is_directed));
    IGRAPH_FINALLY(state_destroy, &f->current);

    IGRAPH_CHECK(state_init(&f->best,igraph_vcount(g),f->k,f->is_directed));
    IGRAPH_FINALLY(state_destroy, &f->best);

    IGRAPH_CHECK(state_init(&f->local_best,igraph_vcount(g),f->k,f->is_directed));
    IGRAPH_FINALLY(state_destroy, &f->local_best);

    IGRAPH_VECTOR_INIT_FINALLY(&f->indegrees, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&f->outdegrees, 0);

    /* Store degrees */
    igraph_degree(g, &f->outdegrees, igraph_vss_all(), IGRAPH_OUT, 0);
	igraph_degree(g, &f->indegrees, igraph_vss_all(), IGRAPH_IN, 0);
    /* Calculate edge list */
    if (params.weight_file) {
        /* pre-specified weights */
        FILE *wf = fopen(params.weight_file, "r");
        char buf[4096];
        double w, maxw;
        long int lines=0;
        igraph_bool_t b=0;

        if (!wf) {
            perror(params.weight_file);
            abort();
        }

        m = 100+igraph_vcount(f->graph); n = igraph_vcount(f->graph);
        f->edges = (edge_t*)calloc(m, sizeof(edge_t));
        if (f->edges == 0) {
            IGRAPH_ERROR("can't create fuzzy clustering object", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(free, f->edges);
        for (i=0; i<igraph_vcount(f->graph); i++) {
            f->edges[i].src=i;
            f->edges[i].dst=i;
            f->edges[i].weight=1;
            f->edges[i].target=1;
        }

        INFO("Reading weight file %s...\n", params.weight_file);
		maxw=0;
        while (fgets(buf, 4095, wf)) {
            lines++;
            if (buf[0] == '\r' || buf[0] == '\n' ||
                buf[0] == '#' || buf[0] == '%') continue; 
            k=sscanf(buf, "%ld %ld %lf", &i, &j, &w);
            if (k < 2) {
                WARNING("Malformed input skipped at %s:%ld\n", params.weight_file, lines);
            } else if (k > 2 && w < 0) {
				WARNING("Negative weight skipped at %s:%ld\n", params.weight_file, lines);
			} else {
                if (i == j) continue;
                if (k == 2) w = 1.0;
				if (w > maxw) maxw = w;
                igraph_are_connected(f->graph, i, j, &b);
                f->edges[n].src=i;
                f->edges[n].dst=j;
                f->edges[n].weight=(igraph_real_t)w;
                f->edges[n].target=b?1:0;
                n++;
                if (n == m) {
                    m *= 2;
                    f->edges = (edge_t*)realloc(f->edges, sizeof(edge_t)*m);
                    if (f->edges == 0) {
                        IGRAPH_ERROR("can't create fuzzy clustering object", IGRAPH_ENOMEM);
                    }
                }
				if (!f->is_directed) {
					f->edges[n].src=j;
					f->edges[n].dst=i;
					f->edges[n].weight=(igraph_real_t)w;
					f->edges[n].target=b?1:0;
					n++;
					if (n == m) {
						m *= 2;
						f->edges = (edge_t*)realloc(f->edges, sizeof(edge_t)*m);
						if (f->edges == 0) {
							IGRAPH_ERROR("can't create fuzzy clustering object", IGRAPH_ENOMEM);
						}
					}
				}
            }
        }
        fclose(wf);
    } else {
        /* Non-optimized case, there will be an entry for every vertex pair */
        n = igraph_vcount(f->graph);
        f->edges = (edge_t*)calloc(n*n, sizeof(edge_t));
        if (f->edges == 0) {
            IGRAPH_ERROR("can't create fuzzy clustering object", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(free, f->edges);

        k = 0;
        m = igraph_ecount(g);
		if (!f->is_directed) m *= 2;
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++, k++) {
                f->edges[k].src=i;
                f->edges[k].dst=j;
				f->edges[k].target=(i==j)?1:0;
                f->edges[k].weight=1;
            }
        }

        if (params.similarity_type == SIMILARITY_ADJACENCY) {
            IGRAPH_CHECK(igraph_eit_create(g, igraph_ess_all(IGRAPH_EDGEORDER_FROM), &eit));
            IGRAPH_FINALLY(igraph_eit_destroy, &eit);
            while (!IGRAPH_EIT_END(eit)) {
                i = IGRAPH_EIT_GET(eit);
                IGRAPH_CHECK(igraph_edge(g, i, &from, &to));
                IGRAPH_EIT_NEXT(eit);
                if (from==to) continue;
                f->edges[(long)(from*n+to)].target += 1;
                if (!f->is_directed) f->edges[(long)(to*n+from)].target += 1;
            }
            igraph_eit_destroy(&eit);
            IGRAPH_FINALLY_CLEAN(1);
        } else if (params.similarity_type == SIMILARITY_JACCARD) {
            igraph_matrix_t similarities;
            IGRAPH_CHECK(igraph_matrix_init(&similarities, 0, 0));
            IGRAPH_FINALLY(igraph_matrix_destroy, &similarities);
            //igraph_similarity_jaccard(f->graph, &similarities, igraph_vss_all(), IGRAPH_ALL, 1);
            k=0;
            for (i=0; i<n; i++)
              for (j=0; j<n; j++, k++)
                f->edges[k].target = MATRIX(similarities, i, j);
            igraph_matrix_destroy(&similarities);
            IGRAPH_FINALLY_CLEAN(1);
        }

        n = igraph_vcount(f->graph) * igraph_vcount(f->graph);
    }
    f->edge_count = n;
    qsort(f->edges, f->edge_count, sizeof(edge_t), edge_sort_aux);
    INFO("Created %ld constraints\n", f->edge_count);

    /* Add modularity weights if necessary */
    if (params.weights_from_modularity) {
        m = igraph_ecount(g);
		if (!f->is_directed) m *= 2;
        for (k=0; k < f->edge_count; k++) {
            f->edges[k].weight = sqr(f->edges[k].target - (VECTOR(f->outdegrees)[f->edges[k].src] * VECTOR(f->indegrees)[f->edges[k].dst]) / ((float)m));
        }
    }

    /* Initialize goal function and iteration function pointer */
    f->goal_function = (goal_function_t*)fuzzy_clustering_simple_goal_function;
    f->iterate = (iteration_function_t*)fuzzy_clustering_simple_step;

    /* Calculate initial similarities */
    fuzzy_clustering_recalculate(f);
    /* Calculate initial goal function */
    f->current.goal = f->goal_function(f);

    /* Copy current state to the best */
    state_copy(&f->best, &f->current);
    state_copy(&f->local_best, &f->current);

    IGRAPH_FINALLY_CLEAN(5);
    return 0;
}

void fuzzy_clustering_recalculate(fuzzy_clustering_t *f) {
    long int i;
	membership_matrix_t *u, *v;

	u = &f->current.u;
	v = f->is_directed ? &f->current.v : &f->current.u;

    /* Recalculate similarities for edges */
    for (i=0; i < f->edge_count; i++) {
        f->edges[i].similarity =
            membership_matrix_get_similarity(u, v, f->edges[i].src, f->edges[i].dst);
    }
   
    /* Recalculate state variables */
    state_recalculate(&f->current);
}

/* PRE: the similarities have already been calculated */
igraph_real_t fuzzy_clustering_simple_goal_function(const fuzzy_clustering_t *f) {
    igraph_real_t result = 0.0, x;
    long int i, n;

    /* Loop through the edge list and calculate the differences from their
     * target similarities, square them and add them together */
    n = f->edge_count;
    for (i=0; i<n; i++) {
        x = f->edges[i].weight*sqr(f->edges[i].target-f->edges[i].similarity);
        result += x;
    }

    return result;
}

igraph_real_t fuzzy_clustering_modularity(const fuzzy_clustering_t *f) {
    igraph_real_t result = 0.0;
	edge_t *e;
    long int i, k, n, m;
	const membership_matrix_t *u, *v;

    m = igraph_ecount(f->graph);
	if (!f->is_directed) m *= 2;
    n = igraph_vcount(f->graph);
    k = 0;
	u = &f->current.u;
	v = f->is_directed ? &(f->current.v) : &(f->current.u);
    n = f->edge_count;
    for (i=0, e=f->edges; i<n; i++, e++) {
        if (e->src == e->dst) continue;
        result += (e->target - VECTOR(f->outdegrees)[e->src] * VECTOR(f->indegrees)[e->dst] / m) * e->similarity;
    }
    result /= m;
    return result; 
}

void fuzzy_clustering_hadamard_power(fuzzy_clustering_t *f, igraph_real_t power) {
    state_hadamard_power(&f->current, power);
    fuzzy_clustering_recalculate(f);
}

void fuzzy_clustering_prune(fuzzy_clustering_t *f, igraph_real_t eps) {
    state_prune(&f->current, eps);
    fuzzy_clustering_recalculate(f);
    f->current.goal = f->goal_function(f);
}

int fuzzy_clustering_simple_step(fuzzy_clustering_t *f, igraph_bool_t prune_now) {
    /* Calculating gradients */

    /* First we must calculate the product of the difference matrix and the
     * membership matrix. The difference matrix (e) itself is the difference
     * of the similarity matrix and the adjacency matrix, multiplied by the
	 * weights.
     *
     * We do the multiplication by looping through the edge list. Since the
     * edge list is sorted by the source vertex ID, it is enough to loop through
     * the list only once.
     */
    igraph_matrix_t grads_u, grads_v, e;
    long int i, j, k, l, n, eidx;
    igraph_real_t nu, g, glen, maxglen;
	membership_matrix_t *u, *v;

	u = &f->current.u;
	v = f->is_directed ? &f->current.v : &f->current.u;

	n = igraph_vcount(f->graph);
    IGRAPH_CHECK(igraph_matrix_init(&e, n, n));
    IGRAPH_FINALLY(igraph_matrix_destroy, &e);
    IGRAPH_CHECK(igraph_matrix_init(&grads_u, n, f->k));
    IGRAPH_FINALLY(igraph_matrix_destroy, &grads_u);
    IGRAPH_CHECK(igraph_matrix_init(&grads_v, n, f->k));
    IGRAPH_FINALLY(igraph_matrix_destroy, &grads_v);

    eidx = 0;
    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {
            if (f->edges[eidx].src == i && f->edges[eidx].dst == j) {
                /* We have an edge specification between i and j */
				MATRIX(e, i, j) = f->edges[eidx].weight * (f->edges[eidx].target - f->edges[eidx].similarity);
				eidx++;
			}
		}
	}

    /*for (i=0; i < n; i++) 
        for (j=0; j < n; j++) 
			for (l=0; l < f->k; l++)
				MATRIX(dc, i, l) -= MATRIX(e, i, j) * MATRIX(*u, j, l);*/

    /*maxglen = 0.0;
    for (i=0; i < n; i++) {
        nu = 0.0;
        for (l=0, k=0; l < f->k; l++) {
            if (MATRIX(f->current.u, i, l) > 0) {
                k++;
                nu += MATRIX(dc, i, l);
            }
        }
        nu += (f->current.self_sims[i]-1);
        nu *= 2.0/k;
        
        memset(mu, 0, sizeof(igraph_real_t)*f->k);
        for (l=0; l < f->k; l++) 
            if (MATRIX(f->current.u, i, l) == 0)
                mu[l] = 2*MATRIX(dc, i, l)-nu;

        glen = 0.0;
        for (l=0; l < f->k; l++) {
            g = 2 * MATRIX(dc, i, l) - mu[l] - nu;
            g += 2 * (f->current.self_sims[i]-1) * MATRIX(f->current.u, i, l);
            if (g < 1e-6 && g > -1e-6) g=0;
            MATRIX(grads, i, l) = -g;
            glen += sqr(g);
        }
        if (glen > maxglen) maxglen = glen;
    }*/
	maxglen = 0.0;
	for (l=0; l < n; l++) {
		glen = 0.0;
		for (k=0; k < f->k; k++) {
			g = 0.0;
			for (i=0; i < n; i++) 
				g -= MATRIX(e, l, i)*(1.0/f->k - MATRIX(*v, i, k));
			MATRIX(grads_u, l, k) = g;
			glen += sqr(g);
		}
		if (glen > maxglen) maxglen = glen;
	}
	if (f->is_directed) {
		for (l=0; l < n; l++) {
			glen = 0.0;
			for (k=0; k < f->k; k++) {
				g = 0.0;
				for (i=0; i < n; i++) 
					g -= MATRIX(e, i, l)*(1.0/f->k - MATRIX(*u, i, k));
				MATRIX(grads_v, l, k) = g;
				glen += sqr(g);
			}
			if (glen > maxglen) maxglen = glen;
		}
	}

    /* printf("Gradients:\n");
    for (i=0; i < n; i++) {
        printf("%4ld -- ", i);
        for (l=0; l < f->k; l++) printf("%.4f ", MATRIX(grads, i, l));
        printf("\n");
    }
    printf("%.4f\n", sqrt(maxglen));
    printf("%.4f\n", f->epsilon); */

    if (prune_now) {
        /* Now, having calculated the gradients, check the current cluster
         * profile matrix. If a vertex is almost excluded from a cluster
         * (its membership is less than the pruning threshold) and its
         * corresponding gradient is negative (meaning that it really should
         * be excluded), set it to zero. */
        j=0;
        for (i=0; i < n; i++) {
          for (k=0; k < f->k; k++) {
            if (MATRIX(f->current.u, i, k) < params.prune_threshold &&
              MATRIX(f->current.u, i, k) > 0 && MATRIX(grads_u, i, k) < 0) {
              MATRIX(f->current.u, i, k) = 0.0;
              j++;
            }
          }
        }
		if (f->is_directed) {
			for (i=0; i < n; i++) {
			  for (k=0; k < f->k; k++) {
				if (MATRIX(f->current.v, i, k) < params.prune_threshold &&
				  MATRIX(f->current.v, i, k) > 0 && MATRIX(grads_v, i, k) < 0) {
				  MATRIX(f->current.v, i, k) = 0.0;
				  j++;
				}
			  }
			}
		}
        if (j>0) INFO("Pruned %ld elements of the membership matrix\n", j);
    }

    if (maxglen > f->epsilon) {
        /* There is a gradient longer than eps, so do a step */
        nu = f->epsilon;
        for (i=0; i < n; i++) {
            k=0;
            for (l=1; l < f->k; l++)
                if (MATRIX(grads_u, i, l) > MATRIX(grads_u, i, k)) k=l;
            if (MATRIX(grads_u, i, k) > nu) 
                MATRIX(f->current.u, i, k) += MATRIX(grads_u, i, k) * f->step_size;
        }
        membership_matrix_normalize(&f->current.u);
		if (f->is_directed) {
			for (i=0; i < n; i++) {
				k=0;
				for (l=1; l < f->k; l++)
					if (MATRIX(grads_v, i, l) > MATRIX(grads_v, i, k)) k=l;
				if (MATRIX(grads_v, i, k) > nu) 
					MATRIX(f->current.v, i, k) += MATRIX(grads_v, i, k) * f->step_size;
			}
			membership_matrix_normalize(&f->current.v);
		}

        /* OK, calculated the gradients, stepped towards the best component
         * at every vertex. Now evaluate our new position */
        fuzzy_clustering_recalculate(f);
        f->current.goal = f->goal_function(f);

        IGRAPH_FINALLY_CLEAN(3);
		igraph_matrix_destroy(&e);
        igraph_matrix_destroy(&grads_u);
        igraph_matrix_destroy(&grads_v);

        return 1;
    }

    IGRAPH_FINALLY_CLEAN(3);
    igraph_matrix_destroy(&e);
    igraph_matrix_destroy(&grads_u);
    igraph_matrix_destroy(&grads_v);

    return 0;
}

void fuzzy_clustering_mutate(fuzzy_clustering_t *f, igraph_real_t amount) {
    long int i, j, n;

    n = igraph_vcount(f->graph);

    for (i=0; i < n; i++) 
        for (j=0; j < f->k; j++) 
            MATRIX(f->current.u, i, j) += ((float)(rand()) / RAND_MAX) * amount;
    membership_matrix_normalize(&f->current.u);
    if (f->is_directed) {
		for (i=0; i < n; i++) 
			for (j=0; j < f->k; j++) 
				MATRIX(f->current.v, i, j) += ((float)(rand()) / RAND_MAX) * amount;
		membership_matrix_normalize(&f->current.v);
	}
    fuzzy_clustering_recalculate(f);
    f->current.goal = f->goal_function(f);
}

int fuzzy_clustering_add_community(fuzzy_clustering_t *f, long int n) {
    f->k += n;
    IGRAPH_CHECK(state_add_community(&f->current, n));
    IGRAPH_CHECK(state_add_community(&f->best, n));
    IGRAPH_CHECK(state_add_community(&f->local_best, n));
    return 0;
}

int fuzzy_clustering_remove_empty_communities(fuzzy_clustering_t *f, double eps) {
	int i, j, n, k;
	igraph_vector_t community_sizes;
	int* remapping;

	k = f->k;
	n = f->current.n;

	IGRAPH_VECTOR_INIT_FINALLY(&community_sizes, f->k);
	IGRAPH_CHECK(membership_matrix_get_community_sizes(&f->current.u, &community_sizes));

	remapping = (int*)calloc(k, sizeof(int));
	IGRAPH_FINALLY(free, remapping);
	j = 0;
	for (i = 0; i < k; i++) {
		if (VECTOR(community_sizes)[i] < eps)
			remapping[i] = -1;
		else {
			remapping[i] = j;
			j++;
		}
	}
	f->k = j;
	igraph_vector_destroy(&community_sizes);
	IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(state_remove_communities(&f->current, remapping));
    IGRAPH_CHECK(state_remove_communities(&f->best, remapping));
    IGRAPH_CHECK(state_remove_communities(&f->local_best, remapping));

	free(remapping);
	IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

int fuzzy_clustering_print(const fuzzy_clustering_t *f, FILE* out) {
    long int i, j, k=0, l=0, n;
    /* static const char* dots = " .oO#"; */
    igraph_real_t *bs_u = 0, *bs_v = 0, *corr = 0, *diff = 0;
    igraph_real_t d_mean, d_std, kinv;
    igraph_real_t delta_u, delta_v, s, current;
    n = igraph_vcount(f->graph);

    bs_u = (igraph_real_t*)calloc(n, sizeof(igraph_real_t));
    if (bs_u == 0) IGRAPH_ERROR("can't print fuzzy clustering", IGRAPH_ENOMEM);
    bs_v = (igraph_real_t*)calloc(n, sizeof(igraph_real_t));
    if (bs_v == 0) IGRAPH_ERROR("can't print fuzzy clustering", IGRAPH_ENOMEM);
    corr = (igraph_real_t*)calloc(n, sizeof(igraph_real_t));
    if (corr == 0) IGRAPH_ERROR("can't print fuzzy clustering", IGRAPH_ENOMEM);
    diff = (igraph_real_t*)calloc(n, sizeof(igraph_real_t));
    if (diff == 0) IGRAPH_ERROR("can't print fuzzy clustering", IGRAPH_ENOMEM);

    fprintf(out, "%% Input file: %s\n", params.filename);
    if (params.weight_file)
        fprintf(out, "%% Weight file: %s\n", params.weight_file);
    fprintf(out, "%% Fuzzy clustering for %ld clusters\n", (long)f->k);
    if (f->dominance_threshold < 0)
        fprintf(out, "%% No dominance threshold\n");
    else
        fprintf(out, "%% Dominance threshold = %.4f\n", (float)f->dominance_threshold);

    /* Calculate mean degree and std */
    d_mean = 0.0; s = 0.0;
    for (i=0; i < n; i++) {
        delta_u = VECTOR(f->outdegrees)[i] - d_mean;
        d_mean += delta_u / (i+1);
        s += delta_u * (VECTOR(f->outdegrees)[i] - d_mean);
    }
    d_std = sqrt(s/(n-1));

	kinv = 1.0/f->k;
    /* Calculate bridgenesses and correlation */
    for (i=0; i < n; i++) {
        for (j=0; j < f->k; j++) {
            current = MATRIX(f->current.u, i, j);
            bs_u[i] += sqr(current-kinv);
        }
		if (f->is_directed) {
			for (j=0; j < f->k; j++) {
				current = MATRIX(f->current.v, i, j);
				bs_v[i] += sqr(current-kinv);
				corr[i] += (MATRIX(f->current.u, i, j)-kinv)*(MATRIX(f->current.v, i, j)-kinv);
				diff[i] += sqr(MATRIX(f->current.u, i, j)-MATRIX(f->current.v, i, j));
			}
			corr[i] /= (sqrt(bs_u[i])*sqrt(bs_v[i]));
			diff[i] /= f->k;
		} else {
			bs_v[i] = bs_u[i];
			corr[i] = 1.0;
			diff[i] = 0.0;
		}
		bs_v[i] = 1.0 - sqrt(f->k / (f->k-1.0)) * sqrt(bs_v[i]);
		if (bs_v[i]<0) bs_v[i]=0; /* correction because of rounding errors */
        bs_u[i] = 1.0 - sqrt(f->k / (f->k-1.0)) * sqrt(bs_u[i]);
        if (bs_u[i]<0) bs_u[i]=0; /* correction because of rounding errors */
    }

    fprintf(out, "%% Degree mean, sd: %.4f +- %.4f\n", d_mean, d_std);

	{
		igraph_vector_t v;
		igraph_vector_init(&v, f->k);
		state_get_community_sizes(&f->current, &v);
		fprintf(out, "%% Cluster sizes: [ %.4f", VECTOR(v)[0]);
		for (i=1; i < f->k; i++) fprintf(out, ", %.4f", VECTOR(v)[i]);
		fprintf(out, " ]\n");
		igraph_vector_destroy(&v);
	}

    fprintf(out, "%% Goal function: %.7f\n", f->current.goal);
    fprintf(out, "%% Modularity: %.7f\n", fuzzy_clustering_modularity(f));

    fprintf(out, "Vertex\t");
	for (i=0; i < f->k; i++) fprintf(out, "u_%ld\t", i);
    if (f->is_directed)
		for (i=0; i < f->k; i++) fprintf(out, "v_%ld\t", i);
    fprintf(out, "Degree\t");
    fprintf(out, "Dom\t");
    if (f->is_directed) fprintf(out, "vDom\t");
    fprintf(out, "Relev\t");
    if (f->is_directed) fprintf(out, "vRelev\t");
	fprintf(out, "Bridge\t");
    if (f->is_directed) fprintf(out, "vBridge\tCorr\tSqDiff\t");
	fprintf(out, "\n");

    for (i=0; i < n; i++) {
        /* Cluster profile */
        fprintf(out, "%ld\t", i);
        for (j = 0; j < f->k; j++)
            fprintf(out, "%.4f\t", MATRIX(f->current.u, i, j));
        if (f->is_directed)
			for (j = 0; j < f->k; j++)
				fprintf(out, "%.4f\t", MATRIX(f->current.v, i, j));
        
        /* Degree */
        fprintf(out, "%ld\t", (long)VECTOR(f->outdegrees)[i]);

        /* Dominant cluster */
        if (f->dominance_threshold < 0) {
            for (j=1, k=0; j < f->k; j++)
                if (MATRIX(f->current.u, i, j) > MATRIX(f->current.u, i, k)) k=j;
			if (f->is_directed) {
				for (j=1, l=0; j < f->k; j++)
					if (MATRIX(f->current.v, i, j)>MATRIX(f->current.v, i, l)) l=j;
			}
        } else {
            for (j=0, k=-1; j < f->k; j++)
                if (MATRIX(f->current.u, i, j) > f->dominance_threshold) {
                    if (k >= 0) { k = -1; break; } else k = j;
                }
			if (f->is_directed) {
				for (j=0, l=-1; j < f->k; j++)
					if (MATRIX(f->current.v, i, j) > f->dominance_threshold) {
						if (l >= 0) { l = -1; break; } else l = j;
					}
			}
        }
        if (k == -1) fprintf(out, "-\t"); else fprintf(out, "%d\t", (int)k);
        if (f->is_directed) {
			if (l == -1) fprintf(out, "-\t"); else fprintf(out, "%d\t", (int)l);
		}

        /* Number of relevant clusters */
        delta_u = 0.0;
        for (j=0; j < f->k; j++) {
            s = MATRIX(f->current.u, i, j);
            if (s > 0) delta_u -= s*log(s);
        }
		delta_u = exp(delta_u);
        fprintf(out, "%.4f\t", delta_u);
        if (f->is_directed) {
			delta_v = 0.0;
			for (j=0; j < f->k; j++) {
				s = MATRIX(f->current.v, i, j);
				if (s > 0) delta_v -= s*log(s);
			}
			delta_v = exp(delta_v);
			fprintf(out, "%.4f\t", delta_v);
		} else delta_v = delta_u;

        /* Bridgeness */
        fprintf(out, "%.4f\t", bs_u[i]);
        if (f->is_directed) {
			fprintf(out, "%.4f\t", bs_v[i]);
			fprintf(out, "%.4f\t", corr[i]);
			fprintf(out, "%.4f\t", diff[i]);
		}

        fprintf(out, "\n");
    }

    free(bs_u); free(bs_v); free(corr); free(diff);

    return 0;
}

void fuzzy_clustering_destroy(fuzzy_clustering_t *f) {
    free(f->edges);
    state_destroy(&f->current);
    state_destroy(&f->best);
    state_destroy(&f->local_best);
    igraph_vector_destroy(&f->outdegrees);
    igraph_vector_destroy(&f->indegrees);
}

void usage(int argc, char* argv[], igraph_bool_t long_help) {
    printf("Fuzzy graph partitioning " VERSION "\n\n");
    printf("Usage: %s [options] filename\n", argv[0]);
    printf("\n");
    printf("filename refers to the file containing the graph to be processed.\n");
    printf("Supported formats: edge list (.txt), Pajek (.net), NCOL (.ncol),\n");
    printf("LGL (.lgl), GML (.gml), GraphML (.graphml)\n\n");

    if (!long_help) {
        printf("Options:\n\n");
        printf("    -h: short help (this message)\n");
        printf("    -H, --help: detailed help\n");
        printf("    -c CLUSTERS: number of clusters\n");
        printf("    -f N: pruning frequency\n");
        printf("    -p [THRESHOLD]: pruning threshold\n");
        printf("    -q: be quiet\n");
        printf("    -t THRESHOLD: dominance threshold\n");
        printf("    -w WEIGHTS: load weights from the given file\n");
        return;
    }

    printf("Options:\n\n");
    printf("    -h:\n");
    printf("        displays a short help message.\n");
    printf("    -H, --help:\n");
    printf("        displays a long help message.\n");
    printf("    -c CLUSTERS, --clusters CLUSTERS:\n");
    printf("        the number of clusters to use. If omitted, the algorithm\n");
    printf("        tries to determine the number of clusters automatically\n");
    printf("        by progressively increasing the number of possible\n");
    printf("        clusters whenever it seems to be stuck in a local\n");
    printf("        minimum of the goal function.\n");
    printf("    -f N, --prune-frequency N:\n");
    printf("       prune the membership matrix every Nth step during\n");
    printf("       optimization as well (not only at the end of the process).\n");
    printf("       Use -p to specify the pruning threshold. This additional\n");
    printf("       pruning is usually not necessary, but it may be used for\n");
    printf("       speeding up the optimization process.\n");
    printf("    -p [THRESHOLD], --prune [THRESHOLD]:\n");
    printf("        when the optimization process ended, prune the cluster\n");
    printf("        profile matrix by setting elements less than the given\n");
    printf("        threshold to zero. The rows of the matrix are then\n");
    printf("        re-scaled to ensure that the row sums are equal to 1.\n");
    printf("        If this option is given, but the threshold is omitted,\n");
    printf("        pruning is turned off. The default pruning level (if\n");
    printf("        this option is not given) is 10 times the number of\n");
    printf("        clusters.\n");
    printf("    -q:\n");
    printf("        be quiet, don't write anything to stdout\n");
    /*printf("    -s SIMILARITY, --similarity SIMILARITY:\n");
    printf("        use the given target similarity assumption for the\n");
    printf("        vertices. SIMILARITY can be one of: adjacency (use the\n");
    printf("        adjacency matrix), jaccard (Jaccard similarity index)\n");*/
	printf("    -S SEED, --seed SEED:\n");
	printf("        use the given random SEED. Use this option to ensure that\n");
	printf("        the results are the same on the same machine no matter how\n");
	printf("        many times you run the executable.\n");
    printf("    -t THRESHOLD, --dominance-threshold THRESHOLD:\n");
    printf("        when determining dominant clusters for each vertex, a\n");
    printf("        cluster is considered dominant if the vertex belongs\n");
    printf("        to it to a greater extent than this threshold. If more\n");
    printf("        than one dominant cluster exists for the vertex, it is\n");
    printf("        considered a bridge. If a negative threshold is given,\n");
    printf("        all vertices will have a dominant cluster: the one to\n");
    printf("        which they belong with the greatest degree of membership\n");
    printf("    -w [FILE], --weight-file [FILE]:\n");
    printf("        load constraint weights from the given file.\n");
    printf("        Each line in the file must contain a whitespace-separated\n");
    printf("        list consisting of the following values: vertex 1, vertex 2,\n");
    printf("        weight (optional, defaults to 1). This option can't be used\n");
    printf("        together with -d. Vertex pairs not specified in the\n");
    printf("        file will not have an assigned constraint. Make sure that\n");
    printf("        every vertex pair appears only once in the file.\n");
}

int process_command_line(int argc, char* argv[], parameters_t *params) {
    int ch, idx;

    static const struct option opts[] = {
        {"clusters", required_argument, 0, 'c'},
        {"help", no_argument, 0, 'H'},
        {"prune-frequency", required_argument, 0, 'f'},
        {"prune", required_argument, 0, 'p'},
        {"quiet", no_argument, 0, 'q'},
		{"seed", required_argument, 0, 'S'},
        {"similarity", required_argument, 0, 's'},
        {"dominance-threshold", required_argument, 0, 't'},
        {"weight-file", optional_argument, 0, 'w'},
        {0, 0, 0, 0}
    };
    params->no_of_clusters = 2;
    params->dominance_threshold = -1;
	params->seed = 0;
	params->has_seed = 0;
    params->prune_threshold = -1;
    params->prune_frequency = 0;
    params->verbosity = 1;
    params->weight_file = 0;
    params->similarity_type = SIMILARITY_ADJACENCY;
    params->adaptive_cluster_count = 1;

    while ((ch = getopt_long(argc, argv, "c:f:hHp:qs:S:t:w::", opts, &idx)) != -1) {
        switch (ch) {
            case 'c':
                /* Number of clusters */
                params->adaptive_cluster_count = 0;
                params->no_of_clusters = atoi(optarg);
                break;
            case 'f':
                /* Set pruning frequency */
                params->prune_frequency = atoi(optarg);
                break;
            case 'p':
                /* Set prune threshold */
                if (optarg)
                    params->prune_threshold = atof(optarg);
                else
                    params->prune_threshold = 0;
                break;
            case 'q':
                /* Be quiet */
                params->verbosity = 0;
                break;
            case 's':
                /* Similarity type */
                if (!strcmp(optarg, "adjacency"))
                  params->similarity_type = SIMILARITY_ADJACENCY;
                else if (!strcmp(optarg, "jaccard"))
                  params->similarity_type = SIMILARITY_JACCARD;
                else
                  FATAL("unknown similarity type");
                break;
			case 'S':
				/* Random seed */
				params->seed = atoi(optarg);
				params->has_seed = 1;
				break;
            case 't':
                params->dominance_threshold = atof(optarg);
                break;
            case 'w':
                if (optarg)
                    params->weight_file = optarg;
                else
                    params->weights_from_modularity = 1;
                break;
            case 'h':
                usage(argc, argv, 0);
                return 1;
            case 'H':
                usage(argc, argv, 1);
                return 1;
            default:
                return 1;
        }
    }

    if (params->prune_threshold == -1)
      params->prune_threshold = 0.05;

    if (params->no_of_clusters < 2) 
        FATAL("number of clusters must be at least 2\n");

    if (argv[optind]) params->filename = argv[optind];
    else FATAL("%s: no filename given\n", argv[0]);
    
    return 0;
}

int main(int argc, char* argv[]) {
    igraph_t g;
    fuzzy_clustering_t clustering;
    igraph_real_t diff;
    FILE* f;
    char* ext;
    int error_code;
    int tries, max_tries=5, max_cluster_tries=1, i;
    igraph_real_t start_step_size=0.25;
    long int steps=0;
    clock_t start_time, end_time;
    igraph_real_t q, prev_q = -666;

    error_code=process_command_line(argc, argv, &params);
    if (error_code) return error_code;

    srand(params.has_seed ? params.seed : time(0));

    /*IGRAPH_CHECK(igraph_small(&g, 5, 0, 0, 1, 0, 2, 1, 2, 2, 3, 2, 4, 3, 4,-1));*/
    if (strcmp(params.filename, "-"))
        f=fopen(params.filename, "r");
    else
        f=stdin;
    if (!f) { perror(params.filename); return 1; }

    /* Identify the extension and use the appropriate loader */
    if (f != stdin) {
        ext = strrchr(params.filename, '.');
        if (!ext) FATAL("%s: unknown file type\n", params.filename);
        if (!strcmp(ext, ".graphml")) {
            INFO("Loading graph from GraphML file: %s\n", params.filename);
            IGRAPH_CHECK(igraph_read_graph_graphml(&g, f, 0));
        } else if (!strcmp(ext, ".gml")) {
            INFO("Loading graph from GML file: %s\n", params.filename);
            IGRAPH_CHECK(igraph_read_graph_gml(&g, f));
        } else if (!strcmp(ext, ".net")) {
            INFO("Loading graph from Pajek file: %s\n", params.filename);
            IGRAPH_CHECK(igraph_read_graph_pajek(&g, f));
        } else if (!strcmp(ext, ".lgl")) {
            INFO("Loading graph from LGL file: %s\n", params.filename);
            IGRAPH_CHECK(igraph_read_graph_lgl(&g, f, 0, 0, IGRAPH_UNDIRECTED));
        } else if (!strcmp(ext, ".ncol")) {
            INFO("Loading graph from NCOL file: %s\n", params.filename);
            IGRAPH_CHECK(igraph_read_graph_ncol(&g, f, 0, 0, 0, IGRAPH_UNDIRECTED));
        } else if (!strcmp(ext, ".txt")) {
            INFO("Loading graph from edge list: %s\n", params.filename);
            IGRAPH_CHECK(igraph_read_graph_edgelist(&g, f, 0, 0));
        } else FATAL("%s: unknown file type\n", params.filename);
        fclose(f);
    } else {
        INFO("Loading graph from edge list: %s\n", params.filename);
        IGRAPH_CHECK(igraph_read_graph_edgelist(&g, f, 0, 0));
    }

    IGRAPH_FINALLY(igraph_destroy, &g);
    if (igraph_is_directed(&g)) {
		INFO("Converting directed graph to undirected\n");
		igraph_to_undirected(&g, IGRAPH_TO_UNDIRECTED_COLLAPSE, 0);
    } else {
        INFO("Simplifying graph\n");
        igraph_simplify(&g, 1, 1, 0);
    }
    INFO("Graph has %ld vertices and %ld edges\n", (long)igraph_vcount(&g), (long)igraph_ecount(&g));
    INFO("Initializing constraints\n");
    IGRAPH_CHECK(fuzzy_clustering_init(&clustering, &g, 2, params.dominance_threshold));
    IGRAPH_FINALLY(fuzzy_clustering_destroy, &clustering);

	tries = (params.adaptive_cluster_count || params.no_of_clusters > clustering.k) ? max_cluster_tries : max_tries;
    start_time = clock();
    while (clustering.step_size>clustering.epsilon) {
        int has_large_gradient = clustering.iterate(&clustering,
		  params.prune_frequency>0 && ((steps+1) % params.prune_frequency==0));

        steps += 1;
        diff = clustering.current.goal - clustering.local_best.goal;
        if (diff >= 0) {
            if (clustering.step_size>clustering.epsilon) clustering.step_size /= 2;
            clustering.successes = 0;
            state_copy(&clustering.current, &clustering.local_best);
            INFO("- [%.4f] %.4f\r", clustering.step_size, clustering.current.goal);
        } else {
            clustering.successes++;
			tries = (params.adaptive_cluster_count || params.no_of_clusters > clustering.k) ? max_cluster_tries : max_tries;
            INFO("+ [%.4f] %.4f\r", clustering.step_size, clustering.current.goal);
            state_copy(&clustering.local_best, &clustering.current);
            if (clustering.successes >= 3) {
                clustering.step_size *= 1.5;
                clustering.successes = 0;
            }
        }
        if (clustering.step_size <= clustering.epsilon || !has_large_gradient) {
            if (tries > 0 && has_large_gradient) {
                clustering.step_size *= 10;
                fuzzy_clustering_mutate(&clustering, 0.2);
                tries--;
            } else {
                if (params.adaptive_cluster_count || params.no_of_clusters > clustering.k) {
					tries = (params.adaptive_cluster_count || params.no_of_clusters > clustering.k) ? max_cluster_tries : max_tries;
					if (params.prune_threshold > 0) {
						fuzzy_clustering_prune(&clustering, params.prune_threshold);
					}
                    q = fuzzy_clustering_modularity(&clustering);
                    INFO("Modularity is now %.4f (%ld clusters)\n", q, clustering.k);
                    if (q-prev_q < 0.0001 && params.adaptive_cluster_count) {
                        state_copy(&clustering.current, &clustering.best);
                    } else {
                        INFO("Increasing number of clusters to %ld\n", clustering.k+1);
                        prev_q = q;
                        state_copy(&clustering.best, &clustering.current);
                        fuzzy_clustering_add_community(&clustering, 1);
                        clustering.step_size = start_step_size;
                        clustering.successes = 0;
                        fuzzy_clustering_mutate(&clustering, 0.1);
                        state_copy(&clustering.local_best, &clustering.current);
                    }
                }
            }
        }
    }
    INFO("Convergence achieved in %ld steps\n", (long)steps);

    if (params.prune_threshold > 0) {
        INFO("Pruning matrix with threshold = %.4f\n", params.prune_threshold);
        fuzzy_clustering_prune(&clustering, params.prune_threshold);
    }
    end_time = clock();

	// it can happen that some of the communities are (almost) empty. We remove
	// them now
	i = clustering.k;
	fuzzy_clustering_remove_empty_communities(&clustering, 1.0);
	if (clustering.k < i) {
		INFO("Removed %ld empty clusters\n", (i - clustering.k));
	}

    WARNING("CPU time used: %.4fs\n", ((double)(end_time-start_time))/CLOCKS_PER_SEC);
    fuzzy_clustering_recalculate(&clustering);
    fuzzy_clustering_print(&clustering, stdout);

    fuzzy_clustering_destroy(&clustering);
    igraph_destroy(&g);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

