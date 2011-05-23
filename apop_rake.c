//#define __USE_POSIX //for strtok_r
#include "apop_internal.h"
#include <gsl/gsl_sort_vector.h>

/* This is the internal documentation for apop_rake(). I assume you've read the usage
   documentation already (if you haven't, it's below, just above apop_rake).
  
  I started with:
Algorithm AS 51 Appl. Statist. (1972), vol. 21, p. 218
original (C) Royal Statistical Society 1972

But at this point, I'm not sure if any of the original code remains at all, because the
sparse method and the full-matrix method differ so much. 

There are two phases to the process: the SQL part and the in-memory table part. The file
here begins with the indexing of the in-memory table. The index on the in-memory
table is somewhat its own structure, with get/set functions and so forth, so it has
its own section, listed first in this file. After that, we get to the raking function
(c_loglin) and its supporting functions, and then the apop_rake function itself. The
SQL work itself is inside apop_rake, and c_loglin is called at the end to do the
raking. */


/* This section indexes a PMF-type apop_data struct. The index is held in a 2-d grid
of nodes. The index's first index is the dimension, and will hold one column for each
column in the original data set. The index's second index goes over all of the values
in the given column.

PS: this was intended to one day become a general-purpose index for PMFs; making that
happen remains on the to-do list.

mnode[i] = a dimension row
mnode[i][j] = a value in a given dimension
mnode[i][j].margin_ptrs = a list of all of the rows in the data set with the given value.
mnode[i][j].margin_ptrs[k] = the kth item for the value.
  */

//// If the index part of this file were split into separate files, the next few lines would be index.h.
typedef struct {
    double val;
    int len;
    int *margin_ptrs;
} mnode_t;

typedef void(*index_apply_f)(mnode_t**, int, void*);

mnode_t **index_generate(const apop_data *in);
void index_free(mnode_t **in);
void index_add_node(mnode_t **mnodes, size_t dim, size_t row, double val);
void index_foreach(mnode_t *index[], index_apply_f f, void *args);
void index_get_element_list(mnode_t **index, int *d);

//End index.h; begin index.c

static int find_val(double findme, mnode_t *nodecol){
    for (int i=0; nodecol[i].val <= findme; i++)
        if (nodecol[i].val == findme)
           return i;
    Apop_assert(0, "I can't find a value that should've been inserted.");
}

static int size;

void index_add_node(mnode_t **mnodes, size_t dim, size_t row, double val){
    int index = find_val(val, mnodes[dim]);
    mnodes[dim][index].margin_ptrs[row] = 1;
}

mnode_t **index_generate(const apop_data *in){
    size_t margin_ct = in->matrix->size2;
    size = in->matrix->size1;
    mnode_t **mnodes = malloc(sizeof(mnode_t*)*(margin_ct+1));
    //allocate every node
    for(size_t i=0; i < margin_ct; i ++){
        Apop_col(in, i, col);
        gsl_vector *vals = apop_vector_unique_elements(col);
		gsl_sort_vector(vals);
        mnodes[i] = malloc(sizeof(mnode_t)*(vals->size+1));
        for(size_t j=0; j < vals->size; j ++)
            mnodes[i][j] = (mnode_t) {.val = gsl_vector_get(vals, j), .len=0, 
                        .margin_ptrs = calloc(size, sizeof(int))};
        mnodes[i][vals->size] = (mnode_t) {.val = GSL_POSINF}; //end-of-array sentinel
        gsl_vector_free(vals);
    }
    mnodes[margin_ct] = NULL; //end-of-array sentinel
    //put data from the matrix into the right pigeonhole
    for(size_t i=0; i < margin_ct; i++)
        for(size_t j=0; j < in->matrix->size1; j++)
            index_add_node(mnodes, i, j, apop_data_get(in, j, i));
    return mnodes;
}

void index_free(mnode_t **in){
	for (int i=0; in[i]; i++){
		//for (int j=0; isfinite(in[i][j].val); j++)
		/*for (int j=0; in[i][j].margin_ptrs; j++)
			free(in[i][j].margin_ptrs);*/
        free(in[i]);
	}
    free(in);
}

/* The next two functions are a recursive iteration over all combinations of values
for a given index (which may be the whole thing or a margin). index_foreach does the
initialization of some state variables; value_loop does the odometer-like recursion. At
each step, value_loop will either increment the current dimension's index, or if the
current index is at its limit, will loop back to zero on this index and then set the
next dimension as active.  
 */
static void value_loop(mnode_t *icon[], int *indices, mnode_t **values, 
            int this_dim, index_apply_f f, int *ctr, void *args){
    while (isfinite(icon[this_dim][indices[this_dim]].val)){
        values[this_dim] = &icon[this_dim][indices[this_dim]];
		if (icon[this_dim+1]){
			indices[this_dim+1]=0;
			value_loop(icon, indices, values, this_dim+1, f, ctr, args);
		} else{
			f(values, *ctr, args);
			(*ctr)++;
		}
        indices[this_dim]++;
	}
}

/* Inputs: the index to be iterated over, the function to apply to each combination,
and a void* with other sundry arguments to the function. I'll apply the function 
to an mnode_t* with a single mnode_t for each dimension.
 */
void index_foreach(mnode_t *index[], index_apply_f f, void *args){
    int j, ctr=0;
    for (j=0; index[j]; j++) ;//Do nothing; just counting elmts.
    mnode_t *values[j+1];
    for (j=0; index[j]; j++) 
        values[j] = &index[j][0];
    values[j] = malloc(sizeof(mnode_t));
    values[j]->val = GSL_POSINF;
    int indices[j];
    memset(indices, 0, sizeof(int)*j);
    value_loop(index, indices, values, 0, f, &ctr, args);
    free(values[j]);
}

/* Check whether an observation falls into all of the given margins.

\param index Is actually a partial index: for each dimension, there should be only one value. Useful for the center of an index_foreach loop.
\param d Should already be allocated to the right size, may be filled with garbage.
\return d will be zero or one to indicate which rows of the indexed data set meet all criteria.
*/
void index_get_element_list(mnode_t **index, int *d){
    memcpy(d, index[0]->margin_ptrs, size *  sizeof(int));
    for(int i=1; isfinite(index[i]->val); i++)
        for(int j=0; j < size; j++)
            //Next two lines are equivalent to d[j]=d[j] && index[i]->..., but clock in as faster.
            if (d[j])
                d[j] = index[i]->margin_ptrs[j];
}

////End index.c

typedef struct {
    const apop_data *indata; 
    apop_data *fit; 
    size_t **elmtlist;
    size_t *elmtlist_sizes;
	gsl_vector *indata_values;
    mnode_t **index;
    size_t ct, al;
} rake_t;

void grow(rake_t *r){
    r->al = (r->al+1)*2;
    r->elmtlist = realloc(r->elmtlist , sizeof(size_t*) * r->al);
    r->elmtlist_sizes = realloc(r->elmtlist_sizes, sizeof(size_t) * r->al);
    r->indata_values = apop_vector_realloc(r->indata_values, r->al);
}

void rakeinfo_free(rake_t r){
    for (int i=0; i < r.ct; i++){
        free(r.elmtlist[i]);
        r.elmtlist[i] = NULL;
    }
    free(r.index);
    free(r.elmtlist_sizes);
    r.index = NULL;
    r.elmtlist_sizes = NULL;
/*    free(r.elmtlist);
    r.elmtlist = NULL;*/
    gsl_vector_free(r.indata_values);
    r.indata_values= NULL;
}

double overall_max_dev;

void scaling(const size_t *elmts, size_t n,  gsl_vector *weights, double in_sum){
    double fit_sum = 0;
    for(size_t i=0; i < n; i ++)
        fit_sum += weights->data[elmts[i]];
    for(size_t i=0; i < n; i ++){
        weights->data[elmts[i]] *= in_sum/fit_sum;
	}
	overall_max_dev = GSL_MAX(overall_max_dev, fabs(in_sum-fit_sum));
}

/* Given one set of values from one margin, do the actual scaling.
 On the first pass, this function takes notes on each margin's element list and total 
 in the original data. Later passes just read the notes and call the scaling() function above.
*/
static void one_set_of_values(mnode_t **icon, const int ctr, void *in){
    rake_t *r = in;
    int size = r->indata->matrix->size1;
    static int *t = NULL;
    if (!t) t = malloc(size * sizeof(int));
	int first_pass = 0;
    double in_sum;
	if (ctr < r->ct)
		in_sum = gsl_vector_get(r->indata_values, ctr);
    else {
        r->ct++;
        if (ctr >= r->al) 
            grow(r);
   		index_get_element_list(icon, t);
        in_sum = 0;
        int n=0, al=0;
        r->elmtlist[ctr] = NULL;
        for(int m=0; m < size; m++)
            if (t[m]){
                in_sum += r->indata->weights->data[m];
                if (n >= al) {
                    al = (al+1)*2;
                    r->elmtlist[ctr] = realloc(r->elmtlist[ctr], al*sizeof(size_t));
                }
                r->elmtlist[ctr][n++] = m;
            }
        r->elmtlist_sizes[ctr] = n;
        r->indata_values->data[ctr] = in_sum;
		first_pass++;
	}
    if (!r->elmtlist_sizes[ctr]) return;
    if (!first_pass && !in_sum)  return;
    scaling(r->elmtlist[ctr], r->elmtlist_sizes[ctr], r->fit->weights, in_sum);
}

/* For each configuration margin, for each combination for that margin, 
   call the above one_set_of_values() function. */
static void main_loop(int config_ct, rake_t *rakeinfo, int k){
	overall_max_dev = GSL_NEGINF;
    for(size_t i=0; i < config_ct; i ++)
		if (k==1)
			index_foreach(rakeinfo[i].index, one_set_of_values, rakeinfo+i);
		else
			for(int m=0; m < rakeinfo[i].ct; m++)
				one_set_of_values(NULL, m, rakeinfo+i);
}

/* Following the FORTRAN, 1 contrast ==> icon. Here, icon will be a
subset of the main index including only the columns pertaining to a given margin. */
void generate_margin_index(mnode_t **icon, const apop_data *margin, mnode_t **mainindex, size_t col){
    Apop_col(margin, col, iconv)
    int ct = 0;
    for (int j=0; mainindex[j]; j++)
        if (gsl_vector_get(iconv, j))
            icon[ct++] = mainindex[j];
    icon[ct] = NULL;
}

void cleanup(mnode_t **index, rake_t rakeinfos[], int contrast_ct){
	for(size_t i=0; i < contrast_ct; i++)
		rakeinfo_free(rakeinfos[i]);
	index_free(index);
}

/*
\param config 	An nvar x ncon matrix; see below. [as in the original, but not squashed into 1-D.]
\param indata 	the actual table. I use a PMF format.
\param fit 		the starting table. Same size as table.
\param maxdev 	maximum deviation; stop when this is met.
\param maxit 	maximum iterations; stop when this is met.
   
Re: the contrast array: Each _column_ is a contrast. Put a one in each col
 involved in the contrast. E.g., for (1,2), (2,3):

 1 0<br>
 1 1<br>
 0 1
 */
static void c_loglin(const apop_data *config, const apop_data *indata, apop_data *fit, double maxdev, int maxit) {
    mnode_t ** index = index_generate(indata);

    /* Make a preliminary adjustment to obtain the fit to an empty configuration list */
    //fit->weights is either all 1 (if no weights_col) or the initial counts from the db.
    double x = apop_vector_sum(indata->weights);
    double y = apop_sum(fit->weights);
    gsl_vector_scale(fit->weights, x/y);

	int contrast_ct =config && config->matrix ? config->matrix->size2 : 0;
    rake_t rakeinfos[contrast_ct];
    for(size_t i=0; i < contrast_ct; i ++){
        Apop_col(config, i, iconv)
        rakeinfos[i] = (rake_t) {
            .indata = indata,
            .fit = fit,
            .index = malloc(sizeof(mnode_t) *(apop_sum(iconv)+1))
            //others are NULL, to be filled in as we go.
        };
        generate_margin_index(rakeinfos[i].index, config, index, i);
    }
    int k;
    for (k = 1; k <= maxit; ++k) {
		if (!(k%100)){printf(".");fflush(NULL);}
        main_loop(contrast_ct, rakeinfos, k);
        if (overall_max_dev < maxdev) { // Normal termination 
            cleanup(index, rakeinfos, contrast_ct);
            return;
        }
    }
    cleanup(index, rakeinfos,contrast_ct);
    Apop_assert_c(k < maxit, , 0, "Maximum number of iterations reached. The max "
            "deviation is %g, which is still more than tolerance = %g", overall_max_dev, maxdev);
}
char *pipe_parse = "[ \n\t]*([^| \n\t]+)[ \n\t]*([|]|$)";

apop_data **generate_list_of_contrasts(char **contras_in, int contrast_ct){
  apop_data** out = malloc(sizeof(apop_data*)* contrast_ct);
	for (int i=0; i< contrast_ct; i++) 
        apop_regex(contras_in[i], pipe_parse, out+i);
	return out;
}

int get_var_index(apop_data *all_vars, char *findme){
	for (int i=0; i< all_vars->textsize[0]; i++)
		if (apop_strcmp(all_vars->text[i][0], findme))
			return i;
	Apop_assert_c(0, -1, 0, "I couldn't find %s in the full list of variables. Returning -1.", findme);
}

/** Fit a log-linear model via iterative proportional fitting, aka raking.

See Wikipedia for an overview of Log linear models, aka
<a href="http://en.wikipedia.org/wiki/Poisson_regression">Poisson regressions</a>. 
One approach toward log-linear modeling is a regression form; let there be four
categories, A, B, C, and D, from which we can produce a model positing, for example,
that cell count is a function of a form like \f$g_1(A) + g_2(BC) + g_3(CD)\f$. In this case, we would
assign a separate coefficient to every possible value of A, every possible value of
(B, C), and every value of (C, D). Raking is the technique that searches for that large
set of parameters.

The combinations of categories that are considered to be relevant are called \em
contrasts, after ANOVA terminology of the 1940s.

The other constraint on the search are structural zeros, which are values that you know
can never be non-zero, due to field-specific facts about the variables. For example, U.S.
Social Security payments are available only to those age 65 or older, so "age <65 and
gets_soc_security=1" is a structural zero.

Because there is one parameter for every combination, there may be millions of parameters
to estimate, so the search to find the most likely value requires some attention to
technique. For over half a century, the consensus method for searching has been raking, which
iteratively draws each category closer to the mean in a somewhat simple manner (this was
first developed circa 1940 and had to be feasible by hand), but which is guaranteed to
eventually arrive at the maximum likelihood estimate for all cells.

Another complication is that the table is invariably sparse. One can easily construct
tables with millions of cells, but the corresponding data set may have only a few
thousand observations.

This function uses the database to resolve the sparseness problem. It constructs a query
requesting all combinations of categories the could possibly be non-zero after raking,
given all of the above constraints. Then, raking is done using only that subset. 
This means that the work is done on a number of cells proportional to the number of data
points, not to the full cross of all categories.

\param table_name The name of the table in the database.  The table should have one
observation per row.
No default.

\param all_vars The full list of variables to search. Provide these as a single,
pipe-delimited string: e.g., "age | sex |income". Spaces are ignored. Default: if you are using SQLite, I
will use all columns in the table; if you are using mySQL, I haven't implemented this yet
and you will have to provide a list.

\param contrasts The contrasts describing your model. Like the \c all_vars input, each
contrast is a pipe-delimited list of variable names. No default.

\param contrast_ct The number of contrasts in the list of contrasts. No default.

\param structural_zeros a SQL clause indicating combinations that can never take a nonzero
value. This will go into a \c where clause, so anything you could put there is OK, e.g.
"age <65 and gets_soc_security=1 or age <15 and married=1". Default: no structural zeros.
\param max_iterations Number of rounds of raking at which the algorithm halts. Default: 1000.

\param tolerance I calculate the change for each cell from round to round;
if the largest cell change is smaller than this, I stop. Default: 1e-5.

\param weights_col This column gives the count of how many observations are represented
by each row. If \c NULL, ech row represents one person. Default: \c NULL.

\param run_number Because I write intermediate tables to the database, I need a way to
distinguish distinct runs should you be threading several runs at once. If you aren't
running several instances simultaneously, don't worry about this; if you are, do supply a
value, since it's hard for the function to supply one in a race-proof manner. Default:
internally-maintained values.

\return An \ref apop_data set where every row is a single combination of variable values
and the \c weights vector gives the most likely value for each cell.

\li This function uses the \ref designated syntax for inputs.

\li The interface is still beta, and subject to change---notably, handling of text
categories will soon be added.

\li Legal: this function is part of a larger work (C) Ben Klemens, but was partially
written by a U.S. government employee during work hours. Some lawyers will tell you that
the code is licensed via the same modified GPL v2 as the main work; others will tell
you that it is public domain.
*/
APOP_VAR_HEAD apop_data * apop_rake(char *table_name, char *all_vars, char **contrasts, int contrast_ct, char *structural_zeros, int max_iterations, double tolerance, char *weights_col, int run_number){
    static int defaultrun = 0;
    char * apop_varad_var(table_name, NULL);
    Apop_assert(table_name, "I need the name of a table in the database that will be the data source.");
    char * apop_varad_var(all_vars, NULL);
    char ** apop_varad_var(contrasts, NULL); //default to all vars?
    int apop_varad_var(contrast_ct, 0);
    char * apop_varad_var(structural_zeros, NULL);
    char * apop_varad_var(weights_col, NULL);
    int apop_varad_var(max_iterations, 1e3);
    double apop_varad_var(tolerance, 1e-5);
    int apop_varad_var(run_number, defaultrun++);
APOP_VAR_ENDHEAD
	apop_data **contras = generate_list_of_contrasts(contrasts, contrast_ct);
    apop_data *all_vars_d;
    if (!all_vars){
        Apop_assert(apop_opts.db_engine!='m', "I need a list of the full set of variable "
                                            "names sent as .all_vars=\"var1 | var2 |...\"");
        //use SQLite's table_info, then shift the second col to the first.
        all_vars_d = apop_query_to_text("PRAGMA table_info(%s)", table_name);
        for (int i=0; i< all_vars_d->textsize[0]; i++)
            apop_text_add(all_vars_d, i, 0, all_vars_d->text[i][1]);
    }
    apop_regex(all_vars, pipe_parse, &all_vars_d);
    int var_ct = all_vars_d->textsize[0];

	char *q;
	apop_query("drop table if exists apop_zerocontrasts_%i", run_number);
	apop_query("drop table if exists apop_contrasts_%i", run_number);
	asprintf(&q, "create table apop_zerocontrasts_%i as select ", run_number);
	//the asprintf(&q, "...", q) leaks; I don't care.
     /*Start with every possible combination:
         " select distinct b.block, qa.qageshort, qs.qsex, r.racep, hi.qspanq, 0"
         " from "
 		" (select distinct block as block from d) as b, "
 		" (select distinct qageshort from d) qa, " ...*/
	for (int i=0; i < var_ct; i++)
		asprintf(&q, "%s t%i.%s, ", q, i, all_vars_d->text[i][0]); 
	asprintf(&q, "%s 0\n from\n", q); 
	char comma = ' ';
	for (int i=0; i < var_ct; i++){
	  	asprintf(&q, "%s %c (select distinct %s as %s from %s) as t%i\n", 
					  q, comma, all_vars_d->text[i][0], all_vars_d->text[i][0], table_name, i);
		comma = ',';
	}
     /*keep only rows that could have nonzero values in the design matrix margins (about 10%):
 		 where b.block||'.'||qa.qageshort||'.' || qs.qsex ||'.'|| r.racep in 
 		 (select block ||'.'|| qageshort ||'.'|| qsex ||'.'|| racep from d) and ...*/
    if (contrast_ct){
        asprintf(&q, "%s\nwhere\n", q);
        char joiner[] = "||'.'||";
        char space[] = " ";
        for (int i=0; i< contrast_ct; i++){
            char *merge = strdup(" ");
            for (int j=0; j< contras[i]->textsize[0]; j++){
                asprintf(&merge, "%s t%i.%s %s", 
                    merge, get_var_index(all_vars_d, contras[i]->text[j][0]),
                    contras[i]->text[j][0], j+1!=contras[i]->textsize[0]? joiner : space);
            }
            asprintf(&q, "%s %s in apop_m%i_%i\n ", q, merge, i, run_number);
            if (i+1 < contrast_ct)
                asprintf(&q, "%s and ", q);

            //While we're here, we should create these mi tables:
            merge = strdup(" ");
            for (int j=0; j< contras[i]->textsize[0]; j++){
                asprintf(&merge, "%s %s %s", 
                    merge, 
                    contras[i]->text[j][0], j+1!=contras[i]->textsize[0]? joiner : space);
            }
            apop_query("drop table if exists apop_m%i_%i", i, run_number);
            apop_query("create table apop_m%i_%i as select distinct %s as concatenated from %s; "
                          "create index apop_mi%i_%i on apop_m%i_%i(concatenated);",
                                i,run_number, merge, table_name, i,run_number, i,run_number);
        }
    }
/* Keep out margins with values for now; join them in below.
 	    except  select block, qageshort, qsex, racep, qspanq, 0 from d"; */
	asprintf(&q, "%s except\nselect ", q);
	for (int i=0; i < var_ct; i++)
		asprintf(&q, "%s %s, ", q, all_vars_d->text[i][0]); 
	apop_query("%s 0 from %s", q, table_name);

    char *format=strdup("w");
    for (int i =0 ; i< var_ct; i++)
        asprintf(&format, "m%s", format);
    /*create table contrasts as 
			  select block, qageshort, qsex, racep, qspanq, count(*) from d
              group by block, qageshort, qsex, racep, qspanq
			      union 
			  select * from zerocontrasts 
		  Then,
              delete from contrasts where [structural_zeros]
     */
	asprintf(&q, "create table apop_contrasts_%i as select ", run_number);
	for (int i=0; i < var_ct; i++)
		asprintf(&q, "%s %s, ", q, all_vars_d->text[i][0]); 
    if (weights_col)
        asprintf(&q, "%s sum(%s) from %s\ngroup by ", q, weights_col, table_name);
    else
        asprintf(&q, "%s count(*) from %s\ngroup by ", q, table_name);
	for (int i=0; i < var_ct; i++)
		asprintf(&q, "%s %s%c ", q, all_vars_d->text[i][0], i+1<var_ct ? ',' : ' '); 
	asprintf(&q, "%s\n  union\nselect * from apop_zerocontrasts_%i ", q, run_number);
	apop_query("%s", q);
	if (structural_zeros)
		 apop_query("delete from apop_contrasts_%i where\n %s", run_number, structural_zeros);

    //apop_contrasts... holds the cells of the grid we actually need. Query them to 
    //an apop_data set and start doing the raking.
    apop_data *fit, *contrast_grid;
    fit = apop_query_to_mixed_data(format, "select * from apop_contrasts_%i", run_number);
    Apop_assert(fit, "This query:\n%s\ngenerated a blank table.", q);

    contrast_grid = apop_data_calloc(0,var_ct, contrast_ct);
	for (int i=0; i< contrast_ct; i++)
		for (int j=0; j< contras[i]->textsize[0]; j++)
			apop_data_set(contrast_grid,  get_var_index(all_vars_d, contras[i]->text[j][0]), i, 1);
	//clean up
	for (int i=0; i< contrast_ct; i++){
		apop_query("drop table apop_m%i_%i", i, run_number);
		apop_data_free(contras[i]);
	}

    apop_data *d =apop_data_copy(fit);
    c_loglin(contrast_grid, d, fit, tolerance, max_iterations);
    apop_data_free(d);
    
	apop_query("drop table apop_zerocontrasts_%i", run_number);
	apop_query("drop table apop_contrasts_%i", run_number);
	apop_data_free(contrast_grid);
	return fit;
}

/*  
    handle strings
    */
