
//#define __USE_POSIX //for strtok_r
#include "apop_internal.h"
#include <stdbool.h>
#include <gsl/gsl_sort_vector.h>
void xprintf(char **q, char *format, ...); //in apop_conversions.c

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

typedef struct {
    double val;
    bool *margin_ptrs, *fit_ptrs;
} mnode_t;

typedef void(*index_apply_f)(mnode_t * const * const, int, void*);

static int find_val(double findme, mnode_t *nodecol){
    for (int i=0; nodecol[i].val <= findme || gsl_isnan(findme); i++)
        if (nodecol[i].val == findme || (gsl_isnan(findme) && gsl_isnan(nodecol[i].val)))
           return i;
    return -1;
}

//returns -1 if the given row/col doesn't exist.
int index_add_node(mnode_t **mnodes, size_t dim, size_t row, double val, bool is_margin){
    int index = find_val(val, mnodes[dim]);
    if (index == -1) return -1;
    if (is_margin) mnodes[dim][index].margin_ptrs[row] = true;
    else           mnodes[dim][index].fit_ptrs[row] = true;
    return 0;
}

mnode_t **index_generate(apop_data const *in, apop_data const *in2){
    size_t margin_ct = in->matrix->size2;
    size_t margin_rows = in->matrix->size1;
    size_t fit_rows = in2->matrix->size1;
    mnode_t **mnodes = malloc(sizeof(mnode_t*)*(margin_ct+1));
    //allocate every node
    for(size_t i=0; i < margin_ct; i ++){
        Apop_col_v(in, i, col);
        gsl_vector *vals = apop_vector_unique_elements(col);
        mnodes[i] = malloc(sizeof(mnode_t)*(vals->size+1));
        for(size_t j=0; j < vals->size; j ++)
            mnodes[i][j] = (mnode_t) {.val = gsl_vector_get(vals, j),
                        .margin_ptrs = calloc(margin_rows, sizeof(bool)),
                        .fit_ptrs = calloc(fit_rows, sizeof(bool))
            };
        mnodes[i][vals->size] = (mnode_t) {.val = GSL_POSINF}; //end-of-array sentinel
        gsl_vector_free(vals);
    }
    mnodes[margin_ct] = NULL; //end-of-array sentinel
    //put data from the matrix into the right pigeonhole
    for(size_t i=0; i < margin_ct; i++){
        for(size_t j=0; j < in->matrix->size1; j++)
            Apop_stopif(index_add_node(mnodes, i, j, apop_data_get(in, j, i), true) == -1, return mnodes,
            0, "I can't find a value, %g, that should've already been inserted.", apop_data_get(in, j, i));
        for(size_t j=0; j < in2->matrix->size1; j++) //these values may not be present, in which case ignore them.
            index_add_node(mnodes, i, j, apop_data_get(in2, j, i), false);
    }
    return mnodes;
}

void index_free(mnode_t **in){
	for (int i=0; in[i]; i++){
		for (int j=0; !isinf(in[i][j].val); j++){
			free(in[i][j].margin_ptrs);
			free(in[i][j].fit_ptrs);
        }
        free(in[i]);
	}
    free(in);
}

/* The next two functions are a recursive iteration over all combinations of values
for a given index (which may be the whole thing or a margin). index_foreach does the
initialization of some state variables; value_loop does the odometer-like recursion. At
each step, value_loop will either increment the current dimension's index, or if the
current index is at its limit, will loop back to zero on this index and then set the
next dimension as active. */
static void value_loop(mnode_t *icon[], int *indices, mnode_t **values, 
            int this_dim, index_apply_f f, int *ctr, void *args){
    while (!isinf(icon[this_dim][indices[this_dim]].val)){
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
to an mnode_t* with a single mnode_t for each dimension.  */
void index_foreach(mnode_t *index[], index_apply_f f, void *args){
    int j, ctr=0;
    for (j=0; index[j]; ) j++;
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

This is for a contrast. We've already calculated the in/out vector for each value on
each margin, and now need to && each dimension into one vector.

\param index Is actually a partial index: for each dimension, there should be only one value. Useful for the center of an index_foreach loop.
\param d Should already be allocated to the right size, may be filled with garbage.
\return d will be zero or one to indicate which rows of the indexed data set meet all criteria.  */
void index_get_element_list(mnode_t *const * index, bool *d, size_t len, bool is_margin){
    memcpy(d, is_margin ? index[0]->margin_ptrs: index[0]->fit_ptrs, len *  sizeof(bool));
    for(size_t i=1; !isinf(index[i]->val); i++){
        bool *ptr_list = is_margin ? index[i]->margin_ptrs: index[i]->fit_ptrs;
        for(size_t j=0; j < len; j++)
             d[j] &= ptr_list[j];
    }
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
    double *maxdev;
} rake_t;

static void rakeinfo_grow(rake_t *r){
    r->al = (r->al+1)*2;
    r->elmtlist = realloc(r->elmtlist , sizeof(size_t*) * r->al);
    r->elmtlist_sizes = realloc(r->elmtlist_sizes, sizeof(size_t) * r->al);
    r->indata_values = apop_vector_realloc(r->indata_values, r->al);
}

static void rakeinfo_free(rake_t r){
    #define free_and_clear(in) free(in), (in) = NULL
    for (int i=0; i < r.ct; i++)
        free_and_clear(r.elmtlist[i]);
    free_and_clear(r.index); //these are just pointers to the main index.
    free_and_clear(r.elmtlist_sizes);
    free_and_clear(r.elmtlist);
    gsl_vector_free(r.indata_values);
    r.indata_values= NULL;
}

static void scaling(size_t const *elmts, size_t const n,  gsl_vector *weights, double const in_sum, double *maxdev){
    double fit_sum = 0, out_sum=0;
    for(size_t i=0; i < n; i ++)
        fit_sum += weights->data[elmts[i]];
    if (!fit_sum) return; //can happen if init table is very different from margins.
    for(size_t i=0; i < n; i ++){
        out_sum +=
        weights->data[elmts[i]] *= in_sum/fit_sum;
    }
    *maxdev = GSL_MAX(fabs(fit_sum - out_sum), *maxdev);
}

/* Given one set of values from one margin, do the actual scaling.
 On the first pass, this function takes notes on each margin's element list and total 
 in the original data. Later passes just read the notes and call the scaling() function above.
*/
static void one_set_of_values(mnode_t *const * const margincons, int ctr, void *in){
    rake_t *r = in;
    size_t marginsize = r->indata->matrix->size1;
    size_t fitsize = r->fit->matrix->size1;
    bool melmts[marginsize];
    bool fitelmts[fitsize];
	bool first_pass = false;
    double in_sum;
	if (ctr < r->ct)
		in_sum = gsl_vector_get(r->indata_values, ctr);
    else {
        r->ct++;
        if (ctr >= r->al || r->al==0) rakeinfo_grow(r);
   		index_get_element_list(margincons, melmts, marginsize, true);
   		index_get_element_list(margincons, fitelmts, fitsize, false);
        in_sum = 0;
        int n=0, al=0;
        r->elmtlist[ctr] = NULL;
        //use margin index to get total for this margin.
        for(int m=0; m < marginsize; m++)
            if (melmts[m]) in_sum += r->indata->weights->data[m];
        //use fit index to get elements involved in this margin
        for(int m=0; m < fitsize; m++)
            if (fitelmts[m]){
                if (n >= al) {
                    al = (al+1)*2;
                    r->elmtlist[ctr] = realloc(r->elmtlist[ctr], al*sizeof(size_t));
                }
                r->elmtlist[ctr][n++] = m;
            }
        r->elmtlist_sizes[ctr] = n;
        r->indata_values->data[ctr] = in_sum;
		first_pass = true;
	}
    if (!r->elmtlist_sizes[ctr]) return;
    if (!first_pass && !in_sum)  return;
    scaling(r->elmtlist[ctr], r->elmtlist_sizes[ctr], r->fit->weights, in_sum, r->maxdev);
}

/* For each configuration margin, for each combination for that margin, 
   call the above one_set_of_values() function. */
static void main_loop(int config_ct, rake_t *rakeinfo, int k){
    for (size_t i=0; i < config_ct; i ++)
		if (k==1) index_foreach(rakeinfo[i].index, one_set_of_values, rakeinfo+i);
		else
			for(int m=0; m < rakeinfo[i].ct; m++)
				one_set_of_values(NULL, m, rakeinfo+i);
}

/* Following the FORTRAN, 1 contrast ==> icon. Here, icon will be a
subset of the main index including only the columns pertaining to a given margin. */
void generate_margin_index(mnode_t **icon, const apop_data *margin, mnode_t **mainindex, size_t col){
    Apop_col_v(margin, col, iconv);
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
static void c_loglin(const apop_data *config, const apop_data *indata, 
                        apop_data *fit, double tolerance, int maxit) {
    mnode_t ** index = index_generate(indata, fit);

    /* Make a preliminary adjustment to obtain the fit to an empty configuration list */
    //fit->weights is either all 1 (for synthetic data) or the initial counts from the db.
    double x = apop_sum(indata->weights);
    double y = apop_sum(fit->weights);
    gsl_vector_scale(fit->weights, x/y);

	int contrast_ct =config && config->matrix ? config->matrix->size2 : 0;
    rake_t rakeinfos[contrast_ct];
    double maxdev=0;
    for(size_t i=0; i < contrast_ct; i ++){
        Apop_col_v(config, i, iconv)
        rakeinfos[i] = (rake_t) {
            .indata = indata, 
            .fit = fit, 
            .maxdev = &maxdev, //one value shared across dimensions
            .index = malloc(sizeof(mnode_t) *(apop_sum(iconv)+1)),
            //others are NULL, to be filled in as we go.
        };
        generate_margin_index(rakeinfos[i].index, config, index, i);
    }
    int k;
    gsl_vector *previous = apop_vector_copy(fit->weights);
    for (k = 1; k <= maxit; ++k) {
        maxdev = 0;
        main_loop(contrast_ct, rakeinfos, k);
        Apop_notify(3, "Data set after round %i of raking.\n", k);
        if (apop_opts.verbose >=3) apop_data_print(fit, .output_pipe=apop_opts.log_file);
        if (maxdev < tolerance) break;// Normal termination 
        gsl_vector_memcpy(previous, fit->weights);
    }
    cleanup(index, rakeinfos,contrast_ct);
    gsl_vector_free(previous);
    Apop_stopif(k == maxit, fit->error='c', 0, "Maximum number of iterations reached.");
}

char *pipe_parse = "[ \n\t]*([^| \n\t]+)[ \n\t]*([|]|$)";

apop_data **generate_list_of_contrasts(char *const *contras_in, int contrast_ct){
  apop_data** out = malloc(sizeof(apop_data*)* contrast_ct);
	for (int i=0; i< contrast_ct; i++) {
        apop_regex(contras_in[i], pipe_parse, out+i);
        apop_text_alloc(out[i], *out[i]->textsize, 1);
    }
	return out;
}

apop_data *get_var_list(char const *margin_table, char const *count_col, char const *init_count_col, char const *all_vars,
                        char * const *varlist, int *var_ct){
    apop_data *all_vars_d=NULL;
    if (!all_vars && !varlist){
        Apop_stopif(apop_opts.db_engine=='m', apop_return_data_error(y),
                    0, "I need a list of the full set of variable "
                       "names sent as .varlist= (char *[]){\"var1\", \"var2\",...}");
        //use SQLite's table_info, then shift the second col to the first.
        all_vars_d = apop_query_to_text("PRAGMA table_info(%s)", margin_table);
        int ctr=0;
        for (int i=0; i< all_vars_d->textsize[0]; i++)
            if (all_vars_d->text[i][1] && (count_col ? strcmp(all_vars_d->text[i][1], count_col) : 1)
                 && (init_count_col ? strcmp(all_vars_d->text[i][1], init_count_col): 1))
                    apop_text_add(all_vars_d, 0, ctr++, all_vars_d->text[i][1]);
        apop_text_alloc(all_vars_d, 1, ctr);
        *var_ct = ctr;
    } else if (!varlist){
        apop_regex(all_vars, pipe_parse, &all_vars_d);
        *var_ct = *all_vars_d->textsize;
        return apop_data_transpose(all_vars_d);
    } else {
        all_vars_d = apop_text_alloc(NULL, 1, *var_ct); 
        for (int i=0; i<*var_ct; i++) apop_text_add(all_vars_d, 0, i, varlist[i]);
    }
    Apop_stopif(!all_vars_d, apop_return_data_error(y), 0, "Trouble getting/parsing the list of variables.");
    return all_vars_d;
}

static int get_var_index(char *const *all_vars, int len, char *findme){
	for (int i=0; i< len; i++)
		if (all_vars[i] && !strcmp(all_vars[i], findme))
			return i;
	Apop_notify(0, "I couldn't find %s in the full list of variables. Returning -1.", findme);
    return -1;
}

void nan_to_zero(double *in){ if (gsl_isnan(*in)) *in=0;}

double nudge_zeros(apop_data *in, void *nudge){
    if (!in->weights->data[0])
        in->weights->data[0] = *(double*)nudge;
    return 0;
}

int find_in_allvars(char const *in, apop_data const *allvars){
    for (int i=0; i< allvars->textsize[1]; i++)
        if (!strcmp(allvars->text[0][i], in)) return i;
    Apop_stopif(1, return -1, 0, "Variable in your contrast list [%s] not in "
                                 "your list of all variables.", in);
}

/* If you are fully synthesizing or nudging zero cells, then calculate the set of 
   cells that could be nonzero (given the margin information). 
   * If using only an init_table, then that's your list of nonzero cells right there.
   * If you have an init_table but want to nudge the zero cells up a bit, then you need
     this, and have to merge with the init_table
   * If you have no init_table, then this list is all the cells.
 
 */
static int setup_nonzero_contrast(char const *margin_table, 
              apop_data const * allvars, int run_number,
              char const *list_of_fields, apop_data *const* contras, int contrast_ct,
              double nudge, char const * structural_zeros, bool have_init_table){
    char *q;
    bool used[allvars->textsize[1]];
    memset(used, 0, sizeof(bool)*allvars->textsize[1]);
	Asprintf(&q, "create table apop_zerocontrasts_%i as select %s, %g from\n", 
            run_number, apop_text_paste(allvars, .between=","), nudge);
    for (int i=0; i < contrast_ct; i++){
        xprintf(&q, "%s%s (select distinct %s  from %s) \n", 
                  q, i>0 ? "natural join" : "", 
                  apop_text_paste(contras[i], .between=","),  margin_table);
        for (int j=0; j<*contras[i]->textsize; j++){
            int val=find_in_allvars(*contras[i]->text[j], allvars);
            Apop_stopif(val==-1, return -1, 0, "Error setting up contrasts");
            used[val]=true; 
        }
    } 
    //make sure all variables are joined in. 
    for (int i=0; i < allvars->textsize[1]; i++)
        if (!used[i]) xprintf(&q, "%s%s (select distinct %s  from %s)\n", 
                  q, (!contrast_ct && i==0) ? "" : "natural join",
                  allvars->text[0][i],  margin_table);
    if (structural_zeros) xprintf(&q, "%s where not (%s)\n", q, structural_zeros);
    if (have_init_table){
        //Keep out margins with values for now; join them in below.
        xprintf(&q, "%s except\nselect %s, %g from %s", q, list_of_fields, nudge, margin_table);
    }
	Apop_stopif(apop_query("%s", q), return 1, 0, "query failed.");
    free(q);
    return 0;
}

/** Fit a log-linear model via iterative proportional fitting, aka raking.

Raking has many uses. The <a href="http://modelingwithdata.org/arch/00000138.htm">Modeling with Data blog</a> presents a series of discussions 
of uses of raking, including some worked examples.

Or see Wikipedia for an overview of Log linear models, aka
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
points, not to the full cross of all categories. Set <tt>apop_opts.verbose</tt> to 2 or greater to show the query on \c stderr.

\li One could use raking to generate `fully synthetic' data: start with observation-level data in a margin table. Begin the raking with a starting data set of all-ones. Then rake until the all-ones set transforms into something that conforms to the margins and (if any) structural zeros. You now have a data set which matches the marginal totals but does not use any other information from the observation-level data. If you do not specify an <tt>.init_table</tt>, then an all-ones default table will be used.

\li Set <tt>apop_opts.verbose=3</tt> to see the intermediate tables at the end of each round of raking.

\li If you want all cells to have nonzero value, then you can do that via pre-processing:
\code
apop_query("update data_table set count_col = 1e-3 where count_col = 0");
\endcode


\param margin_table The name of the table in the database to use for calculating
the margins.  The table should have one observation per row.  No default. (This used
to be called \c table_name; that name is now deprecated.)

\param var_list The full list of variables to search. A list of strings, e.g., <tt>(char *[]){"var1", "var2", ..., "var15"}</tt>

\param var_ct The count of the full list of variables to search.

\param all_vars deprecated.

\param contrasts The contrasts describing your model. Like the \c all_vars input, each
contrast is a pipe-delimited list of variable names. No default.

\param contrast_ct The number of contrasts in the list of contrasts. No default.

\param structural_zeros a SQL clause indicating combinations that can never take a nonzero
value. This will go into a \c where clause, so anything you could put there is OK, e.g.
"age <65 and gets_soc_security=1 or age <15 and married=1". 
Your margin data is not checked for structural zeros.  Default: no structural zeros.

\param max_iterations Number of rounds of raking at which the algorithm halts. Default: 1000.

\param tolerance I calculate the change for each cell from round to round;
if the largest cell change is smaller than this, I stop. Default: 1e-5.

\param count_col This column gives the count of how many observations are represented
by each row. If \c NULL, ech row represents one person. Default: \c NULL.

\param run_number Because I write intermediate tables to the database, I need a way to
distinguish distinct runs should you be threading several runs at once. If you aren't
running several instances simultaneously, don't worry about this; if you are, do supply a
value, since it's hard for the function to supply one in a race-proof manner. Default:
internally-maintained values.

\param init_table The default is to initially set all table elements to one and then
rake from there. This is effectively the `fully synthetic' approach, which uses only
the information in the margins and derives the data set closest to the all-ones data
set that is consistent with the margins. Care is taken to maintan sparsity in this
case.  If you specify an \c init_table, then I will get the initial cell counts from
it. Default: the fully-synthetic approach, using a starting point of an all-ones grid.

\param init_count_col The column in \c init_table with the cell counts.

\param nudge There is a common hack of adding a small value to every zero entry, because
a zero entry will always scale to zero, while a small value could eventually scale
to anything.  Recall that this function works on sparse sets, so I first filter out
those cells that could possibly have a nonzero value given the observations, then I
add <tt>nudge</tt> to any zero cells within that subset.

\param table_name Deprecated; replaced with \c margin_table.

\return An \ref apop_data set where every row is a single combination of variable values
and the \c weights vector gives the most likely value for each cell.
\exception out->error='i' Input was somehow wrong.
\exception out->error='c' Raking did not converge, reached max. iteration count.
\li This function uses the \ref designated syntax for inputs.

\li The interface is still beta, and subject to change---notably, handling of text
categories will soon be added.
*/
#ifdef APOP_NO_VARIADIC
apop_data * apop_rake(char const *margin_table, char * const*var_list, int var_ct, char const *all_vars, char * const *contrasts, int contrast_ct, char const *structural_zeros, int max_iterations, double tolerance, char const *count_col, int run_number, char const *init_table, char const *init_count_col, double nudge, char const* table_name){
#else
apop_varad_head(apop_data *, apop_rake){
    #if __STDC_VERSION__ > 201100L && !defined(__STDC_NO_ATOMICS__)
        static _Atomic(int) defaultrun = 0;
    #else
        static int defaultrun = 0;
    #endif
    char const * apop_varad_var(table_name, NULL); //the deprecated name for margin_table
    char const * apop_varad_var(margin_table, table_name);
    Apop_stopif(!margin_table, apop_return_data_error(i), 0,  
                        "I need the name of a table in the database that will be the data source.");
    Apop_stopif(!apop_table_exists(margin_table), apop_return_data_error(i), 
                        0, "your margin_table, %s, doesn't exist in the database.", margin_table);
    char *const* apop_varad_var(var_list, NULL);
    int apop_varad_var(var_ct, 0);
    char const * apop_varad_var(all_vars, NULL); //deprecated; use var_list/var_ct
    char *const * apop_varad_var(contrasts, NULL); //default to all vars?
    int apop_varad_var(contrast_ct, 0);
    Apop_stopif(contrasts&&!contrast_ct, apop_return_data_error(i),
            0, "you gave me a list of contrasts but not the count. "
            "This is C--I can't count them myself. Please provide the count and re-run.");
    char const * apop_varad_var(structural_zeros, NULL);
    char const * apop_varad_var(count_col, NULL);
    int apop_varad_var(max_iterations, 1e3);
    double apop_varad_var(tolerance, 1e-5);
    int apop_varad_var(run_number, defaultrun++);
    char const * apop_varad_var(init_count_col, NULL);
    char const * apop_varad_var(init_table, NULL);
    Apop_stopif(init_table && !apop_table_exists(init_table), apop_return_data_error(i),
               0, "your init_table, %s, doesn't exist in the database.", init_table);
    if (init_count_col && !init_table) init_table = margin_table;
    double apop_varad_var(nudge, 0);
    return apop_rake_base(margin_table, var_list, var_ct, all_vars, contrasts, contrast_ct, structural_zeros, max_iterations, tolerance, count_col, run_number, init_table, init_count_col, nudge, table_name);
}

 apop_data * apop_rake_base(char const *margin_table, char * const*var_list, int var_ct, char const *all_vars, char * const *contrasts, int contrast_ct, char const *structural_zeros, int max_iterations, double tolerance, char const *count_col, int run_number, char const *init_table, char const *init_count_col, double nudge, char const* table_name){
#endif
	apop_data **contras = generate_list_of_contrasts(contrasts, contrast_ct);
    apop_data *all_vars_d = get_var_list(margin_table, count_col, init_count_col, all_vars, var_list, &var_ct);
    Apop_stopif(all_vars_d->error, return all_vars_d, 0, "Trouble setting up the list of variables.");
    int tt = all_vars_d->textsize[0]; all_vars_d->textsize[0] = 1; //mask all but the first row
    char *list_of_fields = apop_text_paste(all_vars_d, .between=", ");

    if (nudge || !init_table){
        char *tab;
        Asprintf(&tab, "apop_zerocontrasts_%i", run_number);
        apop_table_exists(tab, 'd');
        free(tab);
        Apop_stopif(setup_nonzero_contrast(margin_table, all_vars_d,  
                        run_number, list_of_fields, contras, contrast_ct, (nudge ? nudge : 1), structural_zeros, init_table),
                 apop_return_data_error(q),
                 0, "Couldn't calculate the set of nonzero cells.");
    }

    char *initt=NULL; //handle structural zeros via subquery
    //note that margin data may have invalid rows.
    if (init_table){
        if (structural_zeros) 
             Asprintf(&initt, "(select * from %s where not (%s))", init_table, structural_zeros);
        else Asprintf(&initt, "%s", init_table);
    }
    char *init_q, *pre_init_q = NULL;
    if (init_table){
        char *countstr;
        if (init_count_col) Asprintf(&countstr, "sum(%s) as %s", init_count_col, init_count_col);
        else                Asprintf(&countstr, "count(%s)", **all_vars_d->text);
        Asprintf(&init_q, "select %s, %s from %s group by %s", 
                           list_of_fields, countstr, initt, list_of_fields);
        free(countstr);
    }

    char *marginq, *cc; 
    if (count_col) Asprintf(&cc, "sum(%s)", count_col);
    else           cc = strdup("count(*)");
    Asprintf(&marginq, "select %s, %s  from %s\ngroup by %s", 
                       list_of_fields, cc, margin_table, list_of_fields);
    free(cc); free(initt);

    char *format=strdup("w");
    for (int i =0 ; i< var_ct; i++)
        xprintf(&format, "m%s", format);
    apop_data *d, *contrast_grid;
    d = apop_query_to_mixed_data(format, "%s", marginq);
    Apop_stopif(!d || d->error, apop_return_data_error(q),
            0, "This query:\n%s\ngenerated a blank or broken table.", marginq);
    free(marginq);

    if (pre_init_q) Apop_stopif(apop_query("%s", pre_init_q), apop_return_data_error(q),
            0, "This query:\n%s\ngenerated a blank or broken table.", pre_init_q);

    apop_data *fit;
    if (init_table) {
        fit = (nudge) 
               ? apop_query_to_mixed_data(format, "%s\nunion\nselect * from apop_zerocontrasts_%i ", init_q, run_number)
               : apop_query_to_mixed_data(format, "%s", init_q);
        Apop_stopif(!fit, apop_return_data_error(q), 0, "Query returned a blank table.");
        Apop_stopif(fit->error, apop_return_data_error(q), 0, "Query error.");
    } else {
        fit = apop_query_to_mixed_data(format, "select * from apop_zerocontrasts_%i ", run_number);
        gsl_vector_set_all(fit->weights, nudge ? nudge : 1);
    }
    free(format);
    apop_vector_apply(fit->weights, nan_to_zero);
    if (nudge) apop_map(fit, .fn_rp=nudge_zeros, .param=&nudge);

    contrast_grid = apop_data_calloc(var_ct, contrast_ct);
	for (int i=0; i< contrast_ct; i++)
		for (int j=0; j< contras[i]->textsize[0]; j++)
			apop_data_set(contrast_grid, get_var_index(*all_vars_d->text, all_vars_d->textsize[1], contras[i]->text[j][0]), i, 1);
	
    if (!init_table || nudge)
        for (int i=0; i< contrast_ct; i++) apop_data_free(contras[i]);

    c_loglin(contrast_grid, d, fit, tolerance, max_iterations);
    apop_data_free(d);
    
    all_vars_d->textsize[0] = tt;
    apop_data_free(all_vars_d);
    if (!init_table || nudge) apop_query("drop table apop_zerocontrasts_%i", run_number);
	apop_data_free(contrast_grid);
	return fit;
}
