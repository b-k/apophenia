
/** \file apop_regression.c	Generally, if it assumes something is  Normally distributed, it's here.*/
/* Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/** \defgroup regression  OLS/GLS: The linear projection methods */
/** \defgroup ttest  T-tests: comparing two vectors */
/** \defgroup asst_tests  Various means of hypothesis testing.

 See also the goodness of fit tests in \ref histograms.
 */

#include "apop_internal.h"
#include <search.h> //lsearch; bsearch is in stdlib.

/** For many, it is a knee-jerk reaction to a parameter estimation to test whether each individual parameter differs from zero. This function does that.

\param est  The \ref apop_model, which includes pre-calculated parameter estimates, var-covar matrix, and the original data set.

Returns nothing. At the end of the routine, <tt>est->info->more</tt> includes a set of t-test values: p value, confidence (=1-pval), t statistic, standard deviation, one-tailed Pval, one-tailed confidence.

*/
void apop_estimate_parameter_tests (apop_model *est){
    Nullcheck_p(est, )
    if (!est->data) return;
    apop_data *ep = apop_data_add_page(est->info, apop_data_alloc(est->parameters->vector->size, 2), "<test info>");
    apop_name_add(ep->names, "p value", 'c');
    apop_name_add(ep->names, "confidence", 'c');
    apop_name_stack(ep->names, est->parameters->names, 'r', 'r');
    Get_vmsizes(est->data); //msize1, vsize
    int df = msize1 ? msize1 : vsize;
    df -= est->parameters->vector->size;
    df  = df < 1 ? 1 : df; //some models aren't data-oriented.
    apop_data_add_named_elmt(est->info, "df", df);

    apop_data *one_elmt = apop_data_calloc(1, 1);
    gsl_vector *param_v = apop_data_pack(est->parameters);
    for (size_t i=0; i< est->parameters->vector->size; i++){
        Apop_settings_add_group(est, apop_pm, .index=i);
        apop_model *m = apop_parameter_model(est->data, est);

        double zero = apop_cdf(one_elmt, m);
        apop_model_free(m);
        double conf = 2*fabs(0.5-zero); //parameter is always at 0.5 along a symmetric CDF
        apop_data_set(ep, i, .colname="confidence", .val=conf);
        apop_data_set(ep, i, .colname="p value",    .val=1-conf);
    }
    gsl_vector_free(param_v);
    apop_data_free(one_elmt);
}

//Cut and pasted from the GNU std library documentation, modified to consider NaNs:
static int compare_doubles (const void *a, const void *b) {
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    if (gsl_isnan(*da)) return gsl_isnan(*db) ? 0 : 1;
    if (gsl_isnan(*db)) return -1;
    return (*da > *db) - (*da < *db);
}

typedef const char * ccp;
static int strcmpwrap(const void *a, const void *b){
    const ccp *aa = a;
    const ccp *bb = b;
    return strcmp(*aa, *bb);
}

/** Give me a vector of numbers, and I'll give you a sorted list of the unique elements. 
  This is basically running "select distinct datacol from data order by datacol", but without the aid of the database.

  \param v a vector of items
  \return a sorted vector of the distinct elements that appear in the input.
  \li NaNs appear at the end of the sort order.
  \see apop_text_unique_elements 
*/
gsl_vector * apop_vector_unique_elements(const gsl_vector *v){
    size_t prior_elmt_ctr = 107;
    size_t elmt_ctr = 0;
    double *elmts = NULL;
    for (size_t i=0; i< v->size; i++){
        if (prior_elmt_ctr != elmt_ctr)
            elmts = realloc(elmts, sizeof(double)*(elmt_ctr+1));
        prior_elmt_ctr = elmt_ctr;
        double val = gsl_vector_get(v, i);
        lsearch(&val, elmts, &elmt_ctr, sizeof(double), compare_doubles);
        if (prior_elmt_ctr < elmt_ctr)
            qsort(elmts, elmt_ctr, sizeof(double), compare_doubles);
    }
    gsl_vector *out = apop_array_to_vector(elmts, elmt_ctr);
    free(elmts);
    return out;
}

/** Give me a column of text, and I'll give you a sorted list of the unique
  elements. 
  This is basically running "select distinct * from datacolumn", but without 
  the aid of the database.  

  \param d An \ref apop_data set with a text component
  \param col The text column you want me to use.
  \return An \ref apop_data set with a single sorted column of text, where each unique text input appears once.
  \see apop_vector_unique_elements
*/
apop_data * apop_text_unique_elements(const apop_data *d, size_t col){
  char   **tval;

  //first element for free
  size_t prior_elmt_ctr, elmt_ctr = 1;
  char **telmts = malloc(sizeof(char**)*2);
  telmts[0] = d->text[0][col];

    for (int i=1; i< d->textsize[0]; i++){
        prior_elmt_ctr  = elmt_ctr;
        tval    =  &(d->text[i][col]);
        lsearch (tval, telmts, &elmt_ctr, sizeof(char*), strcmpwrap);
        if (prior_elmt_ctr  < elmt_ctr){
            qsort(telmts, elmt_ctr, sizeof(char*), strcmpwrap);
            telmts = realloc(telmts, sizeof(char**)*(elmt_ctr+1));
        }
    }

    //pack and ship
    apop_data *out = apop_text_alloc(NULL, elmt_ctr, 1);
    for (int j=0; j< elmt_ctr; j++)
        apop_text_add(out, j, 0, telmts[j]);
    free(telmts);
    return out;
}

static char *apop_get_factor_basename(apop_data *d, int col, char type){
    char *name;
    char *catname =   d->names == NULL ? NULL
                    : type == 't' && d->names && d->names->textct > col ? d->names->text[col]
                    : col == -1 && d->names && d->names->vector         ? d->names->vector
                    : col >=0 && d->names && d->names->colct > col      ? d->names->col[col]
                    : NULL;
    if (catname){
        Asprintf(&name, "%s", catname);
        return name;
    }
    if (type == 't'){
        Asprintf(&name, "text column %i", col);
        return name;
    }
    if (col == -1)  Asprintf(&name, "vector");
    else            Asprintf(&name, "column %i", col);
    return name;
}

static char *make_catname (apop_data *d, int col, char type){
    char *name, *subname = apop_get_factor_basename(d, col, type);
    Asprintf(&name, "<categories for %s>", subname);
    free(subname);
    return name;
}

/* Producing dummies consists of finding the index of element i, for all i, then
 setting (i, index) to one.
 Producing factors consists of finding the index and then setting (i, datacol) to index.
 Otherwise the work is basically identical.  
 Also, add a ->more page to the input data giving the translation.
 */
static apop_data * dummies_and_factors_core(apop_data *d, int col, char type, 
                            int keep_first, int datacol, char dummyfactor, 
                            apop_data **factor_list){
    size_t index, elmt_ctr = 0;
    gsl_vector *delmts = NULL;
    char **telmts = NULL;//unfortunately needed for the bsearch.

    //first, create an ordered list of unique elements.
    //Record that list for use in this function, and in a ->more page of the data set.
    char *catname =  make_catname(d, col, type);
    if (type == 't'){
        *factor_list = apop_data_add_page(d, apop_text_unique_elements(d, col), catname);
        elmt_ctr = (*factor_list)->textsize[0];
        //awkward format conversion:
        telmts = malloc(sizeof(char*)*elmt_ctr);
        for (size_t j=0; j< elmt_ctr; j++)
            Asprintf(&(telmts[j]), "%s", (*factor_list)->text[j][0]);
        (*factor_list)->vector = gsl_vector_alloc(elmt_ctr);
        for (size_t i=0; i< (*factor_list)->vector->size; i++)
            apop_data_set(*factor_list, i, -1, i);
    } else {
        if (col==-1)
            delmts = apop_vector_unique_elements(d->vector);
        else{
            Apop_col_v(d, col, to_search);
            delmts = apop_vector_unique_elements(to_search);
        }
        elmt_ctr = delmts->size;
        *factor_list = apop_data_add_page(d, apop_data_alloc(elmt_ctr), catname);
        apop_text_alloc((*factor_list), delmts->size, 1);
        for (size_t i=0; i< (*factor_list)->vector->size; i++){
            //shift to the text, for conformity with the more common text version.
            apop_text_add((*factor_list), i, 0, "%g", gsl_vector_get(delmts, i));
            apop_data_set((*factor_list), i, -1, i);
        }
    }

    //Now go through the input vector, and for row i find the posn of the vector's
    //name in the element list created above (j), then change (i,j) in
    //the dummy matrix to one.
    int s = type == 't' 
            ? d->textsize[0]
            : (col >=0 ? d->matrix->size1 : d->vector->size);
    apop_data *out = (dummyfactor == 'd')
                ? apop_data_calloc(0, s, (keep_first ? elmt_ctr : elmt_ctr-1))
                : d;
    for (size_t i=0; i< s; i++){
        if (type == 'd'){
            double val = apop_data_get(d, i, col);
            index = ((size_t)bsearch(&val, delmts->data, elmt_ctr, sizeof(double), compare_doubles) - (size_t)delmts->data)/sizeof(double);
        } else 
            index   = ((size_t)bsearch(&(d->text[i][col]), telmts, elmt_ctr, sizeof(char**), strcmpwrap) - (size_t)telmts)/sizeof(char**);
        if (dummyfactor == 'd'){
            if (keep_first)
                gsl_matrix_set(out->matrix, i, index,1); 
            else if (index > 0)   //else don't keep first and index==0; throw it out. 
                gsl_matrix_set(out->matrix, i, index-1, 1); 
        } else
            apop_data_set(out, i, datacol, index); 
    }
    //Add names:
    if (dummyfactor == 'd'){
        char *basename = apop_get_factor_basename(d, col, type);
        for (size_t i = (keep_first) ? 0 : 1; i< elmt_ctr; i++){
            char n[1000];
            if (type =='d'){
                sprintf(n, "%s dummy %g", basename, gsl_vector_get(delmts,i));
            } else
                sprintf(n, "%s", telmts[i]);
            apop_name_add(out->names, n, 'c');
        }
    }
    if (delmts)
        gsl_vector_free(delmts);
    if (telmts){
        for (size_t j=0; j< elmt_ctr; j++)
            free(telmts[j]);
        free(telmts);
    }
    free(catname);
    return out;
}

/** A utility to make a matrix of dummy variables. You give me a single
vector that lists the category number for each item, and I'll produce
a matrix with a single one in each row in the column specified.

After that, you have to decide what to do with the new matrix and the original data column. 

\li You can manually join the dummy data set with your main data, e.g.:
\code
apop_data *dummies  = apop_data_to_dummies(main_regression_vars, .col=8, .type='t');
apop_data_stack(main_regression_vars, dummies, 'c', .inplace='y');
\endcode

\li The <tt>.remove='y'</tt> option specifies that I should use \ref apop_data_rm_columns 
to remove the column used to generate the dummies. Implemented only for <tt>type=='d'</tt>.

\li By specifying <tt>.append='y'</tt> or <tt>.append='e'</tt> I will run the above two lines for you. Your \ref apop_data pointer will not change, but its \c matrix element will be reallocated (via \ref apop_data_stack).

\li By specifying <tt>.append='i'</tt>, I will place the matrix of dummies in place,
immediately after the data column you had specified. You will probably use this with
<tt>.remove='y'</tt> to replace the single column with the new set of dummy columns.
Bear in mind that if there are two or more dummy columns (which there probably are if you
are bothering to use this function), subsequent column numbers will change.

\li If <tt>.append='i'</tt> and you asked for a text column, I will append to the end of
the table, which is equivalent to <tt>append='e'</tt>.

\param  d The data set with the column to be dummified (No default.)
\param col The column number to be transformed; -1==vector (default = 0)
\param type 'd'==data column, 't'==text column. (default = 't')
\param  keep_first  if zero, return a matrix where each row has a one in the (column specified MINUS
    ONE). That is, the zeroth category is dropped, the first category
    has an entry in column zero, et cetera. If you don't know why this
    is useful, then this is what you need. If you know what you're doing
    and need something special, set this to one and the first category won't be dropped. (default = 0)
\param append If \c 'e' or \c 'y', append the dummy grid to the end of the original data
matrix. If \c 'i', insert in place, immediately after the original data column. (default = \c 'n')
\param remove If \c 'y', remove the original data or text column. (default = \c 'n')

\return An \ref apop_data set whose \c matrix element is the one-zero
matrix of dummies. If you used <tt>.append</tt>, then this is the main matrix.
Also, I add a page named <tt>"\<categories for your_var\>"</tt> giving a reference table of names and column numbers (where <tt>your_var</tt> is the appropriate column heading).
\exception out->error=='a' allocation error
\exception out->error=='d' dimension error
\li NaNs appear at the end of the sort order.
\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
apop_data * apop_data_to_dummies(apop_data *d, int col, char type, int keep_first, char append, char remove){
#else
apop_varad_head(apop_data *, apop_data_to_dummies){
    apop_data *apop_varad_var(d, NULL)
    Apop_stopif(!d, return NULL, 1, "You sent me a NULL data set for apop_data_to_dummies. Returning NULL.");
    int apop_varad_var(col, 0)
    char apop_varad_var(type, 't')
    int apop_varad_var(keep_first, 0)
    char apop_varad_var(append, 'n')
    char apop_varad_var(remove, 'n')
    if (remove =='y' && type == 't') Apop_notify(1, "Remove isn't implemented for text source columns yet.");
    return apop_data_to_dummies_base(d, col, type, keep_first, append, remove);
}

 apop_data * apop_data_to_dummies_base(apop_data *d, int col, char type, int keep_first, char append, char remove){
#endif
    if (type == 'd'){
        Apop_stopif((col == -1) && d->vector, apop_return_data_error(d),
                                0, "You asked for the vector element "
                                "(col==-1) but the data's vector element is NULL.");
        Apop_stopif((col != -1) && (col >= d->matrix->size2), apop_return_data_error(d),
                                0, "You asked for the matrix element %i "
                                "but the data's matrix element has only %zu columns.", col, d->matrix->size2);
    } else Apop_stopif(col >= d->textsize[1], apop_return_data_error(d),
                                0, "You asked for the text element %i but "
                                    "the data's text element has only %zu elements.", col, d->textsize[1]);
    apop_data *fdummy;
    apop_data *dummies= dummies_and_factors_core(d, col, type, keep_first, 0, 'd', &fdummy);
    //Now process the append and remove options.
    size_t orig_size = d->matrix ? d->matrix->size1 : 0;
    int rm_list[orig_size+1];
    memset (rm_list, 0, (orig_size+1)*sizeof(int)); 
    if (append =='i'){
        apop_data **split = apop_data_split(d, col+1, 'c');
        //stack names, then matrices
        for (int i=0; i < d->names->colct; i++)
            free(d->names->col[i]);
        apop_name_stack(d->names, split[0]->names, 'c');
        for (int k = d->names->colct; k < (split[0]->matrix ? split[0]->matrix->size2 : 0); k++)
            apop_name_add(d->names, "", 'c'); //pad so the name stacking is aligned (if needed)
        apop_name_stack(d->names, dummies->names, 'c');
        apop_name_stack(d->names, split[1]->names, 'c');
        gsl_matrix_free(d->matrix);
        d->matrix = apop_matrix_stack(split[0]->matrix, dummies->matrix, 'c');
        apop_data_free(dummies);
        apop_data_free(split[0]);
        apop_matrix_stack(d->matrix, split[1]->matrix, 'c', .inplace='y');
        apop_data_free(split[1]);
        return d;
    }
    if (remove!='n' && type!='t'){
        rm_list[col]=1;
        apop_data_rm_columns(d, rm_list);
    }
    if (append =='y' || append == 'e' || append ==1 || (append=='i' && type=='t')){
        d = apop_data_stack(d, dummies, 'c', .inplace='y');
        apop_data_free(dummies);
        return d;
    }
    return dummies;
}

/** Convert a column of text or numbers
  into a column of numeric factors, which you can use for a multinomial probit/logit, for example.

  If you don't run this on your data first, \ref apop_probit and \ref apop_logit default to running 
  it on the vector or (if no vector) zeroth column of the matrix of the input \ref apop_data set, because those models need a list of the unique values of the dependent variable.

\param data The data set to be modified in place. (No default. If \c NULL, returns \c NULL and a warning)
\param intype If \c 't', then \c incol refers to text, otherwise (\c 'd'
is a good choice) refers to the vector or matrix. Default = \c 't'.
\param incol The column in the text that will be converted. -1 is the vector. Default = 0.
\param outcol The column in the data set where the numeric factors will be written (-1 means the vector). Default = 0.

For example:
\code
apop_data *d  = apop_query_to_mixed_data("mmt", "select 1, year, color from data");
apop_data_to_factors(d);
\endcode
Notice that the query pulled a column of ones for the sake of saving room for the factors. It reads column zero of the text, and writes it to column zero of the matrix.

Another example:
\code
apop_data *d  = apop_query_to_data("mmt", "select type, year from data");
apop_data_to_factors(d, .intype='d', .incol=0, .outcol=0);
\endcode
Here, the \c type column is converted to sequential integer factors and
those factors overwrite the original data. Since a reference table is
added as a second page of the \ref apop_data set, you can recover the
original values as needed.

\return A table of the factors used in the code. This is an \c apop_data set with only one column of text.
Also, I add a page named <tt>"<categories for your_var>"</tt> giving a reference table of names and column numbers (where <tt>your_var</tt> is the appropriate column heading) use \ref apop_data_get_factor_names to retrieve that table.

\exception out->error=='a' allocation error.
\exception out->error=='d' dimension error.
\li  If the vector or matrix you wanted to write to is \c NULL, I will allocate it for you.
\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
apop_data * apop_data_to_factors(apop_data *data, char intype, int incol, int outcol){
#else
apop_varad_head(apop_data *, apop_data_to_factors){
    apop_data *apop_varad_var(data, NULL)
    Apop_stopif(!data, return NULL, 1, "You sent me a NULL data set. Returning NULL.");
    int apop_varad_var(incol, 0)
    int apop_varad_var(outcol, 0)
    char apop_varad_var(intype, 't')
    return apop_data_to_factors_base(data, intype, incol, outcol);
}

 apop_data * apop_data_to_factors_base(apop_data *data, char intype, int incol, int outcol){
#endif
    if (intype=='t'){
        Apop_stopif(incol >= data->textsize[1], apop_return_data_error(d),
                        0, "You asked for the text column %i but the "
                           "data's text has only %zu elements.", incol, data->textsize[1]);
    } else {
        Apop_stopif((incol == -1) && !data->vector, apop_return_data_error(d),
                0, "You asked for the vector of the data set but there is none.");
        Apop_stopif((incol != -1) && !data->matrix, apop_return_data_error(d),
                0, "You asked for the matrix column %i but the matrix is NULL.", incol);
        Apop_stopif((incol != -1) && (incol >= data->matrix->size2), apop_return_data_error(d),
                0, "You asked for the matrix column %i but "
                            "the matrix has only %zu elements.", incol, data->matrix->size2);
    }
    if (!data->vector && outcol == -1) //allocate a vector for the user.
        data->vector = gsl_vector_alloc(intype=='t' ? data->textsize[0] : data->matrix->size2);
    if (!data->matrix && outcol >= 0) //allocate a matrxi for the user.
        data->matrix = gsl_matrix_calloc(intype=='t' ? data->textsize[0] : data->matrix->size2, outcol+1);
    apop_data *out;
    dummies_and_factors_core(data, incol, intype, 1, outcol, 'f', &out);
    return out;
}

/** Factor names are stored in an auxiliary table with a name like 
<tt>"<categories for your_var>"</tt>. Producing this name is annoying (and prevents us from eventually making it human-language independent), so use this function to get the list of factor names.

\param data The data set. (No default, must not be \c NULL)
\param col The column in the main data set whose name I'll use to check for the factor name list. Vector==-1. (default=0)
\param type If you are referring to a text column, use 't'. (default='d')

\return A pointer to the page in the data set with the given factor names.

\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
apop_data * apop_data_get_factor_names(apop_data *data, int col, char type){
#else
apop_varad_head(apop_data *, apop_data_get_factor_names){
    apop_data *apop_varad_var(data, NULL)
    Apop_stopif(!data, return NULL, 1, "You sent me a NULL data set. Returning NULL.");
    int apop_varad_var(col, 0)
    char apop_varad_var(type, 'd')
    return apop_data_get_factor_names_base(data, col, type);
}

 apop_data * apop_data_get_factor_names_base(apop_data *data, int col, char type){
#endif
    char *name = make_catname (data, col, type);
    apop_data *out = apop_data_get_page(data, name, .match='e');
    free(name);
    return out;
}


/** Deprecated. Use \ref apop_data_to_factors.
  
  Convert a column of text in the text portion of an \c apop_data set
  into a column of numeric elements, which you can use for a multinomial probit, for example.

\param d The data set to be modified in place.
\param datacol The column in the data set where the numeric factors will be written (-1 means the vector, which I will allocate for you if it is \c NULL)
\param textcol The column in the text that will be converted.

For example:
\code
apop_data *d  = apop_query_to_mixed_data("mmt", "select 1, year, color from data");
apop_text_to_factors(d, 0, 0);
\endcode
Notice that the query pulled a column of ones for the sake of saving room for the factors.

\return A table of the factors used in the code. This is an \c apop_data set with only one column of text.
Also, the <tt>more</tt> element is a reference table of names and column numbers.

\exception out->error=='d'  dimension error.
*/
apop_data *apop_text_to_factors(apop_data *d, size_t textcol, int datacol){
    Apop_stopif(textcol >= d->textsize[1],  apop_return_data_error(d),
                0, "You asked for the text element %i but the data's "
                   "text has only %zu elements.", datacol, d->textsize[1]);
    if (!d->vector && datacol == -1) //allocate a vector for the user.
        d->vector = gsl_vector_alloc(d->textsize[0]);
    apop_data *out;
    dummies_and_factors_core(d, textcol, 't', 1, datacol, 'f', &out);
    return out;
}

/** Good ol' \f$R^2\f$.  Let \f$Y\f$ be the dependent variable,
\f$\epsilon\f$ the residual,  \f$n\f$ the number of data points, and \f$k\f$ the number of independent vars (including the constant). Returns an \ref apop_data set with the following entries (in the vector element):

\li  \f$ SST \equiv \sum (Y_i - \bar Y) ^2 \f$
\li  \f$ SSE \equiv \sum \epsilon ^2       \f$
\li  \f$ R^2 \equiv 1 - {SSE\over SST}     \f$
\li  \f$ R^2_{adj} \equiv R^2 - {(k-1)\over (n-k-1)}(1-R^2)     \f$

  Internally allocates (and frees) a vector the size of your data set.

\return: a \f$5 \times 1\f$ apop_data table with the following fields:
\li "R squared"
\li "R squared adj"
\li "SSE"
\li "SST"
\li "SSR"

If the output is in \c sss, use <code>apop_data_get(sss, .rowname="SSE")</code> to get the SSE, and so on for the other items.

\param m    A model. I use the pointer to the data set used for estimation and the info page named \c "<Predicted>". 
The Predicted page should include observed, expected, and residual columns, which I use to
generate the sums of squared errors and residuals, et cetera. All generalized linear
models produce a page with this name and of this form, as do a host of other models. Nothing 
keeps you from finding the \f$R^2\f$ of, say, a kernel smooth; it is up to you to determine 
whether such a thing is appropriate to your given models and situation.

\li <tt>apop_estimate(yourdata, apop_ols)</tt> does this automatically
\li If I don't find a Predicted page, I throw an error on the screen and return \c NULL.
\li The number of observations equals the number of rows in the Predicted page
\li The number of independent variables, needed only for the adjusted \f$R^2\f$, is from the
number of columns in the main data set's matrix (i.e. the first page; i.e. the set of
parameters if this is the \c parameters output from a model estimation). 
\li If your data (first page again) has a \c weights vector, I will find weighted SSE,
SST, and SSR (and calculate the \f$R^2\f$s using those values).

\ingroup regression
  */
apop_data *apop_estimate_coefficient_of_determination (apop_model *m){
  double          sse, sst, rsq, adjustment;
  size_t          indep_ct= m->data->matrix->size2 - 1;
  apop_data       *out    = apop_data_alloc();
    gsl_vector *weights = m->data->weights; //typically NULL.
    apop_data *expected = apop_data_get_page(m->info, "<Predicted>");
    Apop_stopif(!expected, return NULL, 0, "I couldn't find a \"<Predicted>\" page in your data set. Returning NULL.\n");
    size_t obs = expected->matrix->size1;
    Apop_col_tv(expected, "residual", v)
    if (!weights)
        gsl_blas_ddot(v, v, &sse);
    else {
        gsl_vector *v_times_w = apop_vector_copy(weights);
        gsl_vector_mul(v_times_w, v);
        gsl_blas_ddot(v_times_w, v, &sse);
        gsl_vector_free(v_times_w);
    }
    Apop_col_v(expected, 0, vv);
    sst = apop_vector_var(vv, m->data->weights) * (vv->size-1);
    rsq = 1. - (sse/sst);
    adjustment  = ((obs -1.) /(obs - indep_ct)) * (1.-rsq) ;
    apop_data_add_named_elmt(out, "R squared", rsq);
    apop_data_add_named_elmt(out, "R squared adj", 1 - adjustment);
    apop_data_add_named_elmt(out, "SSE", sse);
    apop_data_add_named_elmt(out, "SST", sst);
    apop_data_add_named_elmt(out, "SSR", sst - sse);
    return out;
}

/**  \def apop_estimate_r_squared(in) 
 A synonym for \ref apop_estimate_coefficient_of_determination, q.v. 
 \hideinitializer
 \ingroup regression
 */
