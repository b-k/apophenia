/** \file apop_missing_data.c Some missing data handlers. */
/* Copyright (c) 2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "internal.h"
#include "variadic.h"
#include "likelihoods.h"

static double find_nans(double in){ return isnan(in); }

static void addin(apop_data *predict, size_t i, int j, size_t page){
    int len;
    if (!predict->matrix){
        predict->matrix = gsl_matrix_alloc(1,5); 
        len = 0;
    } else 
        len = predict->matrix->size1;
    apop_matrix_realloc(predict->matrix, len + 1, predict->matrix->size2);
    apop_data_set(predict, .row=len, .colname="row", .val=i);
    apop_data_set(predict, .row=len, .colname="col", .val=j);
    apop_data_set(predict, .row=len, .colname="page", .val=page);
}

static int find_missing(const apop_data *data, apop_data *predict, size_t page, int ct){
    //generate a list of fixed-parameter positions, and their paramvals.
   apop_data * mask = apop_map((apop_data*)data, find_nans, .all_pages='y');
    //find out where the NaNs are
    for (size_t i=0; mask->vector && i< mask->vector->size; i++)
            if (apop_data_get(mask, i, -1))
                addin(predict, i, -1, page);
    for (size_t i=0; mask->matrix && i< mask->matrix->size1; i++)
        for (int j=0; j <mask->matrix->size2; j++)
            if (apop_data_get(mask, i, j))
                addin(predict, i, j, page);
//        apop_assert(mset->ct, 0, 0,'s',"You're asking me to estimate a model where every single parameter is fixed.");
    if (mask->more)
        ct += apop_sum(mask->vector)+ apop_matrix_sum(mask->matrix)
                     + find_missing(mask->more, predict, page+1, ct);
    apop_data_free(mask);
    return ct;
}

#include "output.h"
apop_data *apop_predict_table_prep(apop_data *in, char fill_with_nans){
    apop_data *out = apop_data_alloc(0, 0, 0);
    if (in)
        apop_data_add_page(in, out, "<predict>");
    else 
        sprintf(out->names->title, "<predict>");
    apop_name_add(out->names, "row", 'c');
    apop_name_add(out->names, "col", 'c');
    apop_name_add(out->names, "page", 'c');
    apop_name_add(out->names, "observed", 'c');
    apop_name_add(out->names, "predicted", 'c');
    if (fill_with_nans == 'y')
        find_missing(in, out, 0, 0);
    return out;
}

/** Take a \c predict table and set the entries in the data set to the given predicted
  value. Functions for prediction and imputation use this internally, and append to your
  data a \c predict table of the right form.  For example, \c apop_ml_impute uses
  this internally.
  
  I assume that the ordering of elements in the \c predict table include everything on the
  first page, then everything on the second, et cetera. 

\param data The data set to be filled in. It should have a page named \c
\f$<\f$predict\f$>\f$.
\param predict If your data set doesn't have a \c \f$<\f$predict\f$>\f$ page, then just
provide one this way; else let this be \c NULL;
*/
void apop_data_predict_fill(apop_data *data, apop_data *predict){
    if (!predict)
        predict = apop_data_get_page (data, "<predict>");
    if (!predict) return;
    int this_page_ct = 0;
    apop_data *this_page = data;
    for (int i=0; i < predict->matrix->size1; i++){
        int p = apop_data_get(predict, .row=i, .colname="page");
        if (p != this_page_ct){
            this_page_ct = p; 
            this_page = this_page->more;
        }
        apop_data_set(this_page, .row= apop_data_get(predict, .row=i, .colname="row"),
                                 .col= apop_data_get(predict, .row=i, .colname="col"),
                                 .val= apop_data_get(predict, .row=i, .colname="predicted"));
    }
}

/** If there is an NaN anywhere in the row of numeric data (including the matrix, the vector, and the weights) then delete the row from the data set.

The function returns a new data set with the NaNs removed, so the original data set is left unmolested. You may want to \c apop_data_free the original immediately after this function.

\li If every row has an NaN, then this returns \c NULL.
\li If there is text, it gets pruned as well.
\li If \c inplace = 'y', then I'll free each element of the input data
    set and refill it with the pruned elements. Again, I'll take up (up to)
    twice the size of the data set in memory during the function. If
    every row has an NaN, then your \c apop_data set will have a lot of
    \c NULL elements.
\li I only look at the first page of data (i.e. the \c more element is ignored).
\li This function uses the \ref designated syntax for inputs.

    \param d    The data, with NaNs
    \param inplace If \c 'y', clear out the pointer-to-\ref apop_data that
    you sent in and refill with the pruned data. If \c 'n', leave the
    set alone and return a new data set.
    \return     A (potentially shorter) copy of the data set, without
    NaNs. If <tt>inplace=='y'</tt>, redundant with the input. If the entire data set is
    cleared out, then this will be \c NULL.
*/
APOP_VAR_HEAD apop_data * apop_data_listwise_delete(apop_data *d, char inplace){
    apop_data * apop_varad_var(d, NULL);
    if (!d) return NULL;
    char apop_varad_var(inplace, 'n');
APOP_VAR_ENDHEAD
    Get_vmsizes(d) //defines firstcol, vsize, wsize, msize1, msize2.
    apop_assert(msize1 || vsize, NULL, 0, 'c', 
            "You sent to apop_data_listwise_delete a data set with void matrix and vector. Confused, it is returning NULL.\n");
    //find out where the NaNs are
    int len = vsize ? vsize : msize1;
    int marked[len], not_empty = 0;
    memset(marked, 0, sizeof(int)*len);
    for (int i=0; i< (vsize ? vsize: msize1); i++)
        for (int j=firstcol; j <msize2; j++){
            if (gsl_isnan(apop_data_get(d, i, j))){
                    marked[i] = 1;
                    break;
            }
        }
    for (int i=0; i< wsize; i++)
        if (gsl_isnan(gsl_vector_get(d->weights, i)))
            marked[i] = 1;
    //check that at least something isn't NULL.
    for (int i=0; i< len; i++)
        if (!marked[i]){
            not_empty ++;
            break;
        }
    if (!not_empty)
        return NULL;
    apop_data *out = (inplace=='y'|| inplace=='Y') ? d : apop_data_copy(d);
    apop_data_rm_rows(out, marked);
    return out;
}

//ML imputation

/** \hideinitializer */
#define Switch_back    \
    apop_data *real_data = ml_model->parameters;   \
    apop_model *actual_base = ml_model->more; \
    actual_base->parameters = d; 

static apop_model * i_est(apop_data *d, apop_model *ml_model){
    Switch_back
    return apop_estimate(real_data, *actual_base);
}

#include "mapply.h"
static double i_ll(apop_data *d, apop_model *ml_model){
    Switch_back
    return apop_log_likelihood(real_data, actual_base);
}

static double i_p(apop_data *d, apop_model *ml_model){
    Switch_back
    return apop_p(real_data, actual_base);
}

static apop_model apop_ml_impute_model = {"Internal ML imputation model", .estimate=i_est, .p = i_p, .log_likelihood=i_ll};

/** Impute the most likely data points to replace NaNs in the data, and insert them into 
the given data. That is, the data set is modified in place.

How it works: this uses the machinery for \c apop_model_fix_params. The only difference is 
that this searches over the data space and takes the parameter space as fixed, while basic 
fix params model searches parameters and takes data as fixed. So this function just does the
necessary data-parameter switching to make that happen.

\param  d       The data set. It comes in with NaNs and leaves entirely filled in.
\param  mvn A parametrized \c apop_model from which you expect the data was derived.
if \c NULL, then I'll use the Multivariate Normal that best fits the data after listwise deletion.

\return An estimated <tt>apop_ml_impute_model</tt>. Also, the data input will be filled in and ready to use.
*/
apop_model * apop_ml_impute(apop_data *d,  apop_model* mvn){
    if (!mvn){
        apop_data *list_d = apop_data_listwise_delete(d);
        apop_assert_s(list_d, "Listwise deletion returned no whole rows, "
                            "so I couldn't fit a Multivariate Normal to your data. "
                            "Please provide a pre-estimated initial model.");
        mvn = apop_estimate(list_d, apop_multivariate_normal);
        apop_data_free(list_d);
    }
    apop_model *impute_me = apop_model_copy(apop_ml_impute_model);
    impute_me->parameters = d;
    impute_me->more = mvn;
    apop_model *fixed = apop_model_fix_params(impute_me);
//    Apop_model_add_group(fixed, apop_mle, .want_cov='n', .dim_cycle_tolerance=1);
    Apop_settings_set(fixed, apop_mle, want_cov, 'n');
    apop_model *m = apop_estimate(mvn->parameters, *fixed);
    apop_data_memcpy(d, m->parameters); //A bit inefficient.
    return m;
}


/**
Imputation (in this context) is the process of finding fill-in values for missing data
points. Filling in values from a single imputation and then returning the values will give
you a single complete data set, but statistics you derive from that data set don't reflect
the uncertainty inherent in using artifical, model-derived data rather than actual
observations. This function uses several imputations and a statistic-calculating function
you provide to find the total variance (see details below).


The multiple imputation process involves two steps:

\li Generating imputations via a model of your choosing, and writing down the results in a
\c fill_ins table, described in the parameter list. You've already done this by the time
you call this function.

\li For a given complete data set, generating a summary statistic, such as the mean of a
column, or the ratio of columns one and two; plus the variance of your statistic(s).

This function takes the output from the first step (a list of fill-ins), and uses it to
calculate a list of statistics/variances, as per the second step. After generating this
list of statistics/variances, it ties them together to produce a single best estimate of
the statistic and its full variance.


\param stat A function that takes in a single \ref apop_data set, and calculates
statistics and their covariances. The output should have two pages. The first
will be the statistics themselves; the second will be the covariance matrix. If
I find a page with the name <tt>\<Covariance\></tt> then I will use that; else
the second page. This rule means you can return the \c parameters from most estimated
models.

\param base_data The data, with \c NaNs to be filled in. When calculating the statistics,
I fill in the values in the \c base_data set directly, so it is modified, and in
the end will have the value of the last replication.

\param fill_ins This is the list of values to fill in to the base data set, and it must be
an \ref apop_data set whose matrix includes the following two column names toward the
beginning: \c row and \c col, with an optional \c page. Every column after those two or
three columns is taken to be an imputation that can be used to fill in values. That is, I
will first take the first column after the row/col/page column and plug its values
into the corresponding row/col/pages of \c base_data, calculate the variances, then repeat
with the second column, and so on.

\return The first page is the mean of each replicate's statistics; the second page is the
overall covariance (and will have the page title <tt>\<Covariance\></tt>).
Let \f$S_i\f$ be the covariance for each replicate \f$i\f$; let there be \f$m\f$
replicates; let \f$\mu(\cdot)\f$ indicate the mean; then the overall covariance is the mean of the individual variances plus 

\f$\mu(S_i) + {\rm var}(S_i)/(1+1/m).\f$


\li Multiple pages for input data are not yet implemented.
*/

apop_data * apop_multiple_imputation_variance(apop_data *(*stat)(apop_data *), apop_data *base_data, 
                                                                               apop_data *fill_ins){
    /*Copyright: this function is part of a larger work (C) Ben Klemens, but was partially written
    by a government employee during work hours. Some lawyers will tell you that the code
    is licensed via GPL v2, like the main work; others will tell you that it's public domain.*/

    //Part I: call the statistic-calculating function with each filled-in replicate
	int row_column = apop_name_find(fill_ins->names, "row", 'c');
	int col_column = apop_name_find(fill_ins->names, "col", 'c');
	int page_column = apop_name_find(fill_ins->names, "page", 'c');
	int first_non_address= GSL_MAX(row_column, GSL_MAX(col_column, page_column)) + 1;
	int replicates = fill_ins->matrix->size2 - first_non_address;
    apop_assert_s(replicates, "I couldn't find anything besides row/col/page columns "
            "(which I may or may not have found)")

	apop_data *estimates[replicates];
	for (int i= first_non_address; i< fill_ins->matrix->size2; i++){
		for (int j=0; j < fill_ins->matrix->size1; j++)
			apop_data_set(base_data, 
				.row = apop_data_get(fill_ins, j, row_column),
				.col = apop_data_get(fill_ins, j, col_column),
				.val =  apop_data_get(fill_ins, j, i));
		estimates[i-first_non_address] = stat(base_data);
	}

    //Part II: find the mean of the statistics and the total variance of the cov matrix.
	gsl_vector *vals = gsl_vector_alloc(replicates);
    apop_data *out = apop_data_copy(estimates[0]);
	//take the simple mean of the main data set.
	{ //this limits the scope of the Get_vmsizes macro.
	 Get_vmsizes(estimates[0]); 
     for (int j=0; j < msize2; j++)
         for (int i=0; i < (vsize ? vsize : msize1); i++){
            for (int k=0; k< replicates; k++)
                gsl_vector_set(vals, k, apop_data_get(estimates[k], i, j));
             apop_data_set(out, i, j, apop_vector_mean(vals));
         }
	}
    apop_data *out_var = apop_data_get_page(estimates[0], "<Covariance>");
    int cov_is_labelled = out_var !=NULL;
    if (!cov_is_labelled){
        sprintf(out->more->names->title, "<Covariance>");
        out_var = estimates[0]->more;
    }
	Get_vmsizes(out_var);
    for (int i=0; i < msize1; i++)
        for (int j=0; j < msize2; j++){
            for (int k=0; k< replicates; k++){
                apop_data *this_p = cov_is_labelled ? apop_data_get_page(estimates[k], "<Covariance>")
                                        : estimates[k]->more;
                gsl_vector_set(vals, k, apop_data_get(this_p, i, j));
            }
            apop_data_set(out_var, i, j, apop_vector_mean(vals) + apop_var(vals)/(1+1./replicates));
        }
    return out;	
}
