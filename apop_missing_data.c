/** \file apop_missing_data.c Some missing data handlers. */
/* Copyright (c) 2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include <regex.h>


/** If there is an NaN anywhere in the row of data (including the matrix, the vector, the weights, and the text) then delete the row from the data set.

\li If every row has an NaN, then this returns \c NULL.
\li If \c apop_opts.db_nan is not \c NULL, then I will use that as a regular expression to check the text elements for bad data as well.
\li If \c inplace = 'y', then I'll free each element of the input data
    set and refill it with the pruned elements. I'll still take up (up to)
    twice the size of the data set in memory during the function. If
    every row has an NaN, then your \c apop_data set will end up with
    \c NULL vector, matrix, .... if \c inplace = 'n', then the original data set is left unmolested.
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
    apop_assert_c(msize1 || vsize || d->textsize[0], NULL, 0, 
            "You sent to apop_data_listwise_delete a data set with NULL matrix, NULL vector, and no text. "
            "Confused, it is returning NULL.");
    //find out where the NaNs are
    int len = GSL_MAX(vsize ? vsize : msize1, d->textsize[0]); //still some size assumptions here.
    int not_empty = 0;
    int *marked = calloc(len, sizeof(int));
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
    if (d->textsize[0] && apop_opts.db_nan){
        regex_t    rex;
        int compiled_ok = !regcomp(&rex, apop_opts.db_nan, REG_EXTENDED +  REG_ICASE + REG_NOSUB);
        apop_assert(compiled_ok, "apop_opts.db_nan needs to be a regular expression that "
                                "I can use to check the text element of your data set for "
                                "NaNs, But compiling %s into a regex failed. Or, set "
                                "apop_opts.db_nan=NULL to bypass text checking.", apop_opts.db_nan);
        for(int i=0; i< d->textsize[0]; i++)
            if (!marked[i])
                for(int j=0; j< d->textsize[1]; j++)
                    if (!regexec(&rex, d->text[i][j], 0, 0, 0)){
                        marked[i] ++;
                        break;
                    }
        regfree(&rex);
    }

    //check that at least something isn't NULL.
    for (int i=0; i< len; i++)
        if (!marked[i]){
            not_empty ++;
            break;
        }
    if (!not_empty){
        free(marked);
        return NULL;
    }
    apop_data *out = (inplace=='y'|| inplace=='Y') ? d : apop_data_copy(d);
    apop_data_rm_rows(out, marked);
    free(marked);
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

How it works: this uses the machinery for \ref apop_model_fix_params. The only difference is 
that this searches over the data space and takes the parameter space as fixed, while basic 
fix params model searches parameters and takes data as fixed. So this function just does the
necessary data-parameter switching to make that happen.

\param  d       The data set. It comes in with NaNs and leaves entirely filled in.
\param  mvn A parametrized \ref apop_model from which you expect the data was derived.
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
/*    if(apop_settings_get_group(mvn, apop_mle))
        apop_settings_copy_group(fixed, mvn, "apop_mle");
    else Apop_model_add_group(fixed, apop_mle);
    */Apop_settings_set(fixed, apop_mle, want_cov, 'n');
    Apop_model_add_group(fixed, apop_parts_wanted);
    apop_model *m = apop_estimate(mvn->parameters, *fixed);
    apop_model *filled = apop_model_fix_params_get_base(m);
    apop_data_memcpy(d, filled->parameters); //A bit inefficient.
    return m;
}
