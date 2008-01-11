/** \file apop_OLS.c

  OLS models. Much of the real work is done in regression.c.

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "regression.h"
#include "stats.h"
#include "asst.h"

//shift first col to depvar, rename first col "one".
static void prep_names (apop_model *e){
  int i;
  apop_ls_settings   *p = e->model_settings;
	if (e->data->names->colct > 0) {		
		//apop_name_add(n, n->column[0], 'd');
        apop_name_add(e->expected->names, e->data->names->column[0], 'c');
        apop_name_add(e->expected->names, "predicted", 'c');
        apop_name_add(e->expected->names, "residual", 'c');
        if (e->parameters)
            snprintf(e->parameters->names->title, 100, "Regression of %s", e->data->names->column[0]);
		sprintf(e->data->names->column[0], "1");
        apop_name_add(e->parameters->names, "1", 'r');
        apop_name_add(e->parameters->names, "parameters", 'v');
        for(i=1; i< e->data->names->colct; i++)
            apop_name_add(e->parameters->names, e->data->names->column[i], 'r');
        if (p->want_cov){
            if (e->data->names){
                apop_name_stack(e->covariance->names, e->data->names, 'c');
                apop_name_cross_stack(e->covariance->names, e->data->names, 'r', 'c');
            }
		    sprintf(e->covariance->names->column[0], "1");
		    sprintf(e->covariance->names->row[0], "1");
        }
	}
}

static void ols_prep(apop_data *d, apop_model *m){
    if (!d->vector){
        APOP_COL(d, 0, independent);
        d->vector = apop_vector_copy(independent);
        gsl_vector_set_all(independent, 1);
        if (d->names->colct > 0) {		
            apop_name_add(d->names, d->names->column[0], 'v');
            sprintf(d->names->column[0], "1");
        }
    }
    void *mpt = m->prep; //also use the defaults.
    m->prep = NULL;
    apop_model_prep(d, m);
    m->prep = mpt;
}

/** The assumption that makes a log likelihood possible is that the
errors are normally distributed.

This function is a bit inefficient, in that it calculates the error terms,
which you may have already done in the OLS estimation.

 */
static double ols_log_likelihood (apop_data *d, apop_model *p){ 
  apop_assert(p->parameters, 0, 0,'s', "You asked me to evaluate an un-parametrized model. Returning zero.");
  int         i; 
  long double	ll  = 0; 
  double      sigma, expected, actual;
  gsl_matrix	*data		    = d->matrix;
  gsl_vector  *errors         = gsl_vector_alloc(data->size1);
	for(i=0;i< data->size1; i++){
        APOP_ROW(d, i, datarow);
        gsl_blas_ddot(p->parameters->vector, datarow, &expected);
        if (d->vector){ //then this has been prepped
            actual       = apop_data_get(d,i, -1);
        } else {
            actual       = gsl_matrix_get(data,i, 0);
            expected    += gsl_vector_get(p->parameters->vector,0) * (1 - actual); //data isn't affine.
        }
        gsl_vector_set(errors, i, expected-actual);
    }
    sigma   = sqrt(apop_vector_var(errors));
	for(i=0;i< data->size1; i++){
        ll  += log(gsl_ran_gaussian_pdf(gsl_vector_get(errors, i), sigma));
	} 
    gsl_vector_free(errors);
    return ll;
}

/* $\partial {\cal N}(x\beta - y)/\partial \beta_i = \sum{x_i} \partial {\cal N}(K)/\partial K$ (at $K=x\beta -y$) */
static void ols_score(apop_data *d, gsl_vector *gradient, apop_model *p){ 
  apop_assert_void(p->parameters, 0,'s', "You asked me to evaluate an un-parametrized model. Not changing the gradient");
  size_t         i, j; 
  double      sigma, expected, actual;
  gsl_matrix	*data		    = d->matrix;
  gsl_vector  *errors         = gsl_vector_alloc(data->size1);
  gsl_vector  *normscore      = gsl_vector_alloc(2);
  apop_data  *subdata      = apop_data_alloc(0,1,1);
	for(i=0;i< data->size1; i++){
        APOP_ROW(d, i, datarow);
        gsl_blas_ddot(p->parameters->vector, datarow, &expected);
        if (d->vector){ //then this has been prepped
            actual       = apop_data_get(d,i, -1);
        } else {
            actual       = gsl_matrix_get(data,i, 0);
            expected    += gsl_vector_get(p->parameters->vector,0) * (1 - actual); //data isn't affine.
        }
        gsl_vector_set(errors, i, expected-actual);
    }
    sigma   = sqrt(apop_vector_var(errors));
    apop_model *norm = apop_model_set_parameters(apop_normal, 0.0, sigma);
    gsl_vector_set_all(gradient, 0);
	for(i=0;i< data->size1; i++){
        apop_data_set(subdata, 0, 0, gsl_vector_get(errors, i));
        apop_score(subdata, normscore, norm);
        for(j=0; j< data->size2; j++)
            apop_vector_increment(gradient, j, apop_data_get(d, i, j) * gsl_vector_get(normscore, 0));
	} 
    gsl_vector_free(errors);
    apop_model_free(norm);
}





static void xpxinvxpy(gsl_matrix *data, gsl_vector *y_data, gsl_matrix *xpx, gsl_vector* xpy, apop_model *out){
  apop_ls_settings   *p = out->model_settings;
	if (p->want_cov + p->want_expected_value == 0 ){	
		//then don't calculate (X'X)^{-1}
		gsl_linalg_HH_solve (xpx, xpy, out->parameters->vector);
		return;
	} //else:
  gsl_vector 	*error = gsl_vector_alloc(data->size1);
  gsl_vector 	predicted;
  gsl_matrix	*cov;
  double        s_sq;
	cov	= gsl_matrix_alloc(data->size2, data->size2);
	apop_det_and_inv(xpx, &cov, 0, 1);	    //not yet cov, just (X'X)^-1.
	gsl_blas_dgemv(CblasNoTrans, 1, cov, xpy, 0, out->parameters->vector);      // \beta=(X'X)^{-1}X'Y
	gsl_blas_dgemv(CblasNoTrans, 1, data, out->parameters->vector, 0, error);   // X'\beta ==predicted
	gsl_vector_sub(error,y_data);           //X'\beta - Y == error
    gsl_blas_ddot(error, error, &s_sq);   // e'e
    s_sq    /= data->size1 - data->size2;   //\sigma^2 = e'e / df
	gsl_matrix_scale(cov, s_sq);            //cov = \sigma^2 (X'X)^{-1}
	if (p->want_expected_value){
        gsl_matrix_set_col(out->expected->matrix, 0, y_data);
        gsl_matrix_set_col(out->expected->matrix, 2, error);
        predicted   = gsl_matrix_column(out->expected->matrix, 1).vector;
        gsl_vector_set_zero(&predicted);
        gsl_vector_add(&predicted, y_data);
        gsl_vector_sub(&predicted, error);
    }
    gsl_vector_free(error);
    out->covariance->matrix	= cov;
}

static apop_model * apop_estimate_OLS(apop_data *inset, apop_model *ep){
    apop_assert(inset,  NULL, 0,'s', "You asked me to estimate a regression with NULL data.");
  apop_data         *set;
  apop_model       *epout;
  gsl_vector        *weights    = NULL;
  int               i;
  apop_ls_settings  *olp;
    if (!ep) {
        olp             = apop_ls_settings_alloc(inset, apop_OLS);
        epout           = olp->model;
    } else if (!ep->model_settings) {
        olp             = apop_ls_settings_alloc(inset, *ep);
        epout           = olp->model;
    } else {
        olp             = ep->model_settings;
        epout           = ep;
    }
    epout->data = inset;
    if(epout->parameters)
        apop_data_free(epout->parameters);
    epout->parameters = apop_data_alloc(inset->matrix->size2,0,0);
    set = olp->destroy_data ? inset : apop_data_copy(inset); 
    
    //prep weights.
    if (olp->destroy_data)
        weights = epout->data->weights;  //may be NULL.
    else
        weights = apop_vector_copy(epout->data->weights); //may be NULL.
    if (weights)
        for (i =0; i< weights->size; i++)
            gsl_vector_set(weights, i, sqrt(gsl_vector_get(weights, i)));

  gsl_vector    *y_data     = gsl_vector_alloc(set->matrix->size1); 
  gsl_vector    *xpy        = gsl_vector_calloc(set->matrix->size2);
  gsl_matrix    *xpx        = gsl_matrix_calloc(set->matrix->size2, set->matrix->size2);
    if (olp->want_expected_value)
        epout->expected   = apop_data_alloc(0, set->matrix->size1, 3);
    if (olp->want_cov)
        epout->covariance = apop_data_alloc(0, set->matrix->size2, set->matrix->size2);
    prep_names(epout);
    APOP_COL(set, 0, firstcol);
    gsl_vector_memcpy(y_data,firstcol);
    gsl_vector_set_all(firstcol, 1);     //affine: first column is ones.
    if (weights){
        gsl_vector_mul(y_data, weights);
        for (i = 0; i < set->matrix->size2; i++){
            APOP_COL(set, i, v);
            gsl_vector_mul(v, weights);
        }
    }

    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set->matrix, set->matrix, 0, xpx);   //(X'X)
    gsl_blas_dgemv(CblasTrans, 1, set->matrix, y_data, 0, xpy);       //(X'y)
    xpxinvxpy(set->matrix, y_data, xpx, xpy, epout);
    gsl_vector_free(y_data); gsl_vector_free(xpy);

    epout->llikelihood  = ols_log_likelihood(epout->data, epout);
    if (!olp->destroy_data)
        apop_data_free(set);
    apop_estimate_parameter_t_tests(epout);
    epout->status       = 1;
    return epout;
}

/** Ordinary least squares.

\ingroup models

\param inset The first column is the dependent variable, and the remaining columns the independent. Is destroyed in the process, so make a copy beforehand if you need.

\param epin    An \ref apop_model object. The only
thing we look at is the \c destroy_data element. If this is NULL or
\c destroy_data==0, then the entire data set is copied off, and then
mangled. If \c destroy_data==1, then this doesn't copy off the data set,
but destroys it in place.

\return
Will return an \ref apop_model <tt>*</tt>.
<tt>The_result->parameters</tt> will hold the coefficients; the first
coefficient will be the coefficient on the constant term, and the
remaining will correspond to the independent variables. It will therefore
be of size <tt>(data->size2)</tt>. Do not pre-allocate.

If you asked for it, the covariance matrix, confidence levels, and residuals will also be returned. The confidence intervals give the level of certainty with which we can reject the hypothesis that the given coefficient is zero.

\b sample 

First, you will need a file named <tt>data</tt> in comma-separated form. The first column is the dependent variable; the remaining columns are the independent. For example:
\verbatim
Y, X_1, X_2, X_3
2,3,4,5
1,2,9,3
4,7,9,0
2,4,8,16
1,4,2,9
9,8,7,6
\endverbatim

The program:
\code
#include <apop.h>

int main(void){
apop_data       *data;
apop_model   *est;
    apop_text_to_db("data","d",0,1,NULL);
    data = apop_query_to_data("select * from d");
    est  = apop_estimate(data, apop_OLS);
    apop_model_show(est);
    return 0;
}
\endcode

If you saved this code to <tt>sample.c</tt>, then you can compile it with
\verbatim
gcc sample.c -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endverbatim

and then run it with <tt>./run_me</tt>. Alternatively, you may prefer to compile the program using a \ref makefile .

Feeling lazy? The program above was good form and demonstrated useful
features, but the code below will do the same thing in two lines:

\code
#include <apop.h>
int main(){ apop_model_show(apop_estimate(apop_text_to_data("data", 0, 0), apop_OLS)); }
\endcode


 */
apop_model apop_OLS = {.name="OLS", .vbase = -1, .estimate =apop_estimate_OLS, 
                            .log_likelihood = ols_log_likelihood, .score=ols_score, .prep = ols_prep};

/** The GLS model

  This is basically a wrapper for the GLS regression function, \ref apop_params_GLS.
\ingroup models
*/
//apop_model apop_GLS = {"GLS", -1, apop_params_GLS, NULL, NULL, NULL, NULL, NULL};





//Instrumental variables






static apop_data *prep_z(apop_data *x, apop_data *instruments){
  int       i;
  apop_data *out    = apop_data_copy(x);
    if (instruments->vector)
        for (i=0; i< instruments->vector->size; i++){
            APOP_COL(instruments, i, inv);
            APOP_COL(out, instruments->vector->data[i], outv);
            gsl_vector_memcpy(outv, inv);
        }
    else if (instruments->names->rowct)
        for (i=0; i< instruments->names->rowct; i++){
            int rownumber = apop_name_find(x->names, instruments->names->row[i], 'c');
            apop_assert(rownumber != -1,  NULL, 0, 's', "You asked me to substitute instrument column %i for the data column named %s, but I could find no such name.",  i, instruments->names->row[i]);
            APOP_COL(instruments, i, inv);
            APOP_COL(out, rownumber, outv);
            gsl_vector_memcpy(outv, inv);
        }
    else 
        apop_error(0, 's', "%s: Your instrument matrix has data, but neither a vector element nor row names indicating what columns in the original data should be replaced.\n", __func__);
    return out;
}


static apop_model * apop_estimate_IV(apop_data *inset, apop_model *ep){
  apop_assert(inset, NULL, 0,'s', "You asked me to estimate a regression with NULL data.");
  apop_data         *set, *z;
  apop_model       *epout;
  gsl_vector        *weights    = NULL;
  int               i;
  apop_ls_settings  *olp;
    if (!ep || !ep->model_settings)
        return apop_estimate(inset, apop_OLS);
    olp                 = ep->model_settings;
    olp->want_cov       = 0;//not working yet.
    epout               = apop_model_copy(*ep);
    epout->model_settings = malloc(sizeof(apop_ls_settings));
    memcpy(epout->model_settings, olp, sizeof(apop_ls_settings));
    if (!olp->instruments || !olp->instruments->matrix->size2) 
        return apop_estimate(inset, apop_OLS);
    epout->data = inset;
    if(epout->parameters)
        apop_data_free(epout->parameters);
    epout->parameters = apop_data_alloc(inset->matrix->size2,0,0);
    set = olp->destroy_data ? inset : apop_data_copy(inset); 
    z   = prep_z(inset, olp->instruments);
    
    //prep weights.
    if (olp->destroy_data)
        weights = epout->data->weights;  //may be NULL.
    else
        weights = apop_vector_copy(epout->data->weights); //may be NULL.
    if (weights)
        for (i =0; i< weights->size; i++)
            gsl_vector_set(weights, i, sqrt(gsl_vector_get(weights, i)));

  apop_data    *y_data     = apop_data_alloc(set->matrix->size1, 0, 0); 
    if (olp->want_expected_value)
        epout->expected   = apop_data_alloc(0, set->matrix->size1, 3);
    if (olp->want_cov)
        epout->covariance = apop_data_alloc(0, set->matrix->size1, set->matrix->size1);
    prep_names(epout);
    APOP_COL(set, 0, firstcol);
    gsl_vector_memcpy(y_data->vector,firstcol);
    gsl_vector_set_all(firstcol, 1);     //affine: first column is ones.
    if (weights){
        gsl_vector_mul(y_data->vector, weights);
        for (i = 0; i < set->matrix->size2; i++){
            APOP_COL(set, i, v);
            gsl_vector_mul(v, weights);
        }
    }

    apop_data *zpx    = apop_dot(z, set, 1, 0);
    apop_data *zpy    = apop_dot(z, y_data, 1, 0);
    apop_data *zpxinv = apop_matrix_to_data(apop_matrix_inverse(zpx->matrix));
    epout->parameters = apop_dot(zpxinv, zpy, 0);
    apop_data_free(y_data);
    apop_data_free(zpx); 
    apop_data_free(zpxinv);
    apop_data_free(zpy);

    if (!olp->destroy_data)
        apop_data_free(set);
//    apop_estimate_parameter_t_tests(epout);
    epout->status   = 1;
    return epout;
}

/** Instrumental variable regression
\ingroup models

 Operates much like the \ref apop_estmate_OLS function, but the input
 parameters also need to have a table of substitutions. The vector
 element of the table lists the column numbers to be substituted (the
 dependent var is zero; first independent col is one), and then one
 column for each item to substitute.

If the vector of your apop_data set is NULL, then I will use the row
names to find the columns to substitute. This is generally more robust
and/or convenient.

If the \c instruments data set is somehow NULL or empty, I'll just run OLS.

\code
apop_ls_settings *ivp = apop_ls_settings_alloc(data, apop_IV);
ivp->instruments    = apop_data_alloc(data->matrix->size1, 2);
APOP_COL(ivp->instruments, 0, firstcol);
gsl_vector_memcpy(firstcol, your_data_vector);
APOP_COL(ivp->instruments, 1, secondcol);
gsl_vector_memcpy(firstcol, your_other_data_vector);
apop_name_add(ivp->names, "subme_1", 'r');
apop_name_add(ivp->names, "subme_2", 'r');
apop_estimate(data, ivp->model);
\endcode


\todo This function does some serious internal data copying. It would be
only slightly more human- and labor-intensive to do the linear algebra
without producing the Z matrix explicitly.

 */
apop_model apop_IV = {.name="instrumental variables", .vbase = -1, .estimate =apop_estimate_IV, 
                            .log_likelihood = ols_log_likelihood};
