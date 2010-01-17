/** \file apop_OLS.c

  OLS models. Much of the real work is done in apop_regression.c.*/
/* Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "conversions.h"
#include "settings.h"
#include "stats.h"
#include "asst.h"

void * apop_ls_settings_copy(apop_ls_settings *in){
  apop_ls_settings *out  = malloc(sizeof(*out));
    *out =  *in;
    out->instruments = apop_data_copy(in->instruments);
    return out;
}

void apop_ls_settings_free(apop_ls_settings *in){ free(in); }

/** Initialize the settings for a least-squares--type model. For use
  with \ref Apop_model_add_group.  See \ref apop_ls_settings for the possible elements to set.
  */
apop_ls_settings * apop_ls_settings_init(apop_ls_settings in){
  apop_ls_settings *out  = calloc(1, sizeof(*out));
    apop_varad_setting(in, out, want_cov, 'y');
    apop_varad_setting(in, out, want_expected_value, 'y');
    if (out->want_cov == 1) out->want_cov = 'y';
    if (out->want_expected_value == 1) out->want_expected_value = 'y';
    return out;
}

/*apop_ls_settings * apop_ls_settings_alloc(apop_data *data){
    return apop_ls_settings_init((apop_ls_settings){ }); }*/

//shift first col to depvar, rename first col "one".
static void prep_names (apop_model *e){
  apop_ls_settings   *p = apop_settings_get_group(e, "apop_ls");
    if (e->expected){
        apop_name_add(e->expected->names, (e->data->names->colct ? e->data->names->column[0] : "expected"), 'c');
        apop_name_add(e->expected->names, "predicted", 'c');
        apop_name_add(e->expected->names, "residual", 'c');
    }
	if (e->data->names->colct > 0) {		
        if (e->parameters)
            snprintf(e->parameters->names->title, 100, "Regression of %s", e->data->names->column[0]);
		sprintf(e->data->names->column[0], "1");
        apop_name_add(e->parameters->names, "1", 'r');
        apop_name_add(e->parameters->names, "parameters", 'v');
        for(int i=1; i< e->data->names->colct; i++)
            apop_name_add(e->parameters->names, e->data->names->column[i], 'r');
        if (p->want_cov== 'y'){
            if (e->data->names){
                apop_name_stack(e->covariance->names, e->data->names, 'c');
                apop_name_stack(e->covariance->names, e->data->names, 'r', 'c');
            }
		    sprintf(e->covariance->names->column[0], "1");
		    sprintf(e->covariance->names->row[0], "1");
        }
	}
}

static void ols_shuffle(apop_data *d){
    if (!d->vector){
        APOP_COL(d, 0, independent);
        d->vector = apop_vector_copy(independent);
        gsl_vector_set_all(independent, 1);     //affine; first column is ones.
        if (d->names->colct > 0) {		
            apop_name_add(d->names, d->names->column[0], 'v');
            sprintf(d->names->column[0], "1");
        }
    }
}

static void ols_prep(apop_data *d, apop_model *m){
    ols_shuffle(d);
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
  long double	ll  = 0; 
  long double      sigma, actual, weight;
  double expected;
  int use_weights =  (!strcmp(p->name, "Weighted Least Squares") && d->weights);
  gsl_matrix	*data		    = d->matrix;
  gsl_vector  *errors         = gsl_vector_alloc(data->size1);
	for (size_t i=0;i< data->size1; i++){
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
	for(size_t i=0;i< data->size1; i++){
        weight       = use_weights ? gsl_vector_get(d->weights, i) : 1; 
        ll  += log(gsl_ran_gaussian_pdf(gsl_vector_get(errors, i), sigma)* weight);
	} 
    gsl_vector_free(errors);
    return ll;
}

/* $\partial {\cal N}(x\beta - y)/\partial \beta_i = \sum{x_i} \partial {\cal N}(K)/\partial K$ (at $K=x\beta -y$) */
static void ols_score(apop_data *d, gsl_vector *gradient, apop_model *p){ 
  apop_assert_void(p->parameters, 0,'s', "You asked me to evaluate an un-parametrized model. Not changing the gradient");
  long double      sigma, actual, weight;
  double expected;
  gsl_matrix	*data		    = d->matrix;
  gsl_vector  *errors         = gsl_vector_alloc(data->size1);
  gsl_vector  *normscore      = gsl_vector_alloc(2);
  apop_data  *subdata      = apop_data_alloc(0,1,1);
  int use_weights =  (!strcmp(p->name, "Weighted Least Squares") && d->weights);
	for(size_t i=0;i< data->size1; i++){
        APOP_ROW(d, i, datarow);
        gsl_blas_ddot(p->parameters->vector, datarow, &expected);
        if (d->vector){ //then this has been prepped
            actual       = apop_data_get(d,i, -1);
        } else {
            actual       = gsl_matrix_get(data,i, 0);
            expected    +=  gsl_vector_get(p->parameters->vector,0) * (1 - actual); //data isn't affine.
        }
        gsl_vector_set(errors, i, expected-actual);
    }
    sigma   = sqrt(apop_vector_var(errors));
    apop_model *norm = apop_model_set_parameters(apop_normal, 0.0, sigma);
    gsl_vector_set_all(gradient, 0);
	for(size_t i=0;i< data->size1; i++){
        apop_data_set(subdata, 0, 0, gsl_vector_get(errors, i));
        apop_score(subdata, normscore, norm);
        weight       = use_weights ? gsl_vector_get(d->weights, i) : 1; 
        for(size_t j=0; j< data->size2; j++)
            apop_vector_increment(gradient, j, weight * apop_data_get(d, i, j) * gsl_vector_get(normscore, 0));
	} 
    gsl_vector_free(errors);
    apop_model_free(norm);
}


static void xpxinvxpy(gsl_matrix *data, gsl_vector *y_data, gsl_matrix *xpx, gsl_vector* xpy, apop_model *out){
  apop_ls_settings   *p =  apop_settings_get_group(out, "apop_ls");
	if ((p->want_cov!='y') && (p->want_expected_value != 'y') ){	
		//then don't calculate (X'X)^{-1}
		gsl_linalg_HH_solve (xpx, xpy, out->parameters->vector);
		return;
	} //else:
  gsl_vector 	*error = gsl_vector_alloc(data->size1);
  gsl_vector 	predicted;
  gsl_matrix	*cov;
  double        s_sq;
	cov	= apop_matrix_inverse(xpx);	    //not yet cov, just (X'X)^-1.
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
        gsl_vector_memcpy(&predicted, y_data);
        gsl_vector_add(&predicted, error); //pred = y_data + error
    }
    gsl_vector_free(error);
    if (out->covariance->matrix)
        gsl_matrix_free(out->covariance->matrix);
    out->covariance->matrix	= cov;
}

static apop_model * apop_estimate_OLS(apop_data *inset, apop_model *ep){
    apop_assert(inset,  NULL, 0,'s', "You asked me to estimate a regression with NULL data.");
  apop_data         *set;
  ep->status = 0;
    apop_ls_settings   *olp =  apop_settings_get_group(ep, "apop_ls");
    if (!olp) 
        olp = Apop_model_add_group(ep, apop_ls);
    ep->data = inset;
    set = olp->destroy_data ? inset : apop_data_copy(inset); 
    
    gsl_vector *weights    = olp->destroy_data      //this may be NULL.
                                ? ep->data->weights 
                                : apop_vector_copy(ep->data->weights);
    if (apop_strcmp(ep->name, "Weighted Least Squares") && weights)
        for (size_t i =0; i< weights->size; i++)
            gsl_vector_set(weights, i, sqrt(gsl_vector_get(weights, i)));

  gsl_vector *y_data     = apop_vector_copy(set->vector);
  gsl_vector *xpy        = gsl_vector_calloc(set->matrix->size2);
  gsl_matrix *xpx        = gsl_matrix_calloc(set->matrix->size2, set->matrix->size2);
    if (olp->want_expected_value)
        ep->expected   = apop_data_alloc(0, set->matrix->size1, 3);
    if (olp->want_cov=='y')
        ep->covariance = apop_data_alloc(0, set->matrix->size2, set->matrix->size2);
    prep_names(ep);
    if (weights){
        gsl_vector_mul(y_data, weights);
        for (size_t i = 0; i < set->matrix->size2; i++){
            APOP_COL(set, i, v);
            gsl_vector_mul(v, weights);
        }
    }

    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set->matrix, set->matrix, 0, xpx);   //(X'X)
    gsl_blas_dgemv(CblasTrans, 1, set->matrix, y_data, 0, xpy);       //(X'y)
    xpxinvxpy(set->matrix, y_data, xpx, xpy, ep);
    gsl_vector_free(y_data); 
    gsl_matrix_free(xpx);
    gsl_vector_free(xpy);

    ep->llikelihood  = ols_log_likelihood(ep->data, ep);
    if (!olp->destroy_data)
        apop_data_free(set);
    if (olp->want_cov == 'y')
        apop_estimate_parameter_t_tests(ep);
    ep->status       = 1;
    return ep;
}

apop_data *ols_predict (apop_data *in, apop_model *m){

    if (!in->vector) //in->vector = gsl_vector_alloc(in->matrix->size1);
        ols_shuffle(in);

    //OK, data is now in the right form.
    //find x dot y
    gsl_blas_dgemv (CblasNoTrans, 1, in->matrix, m->parameters->vector, 0, in->vector);
    return in;
}

apop_model apop_ols = {.name="Ordinary Least Squares", .vbase = -1, .estimate =apop_estimate_OLS, 
            .log_likelihood = ols_log_likelihood, .score=ols_score, .prep = ols_prep, .predict=ols_predict};

apop_model apop_wls = {"Weighted Least Squares", .vbase = -1, .estimate = apop_estimate_OLS, 
            .log_likelihood = ols_log_likelihood, .score=ols_score, .prep= ols_prep};


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
        apop_assert(0, 0,  0, 's', "Your instrument matrix has data, but neither a vector element nor row names indicating what columns in the original data should be replaced.");
    return out;
}


static apop_model * apop_estimate_IV(apop_data *inset, apop_model *ep){
  apop_assert(inset, NULL, 0,'s', "You asked me to estimate a regression with NULL data.");
  apop_data         *set, *z;
  int               i;
    apop_ls_settings   *olp =  apop_settings_get_group(ep, "apop_ls");
    if (!olp) 
        olp = Apop_model_add_group(ep, apop_ls);
    olp->want_cov       = 'n';//not working yet.
    if (!olp->instruments || !olp->instruments->matrix->size2) 
        return apop_estimate(inset, apop_ols);
    ep->data = inset;
    if(ep->parameters)
        apop_data_free(ep->parameters);
    ep->parameters = apop_data_alloc(inset->matrix->size2,0,0);
    set = olp->destroy_data ? inset : apop_data_copy(inset); 
    z   = prep_z(inset, olp->instruments);
    
    gsl_vector *weights = olp->destroy_data      //the weights may be NULL.
                             ? ep->data->weights 
                             : apop_vector_copy(ep->data->weights);
    if (weights)
        for (i =0; i< weights->size; i++)
            gsl_vector_set(weights, i, sqrt(gsl_vector_get(weights, i)));

  apop_data    *y_data     = apop_data_alloc(set->matrix->size1, 0, 0); 
    if (olp->want_expected_value)
        ep->expected   = apop_data_alloc(0, set->matrix->size1, 3);
    if (olp->want_cov=='y')
        ep->covariance = apop_data_alloc(0, set->matrix->size1, set->matrix->size1);
    prep_names(ep);
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
    ep->parameters = apop_dot(zpxinv, zpy, 0);
    apop_data_free(y_data);
    apop_data_free(zpx); 
    apop_data_free(zpxinv);
    apop_data_free(zpy);

    if (!olp->destroy_data)
        apop_data_free(set);
//    apop_estimate_parameter_t_tests(epout);
    ep->status   = 1;
    return ep;
}

apop_model apop_iv = {.name="instrumental variables", .vbase = -1, .estimate =apop_estimate_IV, .log_likelihood = ols_log_likelihood};
