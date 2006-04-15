#include <apophenia/headers.h>

apop_model econ_101;

static apop_estimate * econ101_estimate(apop_data * dummy, void *parameters){
apop_estimate           *est;
apop_estimation_params  *mle_params = apop_estimation_params_alloc();
apop_data               *p  = apop_data_from_vector(parameters);
    mle_params->method       = 500;
    mle_params->starting_pt  = NULL;
    mle_params->step_size    = 1e-1;
    mle_params->tolerance    = 1e-8;
    mle_params->verbose      = 0;
    apop_inventory_filter(&(mle_params->uses), econ_101.inventory_filter);
    est = apop_maximum_likelihood(p, econ_101, mle_params);
    est->estimation_params.uses.covariance    = 0;    //MLE finds them, but it's meaningless.
    est->estimation_params.uses.confidence    = 0;    //MLE finds them, but it's meaningless.
    est->estimation_params.uses.log_likelihood  = 0;
    return est;
}

/* The constraint function, including three constraints: x_0>0, x_1>0, and
 the bundle is under budget.  First, we check wether anything binds,
 and if not we return zero immediately.  Both sets of constraints are
 handled in the same way: derive new values that are within bounds,
 and report how far you had to move.
*/

static double budget_constraint(gsl_vector *beta, void * d, 
                                        gsl_vector *returned_beta){
gsl_matrix   *budget = d;
double  price0      = gsl_matrix_get(budget, 0, 0),
        price1      = gsl_matrix_get(budget, 1, 0),
        cash        = gsl_matrix_get(budget, 2, 0),
        x0          = gsl_vector_get(beta, 0),
        x1          = gsl_vector_get(beta, 1),
        tolerance   = 1e-3,
        new_x0, new_x1, penalty;
    if ((x0 * price0 + x1 * price1<= cash) && (x0 > 0) && (x1 > 0))
        return 0;
    //else:
    if (x0 <= 0 || x1 <= 0){
        new_x0  = (new_x0 > 0)? x0 : tolerance;
        new_x1  = (new_x1 > 0)? x1 : tolerance;
        penalty = (fabs(new_x0 - x0) + fabs(new_x1 + x1));
    } else {
        new_x0  = GSL_MAX(0, cash - x1* price1);
        new_x1  = GSL_MIN(x1, cash - new_x0* price0);
        penalty = (GSL_MAX(0,x0 * price0 + x1 * price1- cash));
    }
    gsl_vector_set(returned_beta, 0, new_x0);
    gsl_vector_set(returned_beta, 1, new_x1);
    return penalty;
}

static double econ101_log_likelihood(const gsl_vector *beta, void *d){
gsl_matrix  *params = d;
double      bb0     = gsl_matrix_get(params, 3,0),
            bb1     = gsl_matrix_get(params, 4,0),
            qty0    = gsl_vector_get(beta, 0),
            qty1    = gsl_vector_get(beta, 1);
    return pow(qty0, bb0) * pow(qty1, bb1);
}    

apop_model econ_101 = {"Max Cobb-Douglass subject to a budget constraint", 2,  {
    1,    //parameters
    0,    //cov
    0,    //confidence
    0,    //predicted
    0,    //residuals
    0,    //log_likelihood
    0    //names;
},         
    econ101_estimate, econ101_log_likelihood, NULL, NULL, budget_constraint, NULL};

int main(){
double          param_array[]   =  {1, 3, 38.4, 0.4, 0.6};//see header.
gsl_vector      *params         = apop_array_to_vector(param_array,5);
apop_data       *params_again   = apop_data_from_vector(params);
gsl_vector      *marginals      = gsl_vector_alloc(2);
apop_estimate   *e;
double          x1, x2;
    e   = econ_101.estimate(NULL,  params);
    printf("The optimal quantities:\n");
    apop_estimate_print(e);

    x2  = param_array[2]/(param_array[3]/param_array[4]+1)/param_array[1];
    x1  = (param_array[2] - param_array[1] * x2)/param_array[0];
    printf("\nAnalytically, these should be:\n %g\n %g\n\n", x1, x2);

    apop_fn_for_derivative  = econ_101.log_likelihood;
    apop_numerical_gradient(e->parameters, params_again->matrix, marginals);
    printf("The marginal values:\n");
    sprintf(apop_opts.output_delimiter, "\n");
    apop_vector_show(marginals);
    printf("\nAnalytically, these should be:\n %g\n %g\n", 
            param_array[3]*pow(x1,param_array[3]-1)*pow(x2,param_array[4]),
            param_array[4]*pow(x2,param_array[4]-1)*pow(x1,param_array[3]));
    return 0;
}
