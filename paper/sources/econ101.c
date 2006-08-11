#include <apophenia/headers.h>

apop_model econ_101;

static apop_estimate * econ101_estimate(apop_data * dummy, void *parameters){
apop_estimate           *est;
apop_estimation_params  *mle_params = apop_estimation_params_alloc();
apop_data               *p  = apop_vector_to_data((gsl_vector*)parameters);
double              start[] = {1,1};
    mle_params->method       = 000;
    mle_params->starting_pt  = start;
    mle_params->step_size    = 1e-1;
    mle_params->tolerance    = 1e-8;
    mle_params->verbose      = 0;
    apop_inventory_set(&mle_params->uses, 0);
    mle_params->uses.parameters    = 1;
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
apop_data   *budget = d;
double  price0      = apop_data_get(budget, 0, -1),
        price1      = apop_data_get(budget, 1, -1),
        cash        = apop_data_get(budget, 2, -1),
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

static double econ101_log_likelihood(const gsl_vector *beta, apop_data *d){
double      bb0     = apop_data_get(d, 3,-1),
            bb1     = apop_data_get(d, 4,-1),
            qty0    = gsl_vector_get(beta, 0),
            qty1    = gsl_vector_get(beta, 1);
    return pow(qty0, bb0) * pow(qty1, bb1);
}    

apop_model econ_101 = {"Max Cobb-Douglass subject to a budget constraint", 2,  
    econ101_estimate, econ101_log_likelihood, NULL, NULL, budget_constraint, NULL};

