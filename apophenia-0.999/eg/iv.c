/* Instrumental variables are often used to deal with variables measured with noise, so
this example produces a data set with a column of noisy data, and a separate instrument
measured with greater precision, then sets up and runs an instrumental variable regression.

To guarantee that the base data set has noise and the instrument is cleaner, the
procedure first generates the clean data set, then copies the first column to the
instrument set, then the \c add_noise function inserts Gaussian noise into the base
data set. Once the base set and the instrument set have been generated, the setup for
the IV consists of adding the relevant names and using \ref Apop_model_add_group to add a
\c lm (linear model) settings group with a <tt>.instrument= instrument_data</tt> element.

In fact, the example sets up a sequence of IV regressions, with more noise each
time. This sample is part of Apophenia's test suite, and so checks that the coefficients
are correct along the way.
*/

#include <apop.h>
#define Diff(L, R, eps) Apop_stopif(fabs((L)-(R)>=(eps)), return, 0, "%g is too different from %g (abitrary limit=%g).", (double)(L), (double)(R), eps);

int datalen =1e4;

//generate a vector that is the original vector + noise
void add_noise(gsl_vector *in, gsl_rng *r, double size){
    apop_model *nnoise = apop_model_set_parameters(apop_normal, 0, size);
    for (int i=0; i< in->size; i++){
        double noise;
        apop_draw(&noise, r, nnoise);
        *gsl_vector_ptr(in, i) += noise;
    }
    apop_model_free(nnoise);
}

void test_for_unbiased_parameter_estimates(apop_model *m, double tolerance){
        Diff(apop_data_get(m->parameters, .row=0,.col=-1), -1.4, tolerance);
        Diff(apop_data_get(m->parameters, .row=1,.col=-1), 2.3, tolerance);
}

int main(){
    gsl_rng *r = apop_rng_alloc(234);

    apop_data *data = apop_data_alloc(datalen, 2);
    for(int i=0; i< datalen; i++){
        apop_data_set(data, i, 1, 100*(gsl_rng_uniform(r)-0.5));
        apop_data_set(data, i, 0, -1.4 + apop_data_get(data,i,1)*2.3);
    }
    apop_name_add(data->names, "dependent", 'c');
    apop_name_add(data->names, "independent", 'c');
    #ifndef Testing
    apop_model *oest = apop_estimate(data, apop_ols);
    apop_model_show(oest);
    #endif

    //the data with no noise will be the instrument.
    Apop_col_v(data, 1, col1);
    apop_data *instrument_data = apop_data_alloc(data->matrix->size1, 1);
    Apop_col_v(instrument_data, 0, firstcol);
    gsl_vector_memcpy(firstcol, col1);
    apop_name_add(instrument_data->names, "independent", 'c');
    Apop_model_add_group(apop_iv, apop_lm, .instruments = instrument_data);

    //Now add noise to the base data four times, and estimate four IVs.
    int tries = 4;
    apop_model *ests[tries];
    for (int nscale=0; nscale<tries; nscale++){
        add_noise(col1, r, nscale==0 ? 0 : pow(10, nscale-tries));
        ests[nscale] = apop_estimate(data, apop_iv);
        #ifndef Testing
        if (nscale==tries-1){ //print the one with the largest error.
            printf("\nnow IV:\n");
            apop_model_show(ests[nscale]);
        }
        #endif
    }

    /* Now test. The parameter estimates are unbiased.
       As we add more noise, the covariances expand.
       Test that the ratio of one covariance matrix to the next
       is less than one, though these are typically very much
       smaller than one (as the noise is an order of magnitude 
       larger in each case), and the ratios will be identical
       for each j, k below. */
    test_for_unbiased_parameter_estimates(ests[0], 1e-6);
    for (int i=1; i<tries; i++){
        test_for_unbiased_parameter_estimates(ests[i], 1e-3);

        gsl_matrix *cov = apop_data_get_page(ests[i-1]->parameters, "<Covariance>")->matrix;
        gsl_matrix *cov2 = apop_data_get_page(ests[i]->parameters, "<Covariance>")->matrix;
        gsl_matrix_div_elements(cov, cov2);
        for (int j =0; j< 2; j++)
            for (int k =0; k< 2; k++)
                assert(gsl_matrix_get(cov, j, k) < 1);
    }
}
