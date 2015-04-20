#include <apop.h>
#define Diff(left, right, eps) Apop_stopif(fabs((left)-(right))>(eps), abort(), 0, "%g is too different from %g (abitrary limit=%g).", (double)(left), (double)(right), eps)


/* The entropy function, like some other functions (including apop_update) has a lookup
 table for known models like the Normal distribution. If the input model has \c log_likelihood, \c p, and \c draw functions that are the ones found in \ref apop_nomrmal, then use a known calculation to report entropy; else report based on random draws from the model.

If we make a copy of the \ref apop_normal model and replace the log likelihood with a new function that produces identical values, the lookup table will not find the modified model, and the calculation via random draws will be done. Of course, the final entropy as calculated using both methods should differ only by a small amount.
*/
long double mask(apop_data *d, apop_model *m){
    return apop_normal->log_likelihood(d, m);
}

int main(){
    for (double i=0.1; i< 10; i+=.2){
        apop_model *n = apop_model_set_parameters(apop_normal, 8, i);
        long double v= apop_model_entropy(n);
        n->log_likelihood = mask;
        long double w= apop_model_entropy(n, 50000);
        Diff(v, w, 5e-2);
    }
}
