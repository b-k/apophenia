/* A round trip: generate Zipf-distributed draws, summarize them to a single list of
rankings, then expand the rankings to a list of single entries. The sorted list at the end
of this should be identical to the (sorted) original list. */
#include <apop.h>

int main(){
    gsl_rng *r = apop_rng_alloc(2342);
    int i, length = 1e4;
    apop_model *a_zipf = apop_model_set_parameters(apop_zipf, 3.2);
    apop_data *draws = apop_data_alloc(length);
    for (i=0; i< length; i++)
        apop_draw(apop_data_ptr(draws, i, -1), r, a_zipf);
    apop_data *by_rankings = apop_data_rank_compress(draws);
    //The first row of the matrix is suitable for plotting.
    //apop_data_show(by_rankings);
    assert(apop_matrix_sum(by_rankings->matrix) == length);

    apop_data *re_expanded = apop_data_rank_expand(by_rankings);
    gsl_sort_vector(draws->vector);
    gsl_sort_vector(re_expanded->vector);
    assert(apop_vector_distance(draws->vector, re_expanded->vector) < 1e-5);
}
