#ifdef Datadir
#define DATADIR Datadir
#else
#define DATADIR "."
#endif

#include <apop.h>
#include <unistd.h>
#ifdef Testing
#include "sort_tests.c" //For Apophenia's test suite, some tedious checks that the sorts worked
#endif

//get_distance is for the sort-by-Euclidian distance example below.
double get_distance(gsl_vector *v) {return apop_vector_distance(v);}

int main(){
    apop_text_to_db( DATADIR "/" "amash_vote_analysis.csv" );
    apop_data *d = apop_query_to_mixed_data("mntmtm", "select 1,id,party,contribs/1000.0,vote,ideology from amash_vote_analysis " );

    //use the default order of columns for sorting
    apop_data *sorted = apop_data_sort(d, .inplace='n');
#ifndef Testing
    apop_data_print(sorted);
#else
    check_sorting1(sorted);
#endif

    //set up a specific column order
    apop_data *perm = apop_data_copy(Apop_r(d, 0));
    perm->vector = NULL;
    apop_data_fill(perm, 5, 3, 4);
    apop_text_set(perm, 0, 0, "2");
    apop_text_set(perm, 0, 1, "1");

    apop_data_sort(d, perm);
#ifndef Testing
    apop_data_print(d);
#else
    check_sorting2(d);
#endif

    //sort a list of names
    apop_data *blank = apop_data_alloc();
    apop_data_add_names(blank, 'r', "C", "E", "A");
    apop_data_sort(blank);
    assert(*blank->names->row[0] == 'A');
    assert(*blank->names->row[1] == 'C');
    assert(*blank->names->row[2] == 'E');

    //take each row of the matrix as a vector; store the Euclidian distance to the origin in the vector;
    //sort in descending order.
    apop_data *rowvectors = apop_text_to_data( DATADIR "/" "test_data" );
    apop_map(rowvectors, .fn_v=get_distance, .part='r', .inplace='y');
    apop_data *arow = apop_data_copy(Apop_r(rowvectors, 0));
    arow->matrix=NULL; //sort only by the distance vector
    apop_data_sort(rowvectors, arow, .asc='d');
#ifndef Testing
    apop_data_print(rowvectors);
#else
    double prev = INFINITY;
    for (int i=0; i< rowvectors->vector->size; i++){
        double this = apop_data_get(rowvectors, i, -1);
        assert(this < prev);
        prev = this;
    }
#endif
}
