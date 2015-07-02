#include <apop.h>

/* We promise users that they can copy factors from one data set to the next. But if the
second data set has values that do not appear in the first data set, then they have
to be added to the factor list.
*/

int main(){
    apop_data *d1 = apop_text_alloc(NULL, 5, 1);
    apop_data *d2 = apop_text_alloc(NULL, 5, 1);
    apop_text_fill(d1, "A", "B", "C", "B", "B");
    apop_text_fill(d2, "B", "A", "D", "B", "B");
    apop_data_to_factors(d1);
    apop_data_show(d1);
    d2->more = apop_data_copy(apop_data_get_factor_names(d1, 0, 't'));
    printf("-----\n");
    apop_data_to_dummies(d2, .append='y');
    apop_data_show(d2);


    //some spot checks.
    assert(apop_data_get(d1, 2)==2);
    assert(apop_data_get(d2, 0, 0)==1);
    assert(apop_data_get(d2, 2, 0)==0);
    assert(apop_data_get(d2, 2, 1)==0);
    assert(apop_data_get(d2, 2, 2)==1);
    assert(apop_data_get(d2, 3, 0)==1);
}
