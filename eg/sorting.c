#include <apop.h>

#ifndef Testing
#define showdata apop_data_show(d);
#else
#define showdata 
#endif

int main(){
    apop_data *d = apop_text_to_data("test_data2", 0 ,1);
    Apop_stopif(!d, exit(1), 0, "Please link or copy test_data2 from the tests directory "
                                "to the directory in which you are running this test.");
    for (int i=0; i< d->matrix->size1; i++){
        char *rowname=NULL;
        asprintf(&rowname, "row %i", i);
        apop_name_add(d->names, rowname, 'r');
        free(rowname);
    }

    apop_name_add(d->names, "sorted ascending by first column", 'h');
    apop_data_sort(d, 0, 'a');
    for (int i=1; i< d->matrix->size2; i++){
        assert(apop_data_get(d, i,0) >= apop_data_get(d, i-1,0));
    }
    assert(apop_data_get(d, 0,1)== 32 || apop_data_get(d, 0,1)== 9);
    showdata

    apop_name_add(d->names, "\nsorted descending by first column", 'h');
    apop_data_sort(d, 0, 'd');
    for (int i=1; i< d->matrix->size2; i++){
        assert(apop_data_get(d, i,0) <= apop_data_get(d, i-1,0));
    }
    assert(apop_data_get(d, 0,1)== 55); //(35, 55)
    showdata

    apop_name_add(d->names, "\nsorted descending by second column", 'h');
    apop_data_sort(d, 1, 'd');
    for (int i=1; i< d->matrix->size2; i++){
        assert(apop_data_get(d, i,1) <= apop_data_get(d, i-1,1));
    }
    assert(apop_data_get(d, 0,0) == 19); //(19, 60)
    showdata
}
