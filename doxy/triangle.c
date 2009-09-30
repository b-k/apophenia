#include <apop.h>

int main(){
    apop_data *d = apop_data_alloc(0, 8, 3);
    apop_data_fill(d,   1, 0, 0,
                        .8, .1, 0,
                        .9, 0, .1,
                        12, 4, 1,
                        0, 1, 0,
                        1, 2, 2,
                        2, 1, 2,
                        2, 2, 1,
                        );
    apop_name_add(d->names, "first", 'c');
    apop_name_add(d->names, "second", 'c');
    apop_name_add(d->names, "third", 'c');
    apop_plot_triangle(d, "out.gnup");
}
