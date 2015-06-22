#include <apop.h>

void with_fixed_numbers(){
    apop_data *a =apop_data_alloc(2,2,2);
    double    eight   = 8.0;
    apop_data_fill(a, 8, 2.2, eight/2,
                      0, 6.0, eight);
    apop_data_show(a);
}

void with_a_list(){
  apop_data *a =apop_data_alloc(2,2,2);
  double    eight   = 8.0;
  double list[] = {8, 2.2, eight/2,
                   0, 6.0, eight};
    apop_data_fill_base(a, list);
    apop_data_show(a);
}

int main(){
    with_fixed_numbers();
    printf("-----\n");
    with_a_list();
}
