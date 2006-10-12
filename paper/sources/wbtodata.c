#include <apophenia/headers.h>

void calc_gdp_per_cap(apop_data *d){
  gsl_vector v0  = gsl_matrix_column(d->matrix, 0).vector;
  gsl_vector v1  = gsl_matrix_column(d->matrix, 1).vector;
    d->vector   = gsl_vector_alloc((&v0)->size);
    gsl_vector_memcpy(d->vector, &v1);
    gsl_vector_div(d->vector, &v0);
    apop_name_add(d->names, "GDP per cap", 'v');
}

int main(){
  apop_data *d;
    apop_db_open("data-wb.db");
    strcpy(apop_opts.db_name_column, "country");
    d   = apop_query_to_data("select pop.country as country, \
        pop.population as pop, gdp.GDP as GDP\
        from pop, gdp \
        where pop.country == gdp.country");
    calc_gdp_per_cap(d);
    apop_data_show(d);
    return 0;
}
