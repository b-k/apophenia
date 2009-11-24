#include <apop.h>

int main(){
apop_data   *data, *data_copy;
apop_model  *est;
FILE        *f;
char        outfile[]   = "scatter.gplot";

    apop_db_open("data-metro.db");
    data = apop_query_to_data("select riders, year from riders where station like 'Silver%%'");
    apop_db_close();

    //The regression destroys your data, so copy it first.
    data_copy   = apop_data_copy(data);

    //Run OLS, display results on terminal
    est  = apop_estimate(data, apop_OLS);
    apop_model_show(est);

    //Prep the file with a header, then call the function.
    f    = fopen(outfile, "w");
    fprintf(f,"set term postscript;\n set output \"scatter.eps\"\n set yrange [0:*]\n");
    apop_plot_line_and_scatter(data_copy, est, .output_pipe=f);
    fclose(f);

    return 0;
}
