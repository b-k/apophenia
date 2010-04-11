/* This program draws ten random data points, and then produces two kernel density
estimates: one based on the Normal distribution and one based on the Uniform.

It produces three outputs:
--stderr shows the random draws
--kerneldata is a file written with plot data for both KDEs
--stdout shows instructions to gnuplot, so you can pipe:
./kernel | gnuplot -persist

Most of the code is taken up by the plot() and draw_some_data() functions, which are
straightforward. Notice how plot() pulls the values of the probability distributions 
at each point along the scale.

The set_uniform_edges function sets the max and min of a Uniform distribution so that the
given point is at the center of the distribution. It will be used to center the
distribution around each point in the original data set.

The first KDE uses the defaults, which are based on a Normal distribution with std dev 1;
the second explicitly sets the .kernel and .set_fn for a Uniform.
*/

#include <apop.h>

void set_uniform_edges(apop_data_row r, apop_model *unif){
    apop_data_set(unif->parameters, 0, -1, r.matrix_row.data[0]-0.5);
    apop_data_set(unif->parameters, 1, -1, r.matrix_row.data[0]+0.5);
}

void plot(apop_model *k, apop_model *k2){
    apop_data *onept = apop_data_alloc(0,1,1);
    FILE *outtab = fopen("kerneldata", "w");
    for (double i=0; i<20; i+=0.01){
        apop_data_set(onept,0,0, i);
        fprintf(outtab, "%g %g %g\n", i, apop_p(onept, k), apop_p(onept, k2));
    }
    fclose(outtab);
    printf("plot 'kerneldata' using 1:2\n"
           "replot 'kerneldata' using 1:3\n");
}

apop_data *draw_some_data(){
    gsl_rng *r = apop_rng_alloc(7);
    apop_data *d  = apop_data_alloc(0,10,1); 
    for (int i=0; i<10; i++)
        apop_data_set(d, i, 0, gsl_rng_uniform(r)*20);
    apop_data_print(d, .output_pipe=stderr);
    return d;
}

int main(){
    apop_data *d = draw_some_data();    
    apop_model *k = apop_estimate(d, apop_kernel_density);
    apop_model *k2 = apop_model_copy(apop_kernel_density);
    apop_model_add_group(k2, apop_kernel_density, .base_data=d,
                                    .set_fn = set_uniform_edges,
                                    .kernel = &apop_uniform);
    plot(k, k2);
}
