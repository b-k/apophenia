#include <apophenia/headers.h>

typedef double (*dfn) (double);

char gnuplot[] = "/usr/bin/gnuplot -persist";

double sample_function (double in){
    return log(in)+ sin(in);
}

void plot_a_fn(double min, double max, dfn plotme){
  double    i;
  double    val;
    FILE *f = popen(gnuplot, "w");
    if (!f)
        printf("Couldn't find Gnuplot.\n");
    fprintf(f, "set key off\n plot '-' with lines\n");
    for (i=min; i<max; i+= (max-min)/100.0){
        val = plotme(i);
        fprintf(f, "%g\t%g\n", i, val);
    }
    fclose(f);
}

int main(){
    plot_a_fn(0, 15, sample_function);
    return 0;
}
