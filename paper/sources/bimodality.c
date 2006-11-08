#include <apophenia/headers.h>
double  modect_scale, modect_min, modect_max,
        boot_iterations = 400,
        max_k           = 30,
        pauselength     = 0.4,
        resolution      = 400,
        fixed_max_h     = 1e6,
        tolerance       = 1e-2;
char    outfile[]   = "kernelplot";
        infile[]    = "men.1834.nosport.txt";
gsl_rng *r;

void plot(apop_data *d){
  FILE  *f  = fopen(outfile, "a");
    apop_opts.output_type   = 'p';
    apop_opts.output_pipe   = f;
    fprintf(f, "plot '-' with lines\n");
    apop_data_print(d, NULL);
    fprintf(f, "e\npause %g\n", pauselength);
    fclose(f);
}

int countmodes(gsl_vector *data, double h, int plotme){
  int       len     =(modect_max-modect_min)/modect_scale +1;
  apop_data *fmap   = apop_data_calloc(len, 2);
  double    t, sum;
  int       k, j    = 0;
  int     modect    =0;
	for (t = modect_min; t<modect_max; t+=modect_scale){
        apop_data_set(fmap, j, 0, t);
        sum = 0;
	    for (k = 0; k< data->size; k++)
            sum += gsl_ran_gaussian_pdf((t-gsl_vector_get(data,k))/h,1)/(data->size*h);
        apop_data_set(fmap, j++, 1, sum);
    }
	for (t = 2; t< len; t++){
        if(apop_data_get(fmap,t,1)<=apop_data_get(fmap,t-1,1) 
               && apop_data_get(fmap,t-1,1)>apop_data_get(fmap,t-2,1))
            modect++;
    }
    if (plotme) plot(fmap);
    apop_data_free(fmap);
    return modect;
}

double find_smallest_h_for_k(gsl_vector *data, int goal_k, double min, double max){
  double try_me = min + (max-min)/2;
  int   modect  = countmodes(data, try_me, 0);
    if (modect > goal_k)  //h is too small
        return find_smallest_h_for_k(data, goal_k, try_me, max);
    else {
        if (modect==goal_k && (max-min)/2 < tolerance)
            return try_me;    //stopping condition here.
        else
            return find_smallest_h_for_k(data, goal_k, min, try_me);
    }
}

gsl_vector * produce_k_to_h_table(gsl_vector *data){
  gsl_vector *out    = gsl_vector_alloc(max_k-1);
  int        i;
  double     val     = find_smallest_h_for_k(data, 1, 0, fixed_max_h);
    gsl_vector_set(out, 0, val);
    for(i=2; i< max_k; i++){
        val     = find_smallest_h_for_k(data, i, 0, gsl_vector_get(out, i-2));
        gsl_vector_set(out, i-1, val);
    }
    return out;
}

apop_data *produce_p_table(gsl_vector *data,gsl_vector *ktohmap){
  apop_data     *ptab = apop_data_alloc(max_k-1, 2);
  int           el, i, j, ct  = data->size;
  gsl_vector    *boots = gsl_vector_alloc(ct);
    for (el=0; el< max_k-1; el++)
        apop_data_set(ptab, el, 0, el+1);
    for (i=0; i<boot_iterations; i++){
        for (j=0; j< ct; j++)
            gsl_vector_set(boots, j, gsl_vector_get(data,gsl_rng_uniform(r) *ct));
        for (el=0; el< max_k-1; el++)
            if(countmodes(boots,gsl_vector_get(ktohmap,el), 0)<=el+1)
                apop_matrix_increment(ptab->matrix, el, 1, 1./boot_iterations);
    }
    gsl_vector_free(boots);
    return ptab;
}

void setvars(gsl_vector *data){
    double m1=gsl_vector_max(data);
    double m2=gsl_vector_min(data);
    modect_scale=(m1-m2)/resolution;
    modect_min=m2-(m1-m2)/10;
    modect_max=m1+(m1-m2)/10;
}

int main(){
  apop_data     *set, *p;
  gsl_vector    *m;
  double        h, min, max;
    r   = apop_rng_alloc(342);
    remove(outfile);
    set = apop_text_to_data(infile, 0,0);
    APOP_COL(set, 0, data);
    setvars(data);
    m = produce_k_to_h_table(data);

    min =gsl_vector_get(m, max_k-2) *0.88;
    max =gsl_vector_get(m, 0);
    for (h= min; h< max; h+= (max-min)/100)
        countmodes(data, h, 1); //just for plotting.

    printf("the confidence of rejecting a mode count:\n");
    p = produce_p_table(data,m);
    apop_data_show(p);
    return 0;
}
