#include <gsl/gsl_rng.h>

void boot_draw(gsl_vector *data, int draws, gsl_vector *out_list){
int    i,j,draw, ok, drawn_indices[draws];

//Allocate and set up a random number generator:
gsl_rng *rn;
      gsl_rng_env_setup();
      rn     =gsl_rng_alloc(gsl_rng_default);

//Keep drawing values until we get enough draws.
      for(i=0;i< draws; i++){
            draw   = (int) gsl_rng_uniform_int(rn, data->size);
            ok     = 1;
            for(j=0;j<i;j++)
                  if(drawn_indices[j]==draw){
                        //then throw out the draw
                        i--;
                        ok--;
                        break;
                  }
            if (ok){
                  gsl_vector_set(out_list, i, gsl_vector_get(data, draw));
                  drawn_indices[i] = draw;
            }
}
