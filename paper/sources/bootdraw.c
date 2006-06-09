#include <gsl/gsl_rng.h>

void draw_wo_replacement(gsl_vector *data, int draws, gsl_vector *out_list){
const int   draws = 2000;
int    i,j,draw, ok, drawn_indices[draws];
gsl_rng *rn = apop_rng_alloc(7);

//Keep drawing values until we get enough draws.
      for(i=0;i< draws; i++){
            draw   = (int) gsl_rng_uniform_int(rn, data->size);
            ok     = 1;
            for(j=0; j<i; j++)
                  if(drawn_indices[j]==draw){
                        //then throw out the draw
                        i--;
                        ok--;
                        break;  //once we find what we want, quit the loop.
                  }
            if (ok){
                  gsl_vector_set(out_list, i, gsl_vector_get(data, draw));
                  drawn_indices[i] = draw;
            }
}
