#include <apophenia/headers.h>

int main(){
  double    sigma_min   = 0.1;
  double    sigma_max   = 1;
  double    sigma_incr  = 0.1;
  double    val_min     = 0.1;
  double    val_max     = 2.1;
  double    val_incr    = 0.1;
  double sigma,x;
    for (sigma=sigma_min; sigma< sigma_max; sigma+=sigma_incr*2)
        printf("\t");
    printf("\t    Sigma\n\t");
    for (sigma=sigma_min; sigma< sigma_max; sigma+=sigma_incr)
        printf("\t%4g", sigma);
    printf("\nValue\t\t");
    for (sigma=sigma_min; sigma< sigma_max; sigma+=sigma_incr)
        printf("--------");
    printf("\n");
    for (x=val_min; x< val_max; x+=val_incr){
        printf("%.1f\t|\t", x);
        for (sigma=sigma_min; sigma< sigma_max; sigma+=sigma_incr)
            printf("%.4f\t",  gsl_cdf_gaussian_P(x,sigma));
        printf("\n");
    }
    return 0;
}
