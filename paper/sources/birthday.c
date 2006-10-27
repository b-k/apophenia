#include <math.h>
#include <stdio.h>

int main(){
  double    no_match    = 1;
  double    matches_me;
  int       ct;
    printf("people\tmatches me\tany match\n");
    for (ct=2; ct<=40; ct ++){
        matches_me   = 1- pow(364/365., ct-1);
        no_match    *= (1 - (ct-1)/365.);
        printf("%i\t%.3g%%\t\t%.3g%%\n", ct, matches_me*100, 100*(1-no_match));
    }
    return 0;
}
