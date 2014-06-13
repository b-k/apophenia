#include <apop.h>

int main(){
    double out = apop_generalized_harmonic(270, 0.0);
	assert (out == 270);
	out	= apop_generalized_harmonic(370, -1.0);
	assert (out == 370*371/2);
	out	= apop_generalized_harmonic(12, -1.0);
	assert (out == 12*13/2);
}
