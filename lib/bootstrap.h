
gsl_vector * bootstrap(gsl_matrix * data, gsl_vector *(*boot_fn)(gsl_matrix *, void *, void* , void*), int boot_iterations,
	void *params_1, void *params_2, void *params_3);
