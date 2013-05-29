check_PROGRAMS= test_apop distribution_tests error_test nist_tests rake_test \
../eg/apop_map_row ../eg/binning ../eg/db_fns ../eg/dot_products  			 \
../eg/fake_logit ../eg/fix_params ../eg/iv ../eg/ks_tests ../eg/jacobian ../eg/normalizations \
../eg/pmf_test ../eg/some_cdfs ../eg/sorting ../eg/stack_models 		 \
../eg/test_distances ../eg/test_fisher ../eg/test_harmonic ../eg/test_kl_divergence \
../eg/test_pruning ../eg/test_ranks ../eg/test_regex ../eg/test_strip_dots 	 \
../eg/test_updating
TESTS=$(check_PROGRAMS)

LDADD=../libapophenia.la
AM_CFLAGS = -DTesting $(CFLAGS) -I$(top_build_prefix)/$(top_builddir) 

#AM_TESTS_ENVIRONMENT= cd $(srcdir)/tests &&
EXTRA_DIST = test.py Readme numacc4.dat pontius.dat test_data test_data2 \
		data-mixed test_data_nans wampler1.dat printing_sample test_data_fixed_width
