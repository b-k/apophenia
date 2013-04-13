check_PROGRAMS= ../eg/stack_models ../eg/dot_products ../eg/binning ../eg/test_ranks \
		../eg/test_harmonic ../eg/test_distances test_apop rake_test error_test \
		distribution_tests nist_tests ../eg/some_cdfs ../eg/test_fisher \
		../eg/test_strip_dots ../eg/test_kl_divergence ../eg/iv ../eg/pmf_test \
		../eg/test_pruning ../eg/test_updating ../eg/apop_map_row ../eg/test_regex \
		../eg/normalizations ../eg/fake_logit
TESTS=$(check_PROGRAMS)

LDADD=../libapophenia.la
AM_CFLAGS = -DTesting $(CFLAGS) -I$(top_build_prefix)/$(top_builddir) 

#AM_TESTS_ENVIRONMENT= cd $(srcdir)/tests &&
EXTRA_DIST = test.py Readme numacc4.dat pontius.dat test_data test_data2 \
		data-mixed test_data_nans wampler1.dat printing_sample test_data_fixed_width
