check_PROGRAMS= db_tests distribution_tests error_test nist_tests rake_test sort_example test_apop update_via_rng \
../eg/apop_map_row ../eg/binning ../eg/dconstrain ../eg/draw_to_db ../eg/db_fns ../eg/dot_products  			 \
../eg/fake_logit ../eg/fix_params ../eg/hills2 ../eg/iv ../eg/jacobian ../eg/ks_tests ../eg/normalizations \
../eg/parameterization ../eg/pmf_test ../eg/simple_subsets ../eg/some_cdfs ../eg/stack_models 		 \
../eg/test_distances ../eg/test_fisher ../eg/test_harmonic ../eg/test_kl_divergence \
../eg/test_pruning ../eg/test_ranks ../eg/test_regex ../eg/test_updating

TESTS=utilities_test $(check_PROGRAMS)

LDADD=../libapophenia.la
#AM_CFLAGS = -DTesting $(CFLAGS) -I$(top_build_prefix)/$(top_builddir) -DDatadir=\"$(srcdir)/\"
AM_CFLAGS = -DTesting $(CFLAGS) -DDatadir=\"$(srcdir)/\"

EXTRA_DIST = utilities_test Readme data numacc4.dat pontius.dat amash_vote_analysis.csv faith.data test_data test_data2 \
		data-mixed test_data_nans wampler1.dat printing_sample test_data_fixed_width sort_tests.c
