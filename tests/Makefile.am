
if EXTENDED_TESTS
EXTRA_TESTS = distribution_tests \
	lognormal_test \
	rake_test \
	test_kernel_ll \
	update_via_rng \
	$(top_builddir)/eg/cross_models \
	$(top_builddir)/eg/dconstrain \
	$(top_builddir)/eg/entropy_model \
	$(top_builddir)/eg/f_test \
	$(top_builddir)/eg/fix_params \
	$(top_builddir)/eg/hills2 \
	$(top_builddir)/eg/jack \
	$(top_builddir)/eg/jacobian \
	$(top_builddir)/eg/ml_imputation \
	$(top_builddir)/eg/pmf_test \
	$(top_builddir)/eg/some_cdfs \
	$(top_builddir)/eg/test_kl_divergence \
	$(top_builddir)/eg/test_ranks \
	$(top_builddir)/eg/test_updating \
	$(top_builddir)/eg/transform
else
	EXTRA_TESTS =
endif

check_PROGRAMS= \
	db_tests \
	error_test \
	nist_tests \
	sort_example \
	test_apop \
	$(top_builddir)/eg/apop_map_row \
	$(top_builddir)/eg/binning \
	$(top_builddir)/eg/draw_to_db \
	$(top_builddir)/eg/db_fns \
	$(top_builddir)/eg/dot_products \
	$(top_builddir)/eg/entropy_vector \
	$(top_builddir)/eg/fake_logit \
	$(top_builddir)/eg/iv \
	$(top_builddir)/eg/ks_tests \
	$(top_builddir)/eg/parameterization \
	$(top_builddir)/eg/simple_subsets \
	$(top_builddir)/eg/test_distances \
	$(top_builddir)/eg/test_fisher \
	$(top_builddir)/eg/test_harmonic \
	$(top_builddir)/eg/test_pruning \
	$(top_builddir)/eg/test_regex $(EXTRA_TESTS)

TESTS = \
	utilities_test \
	$(check_PROGRAMS)

AM_CFLAGS = \
	-DTesting \
	-DDatadir=\"$(top_srcdir)/tests/\" \
	-I$(top_srcdir)/tests \
	-I$(top_srcdir) \
	$(GSL_CFLAGS)

AM_LDFLAGS = \
	$(top_builddir)/libapophenia.la \
	$(GSL_LIBS)

DATA_DIST = \
	data \
	data-mixed \
	printing_sample \
	test_data \
	test_data2 \
	test_data_nans\
	test_data_fixed_width \
	amash_vote_analysis.csv \
	numacc4.dat \
	pontius.dat \
	wampler1.dat \
	faith.data \
	sort_tests.c

EXTRA_DIST = \
	Readme \
	$(DATA_DIST)

CLEANFILES = \
	ff.db \
	runs.db \
	td.db \
	draws-k \
	draws-k2 \
	draws-mvN \
	draws-N \
	draws-std_multinormal \
	draws-std_normal \
	the_data.txt \
	print_test.out \
	xxx
