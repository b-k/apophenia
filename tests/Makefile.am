check_PROGRAMS=test_apop rake_test error_test distribution_tests nist_tests
TESTS=$(check_PROGRAMS)

LDADD=../libapophenia.la
AM_CFLAGS = $(CFLAGS) -I$(top_build_prefix)/$(top_builddir)

#AM_TESTS_ENVIRONMENT= cd $(srcdir)/tests &&
EXTRA_DIST = test.py Readme numacc4.dat pontius.dat test_data_mixed test_data test_data2 \
		data-mixed test_data_nans wampler1.dat printing_sample test_data_fixed_width
