
noinst_LTLIBRARIES = libapopmodel.la

libapopmodel_la_SOURCES = \
	apop_bernoulli.c \
	apop_beta.c \
	apop_dirichlet.c \
	apop_exponential.c \
	apop_gamma.c \
	apop_kerneld.c \
	apop_loess.c \
	apop_multinomial.c \
	apop_multivariate_normal.c \
	apop_normal.c \
	apop_ols.c \
	apop_pmf.c \
	apop_poisson.c \
	apop_probit.c \
	apop_t.c \
	apop_uniform.c \
	apop_wishart.c \
	apop_yule.c \
	apop_zipf.c

libapopmodel_la_CFLAGS = \
	-I $(top_srcdir) \
	$(PTHREAD_CFLAGS) \
	$(OPENMP_CFLAGS) \
	$(SQLITE3_CFLAGS) \
	$(GSL_CFLAGS)
