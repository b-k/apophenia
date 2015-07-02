
noinst_LTLIBRARIES = libapoptransform.la

libapoptransform_la_SOURCES = \
	apop_dconstrain.c\
	apop_fix_params.c \
	apop_coordinate_transform.c \
	apop_mixture.c \
	apop_cross.c

libapoptransform_la_CFLAGS = \
	-I $(top_srcdir) \
	$(PTHREAD_CFLAGS) \
	$(OPENMP_CFLAGS) \
	$(SQLITE3_CFLAGS) \
	$(GSL_CFLAGS)
