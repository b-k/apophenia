ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTIONS = \
	dist-xz \
	dist-bzip2 \
	dist-zip

AM_CFLAGS = -g -Wall -O3

## Library versioning (C:R:A == current:revision:age)
## 0.999b 0:0:0
## 0.999c 1:0:0
LIBAPOPHENIA_LT_VERSION = 1:0:0

SUBDIRS = transform model . cmd tests docs eg


include_HEADERS = apop.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA= apophenia.pc

lib_LTLIBRARIES = libapophenia.la

libapophenia_la_LD_VERSION_SCRIPT=
if HAVE_LD_VERSION_SCRIPT
libapophenia_la_LD_VERSION_SCRIPT+= -Wl,--version-script=$(top_srcdir)/apophenia.map
endif

SUBLIBS = \
	libapopkernel.la \
	transform/libapoptransform.la \
	model/libapopmodel.la

libapophenia_la_SOURCES = \
	asprintf.c

noinst_LTLIBRARIES = libapopkernel.la

noinst_HEADERS = apop_internal.h

libapopkernel_la_SOURCES = \
	apop_arms.c \
	apop_asst.c \
	apop_bootstrap.c \
	apop_conversions.c \
	apop_data.c \
	apop_db.c \
	apop_fexact.c \
	apop_hist.c \
	apop_linear_algebra.c \
	apop_linear_constraint.c \
	apop_mapply.c \
	apop_mcmc.c \
	apop_missing_data.c \
	apop_mle.c apop_model.c \
	apop_name.c \
	apop_output.c \
	apop_rake.c \
	apop_regression.c \
	apop_settings.c \
	apop_smoothing.c \
	apop_sort.c \
	apop_stats.c \
	apop_tests.c \
	apop_update.c	\
	apop_vtables.c

apop_db_INCLUDES = \
	apop_db_mysql.c \
	apop_db_sqlite.c

apop_db.c: $(apop_db_INCLUDES)

libapopkernel_la_CFLAGS = \
	$(PTHREAD_CFLAGS) \
	$(OPENMP_CFLAGS) \
	$(MYSQL_CFLAGS) \
	$(SQLITE3_CFLAGS) \
	$(GSL_CFLAGS)

libapophenia_la_LDFLAGS = \
	-version-info $(LIBAPOPHENIA_LT_VERSION) \
	$(libapophenia_la_LD_VERSION_SCRIPT)

libapophenia_la_LIBADD = \
	$(SUBLIBS) \
	$(MYSQL_LDFLAGS) \
	$(SQLITE3_LDFLAGS) \
	$(GSL_LIBS) \
	$(PTHREAD_LIBS) \
	$(LIBM)

EXTRA_DIST = \
	COPYING2 \
	rpm.spec \
	apophenia.pc.in \
	apophenia.map

EXTRA_DIST += \
	$(apop_db_INCLUDES)

## compatibility
doc:
	-$(MAKE) -C docs doc
