
#AM_CFLAGS = -g -O3 -Wall -I./
AM_CFLAGS = -g -Wall -I./
LIBADD = -lgsl -lgslcblas -lsqlite3
#AUTOMAKE_OPTIONS = no-dependencies

#The library
LDADD = libapophenia.la
AM_LDFLAGS =  -lgsl -lgslcblas -lsqlite3
#apop_text_to_db_LDADD = libapophenia.la
#apop_text_to_db_LDFLAGS =  -lgsl -lgslcblas -lsqlite3
lib_LTLIBRARIES = libapophenia.la
libapophenia_la_SOURCES = conversions.c db.c likelihoods.c linear_algebra.c \
			  regression.c stats.c estimate.c name.c \
asnprintf.c  printf-args.c   vasnprintf.c  asprintf.c   printf-parse.c  vasprintf.c   
pkginclude_HEADERS = conversions.h db.h likelihoods.h linear_algebra.h \
		     regression.h stats.h headers.h estimate.h name.h \
alloca_.h printf-parse.h  vasprintf.h printf-args.h vasnprintf.h    xsize.h


#The programs
bin_PROGRAMS = apop_text_to_db apop_db_to_crosstab apop apop_merge_dbs
apop_text_to_db_SOURCES = apop_text_to_db.c
apop_db_to_crosstab_SOURCES = apop_db_to_crosstab.c
apop_merge_dbs_SOURCES = apop_merge_dbs.c
apop_SOURCES = apop.c

#The documentation
info_TEXINFOS = apophenia.texi
AM_MAKEINFOHTMLFLAGS = --cssinclude = typical.css
