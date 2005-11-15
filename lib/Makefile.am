
#AM_CFLAGS = -g -O3 -Wall -I./
AM_CFLAGS = -g -Wall -I./
LIBADD = -lgsl -lgslcblas -lsqlite3
#AUTOMAKE_OPTIONS = no-dependencies

#The library
LDADD = libapophenia.la
AM_LDFLAGS =  -lgsl -lgslcblas -lsqlite3 #-lefence
#apop_text_to_db_LDADD = libapophenia.la
#apop_text_to_db_LDFLAGS =  -lgsl -lgslcblas -lsqlite3
lib_LTLIBRARIES = libapophenia.la
libapophenia_la_SOURCES = conversions.c db.c likelihoods.c linear_algebra.c \
			asst.c distributions.c  regression.c apop_stats.c estimate.c name.c apop_output.c \
asnprintf.c  printf-args.c   vasnprintf.c  asprintf.c   printf-parse.c  vasprintf.c   bootstrap.c\
wnasrt.c  wnbridge.c  wncnjdn2.c  wnconjg.c  wnmatrices.c  wnmbtr.c wnmem.c  wnmemg.c  wnmemn.c  wnvect.c

pkginclude_HEADERS = conversions.h db.h likelihoods.h linear_algebra.h \
		     distributions.h regression.h stats.h headers.h estimate.h name.h output.h \
alloca_.h printf-parse.h  vasprintf.h printf-args.h vasnprintf.h    xsize.h bootstrap.h \
wnaux.h wnlib.h


#The programs
bin_PROGRAMS = apop_text_to_db apop_db_to_crosstab apop apop_merge_dbs
apop_text_to_db_SOURCES = apop_text_to_db.c
apop_db_to_crosstab_SOURCES = apop_db_to_crosstab.c
apop_merge_dbs_SOURCES = apop_merge_dbs.c
apop_SOURCES = apop.c

#The documentation

#the texinfo version, which I didn't like in the end. 
#info_TEXINFOS = apophenia.texi
#AM_MAKEINFOHTMLFLAGS = --cssinclude = typical.css

#the doxygen version. You've gotta call it yourself.
EXTRA_DIST = $(man_MANS)
install_man3_MANS = man/man3/*

doc: 
	doxygen .doxygen
	cp flake.gif html/

test:
	gcc -o test -lapophenia -lgsl -lgslcblas -lsqlite3 test.c -g
	./test
	./test_script
