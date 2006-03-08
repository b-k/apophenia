
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
libapophenia_la_SOURCES = apop_conversions.c apop_db.c apop_mle.c apop_linear_algebra.c 		\
			apop_asst.c apop_data.c apop_regression.c apop_stats.c apop_bootstrap.c 			\
			apop_name.c apop_output.c apop_estimate.c 											\
			vasprintf/asnprintf.c vasprintf/printf-args.c 										\
			vasprintf/vasnprintf.c vasprintf/asprintf.c vasprintf/printf-parse.c 				\
			vasprintf/vasprintf.c  																\
			model/apop_gamma.c model/apop_gamma_rank.c model/apop_normal.c\
			model/apop_probit.c model/apop_exponential.c \
			model/apop_poisson.c	model/apop_exponential_rank.c	\
			model/apop_waring.c model/apop_waring_rank.c model/apop_yule.c model/apop_zipf_rank.c \
			model/apop_zipf.c model/apop_OLS.c 													\
			wn/wnasrt.c  wn/wnbridge.c  wn/wncnjdn2.c  wn/wnconjg.c  wn/wnmatrices.c			\
			wn/wnmbtr.c wn/wnmem.c  wn/wnmemg.c  wn/wnmemn.c  wn/wnvect.c

pkginclude_HEADERS = 	conversions.h db.h likelihoods.h linear_algebra.h bootstrap.h 	\
		     	regression.h stats.h headers.h types.h output.h 						\
		     	model/model.h 															\
			vasprintf/alloca_.h vasprintf/printf-parse.h  vasprintf/vasprintf.h 		\
			vasprintf/printf-args.h vasprintf/vasnprintf.h vasprintf/xsize.h 			\
			wn/wnaux.h wn/wnlib.h


#The programs
bin_PROGRAMS = apop_text_to_db apop_db_to_crosstab apop apop_merge_dbs
apop_text_to_db_SOURCES = cmd_apop_text_to_db.c
apop_db_to_crosstab_SOURCES = cmd_apop_db_to_crosstab.c
apop_merge_dbs_SOURCES = cmd_apop_merge_dbs.c
apop_SOURCES = cmd_apop.c

#The documentation

#the texinfo version, which I didn't like in the end. 
#info_TEXINFOS = apophenia.texi
#AM_MAKEINFOHTMLFLAGS = --cssinclude = typical.css

#the doxygen version. You've gotta call it yourself.
EXTRA_DIST = $(man_MANS)
install_man3_MANS = man/man3/*

doc: 
	doxygen doxy/doxyconfig
	cp doxy/flake.gif doxy/search.gif html/
	rsync -avPe ssh html/ $(SFSSH):ap/htdocs/doc/

test:
	gcc -o test -lapophenia -lgsl -lgslcblas -lsqlite3 test.c -g
	./test
	./test_script
