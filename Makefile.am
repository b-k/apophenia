
#The library
#AM_LDFLAGS =  -lgsl -lgslcblas -lsqlite3 -L/usr/lib/mysql -lmysqlclient -lz 
AM_LDFLAGS =  -no-undefined
AM_CFLAGS = -g -Wall -I./ -I/usr/include/mysql -pthread
#apop_text_to_db_LDFLAGS =  -lgsl -lgslcblas -lsqlite3
lib_LTLIBRARIES = libapophenia.la
libapophenia_la_SOURCES = apop_conversions.c apop_db.c 		\
			apop_linear_algebra.c apop_mle.c apop_hist.c 	\
			apop_asst.c apop_data.c apop_regression.c 		\
			apop_fexact.c apop_stats.c apop_bootstrap.c 	\
			apop_name.c apop_output.c apop_model.c 			\
			apop_mapply.c apop_missing_data.c				\
			apop_linear_constraint.c apop_sort.c			\
			apop_model_fix_params.c	 apop_update.c			\
			vasprintf/asnprintf.c vasprintf/printf-args.c 	\
			vasprintf/vasnprintf.c vasprintf/asprintf.c 	\
			vasprintf/printf-parse.c 						\
			vasprintf/vasprintf.c  							\
			model/apop_multivariate_normal.c				\
			model/apop_normal.c model/apop_OLS.c 		 	\
			model/apop_WLS.c	model/apop_logit.c 			\
			model/apop_exponential.c model/apop_exponential_rank.c			\
			model/apop_gamma.c model/apop_gamma_rank.c 		\
			model/apop_poisson.c	model/apop_probit.c 	\
			model/apop_uniform.c model/apop_beta.c			\
			model/apop_waring.c model/apop_waring_rank.c 	\
			model/apop_yule.c 	model/apop_yule_rank.c 		\
			model/apop_zipf.c model/apop_zipf_rank.c 		\
			model/apop_bernoulli.c model/apop_binomial.c 	\
			model/apop_histogram.c model/apop_kernel.c 

pkginclude_HEADERS = 	conversions.h db.h likelihoods.h linear_algebra.h \
				bootstrap.h regression.h stats.h headers.h types.h output.h \
				histogram.h	mapply.h missing_data.h			\
		     	model/model.h asst.h 						\
			vasprintf/alloca_.h vasprintf/printf-parse.h  	\
			vasprintf/vasprintf.h 		vasprintf/xsize.h 	\
			vasprintf/printf-args.h vasprintf/vasnprintf.h 	


#The programs
bin_PROGRAMS = apop_text_to_db apop_db_to_crosstab apop_merge_dbs apop_plot_query
apop_text_to_db_SOURCES = cmd_apop_text_to_db.c
apop_db_to_crosstab_SOURCES = cmd_apop_db_to_crosstab.c
apop_merge_dbs_SOURCES = cmd_apop_merge_dbs.c
apop_plot_query_SOURCES = cmd_apop_plot_query.c

apop_text_to_db_LDADD = libapophenia.la
apop_db_to_crosstab_LDADD = libapophenia.la
apop_merge_dbs_LDADD = libapophenia.la
apop_plot_query_LDADD = libapophenia.la

EXTRA_DIST = doxy/doxyconfig doxy/flake.gif doxy/foot.html \
		doxy/head.html doxy/search.gif doxy/typical.css	documentation.h \
		tests/Makefile tests/Readme tests/nist_tests.c tests/numacc4.dat \
		tests/pontius.dat tests/test.c tests/test_data tests/test_data2 \
		tests/test_data_nans tests/wampler1.dat \
		apop_db_mysql.c apop_db_sqlite.c	apop_findzeros.c	\
		apophenia/*.h apop.h	\
		COPYING2

#The documentation

#the texinfo version, which I didn't like in the end. 
#info_TEXINFOS = apophenia.texi
#AM_MAKEINFOHTMLFLAGS = --cssinclude = typical.css

#the doxygen version. You've gotta call it yourself.

doc: 
	doxygen doxy/doxyconfig
	cp doxy/flake.gif doxy/search.gif html/
	rsync -avPe ssh html/ $(SFSSH):ap/htdocs/doc/
	sudo cp man/man3/* /usr/share/man/man3/

test:
	cd tests; make; ./test;

install-data-local:
	cp apop.h $(includedir)
