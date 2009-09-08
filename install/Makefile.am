
AM_LDFLAGS =  -no-undefined $(PTHREAD_CFLAGS) $(PTHREAD_LIBS) $(MYSQL_LDFLAGS) $(SQLITE_LDFLAGS)
AM_CFLAGS = -g -Wall -O3 -I./ $(PTHREAD_CFLAGS) $(PTHREAD_LIBS) $(MYSQL_CFLAGS) $(SQLITE_CFLAGS) 

CC = $(PTHREAD_CC)
lib_LTLIBRARIES = libapophenia.la
libapophenia_la_SOURCES = \
            apop_conversions.c apop_db.c 	\
			apop_linear_algebra.c apop_mle.c apop_hist.c 	\
			apop_asst.c apop_data.c apop_regression.c 		\
			apop_fexact.c apop_stats.c apop_bootstrap.c 	\
			apop_name.c apop_output.c apop_model.c 			\
			apop_mapply.c apop_missing_data.c				\
			apop_linear_constraint.c apop_arms.c 			\
			apop_model_fix_params.c	 apop_update.c			\
			apop_tests.c apop_smoothing.c apop_settings.c   \
			asprintf.c \
			model/apop_multinomial_probit.c	model/apop_pmf.c \
			model/apop_multivariate_normal.c 				\
			model/apop_normal.c model/apop_OLS.c 		 	\
			model/apop_exponential.c model/apop_histogram.c \
			model/apop_gamma.c model/apop_zipf.c 		    \
			model/apop_poisson.c model/apop_t_f_chi.c		\
			model/apop_uniform.c model/apop_beta.c			\
			model/apop_waring.c model/apop_yule.c 	 		\
			model/apop_bernoulli.c model/apop_binomial.c 	

pkginclude_HEADERS = 	conversions.h db.h likelihoods.h 	\
				arms.h mapply.h stats.h types.h		 		\
				output.h model.h asst.h settings.h  		\
				variadic.h	linear_algebra.h 

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA= apophenia.pc

#SWIG
if HAVE_PYTHON
python_PYTHON = apop.py

apop.py: apop.i
	swig -python -globals avar apop.i ; \
	sed  -e "s/:print/:_print/" -e "s/def print/def _print/" < apop.py > apop.tmp.py; \
	mv apop.tmp.py apop.py; \
	$(CC) -c -fpic -std=gnu99  apop_wrap.c -I$(pkgincludedir) $(PYTHON_CPPFLAGS) ; \
	$(CC) -shared apop_wrap.o $(PYTHON_CPPFLAGS) -o _apop.so -lapophenia -lgsl -lgslcblas $(MYSQL_LDFLAGS) $(SQLITE_LDFLAGS) ;

install-exec-hook: apop.py
	@if test -e _apop.so ; then \
		cp _apop.so ${pythondir} ; \
	fi
else 
apop.py: apop.i
	@echo "NOT setting up the Python interface. Maybe install a Python$(PYTHON_VERSION)-dev package, or run ./configure --enable-python ?" ; \
	touch apop.py ; 
endif

#The programs
bin_PROGRAMS = apop_text_to_db apop_db_to_crosstab apop_merge_dbs apop_plot_query apop_lookup
apop_text_to_db_SOURCES = cmd_apop_text_to_db.c
apop_db_to_crosstab_SOURCES = cmd_apop_db_to_crosstab.c
apop_merge_dbs_SOURCES = cmd_apop_merge_dbs.c
apop_plot_query_SOURCES = cmd_apop_plot_query.c
apop_lookup_SOURCES = cmd_apop_lookup.c

apop_text_to_db_LDADD = libapophenia.la
apop_db_to_crosstab_LDADD = libapophenia.la
apop_merge_dbs_LDADD = libapophenia.la
apop_plot_query_LDADD = libapophenia.la
apop_lookup_LDADD = libapophenia.la

EXTRA_DIST = apophenia.pc doxy/doxyconfig doxy/right.png doxy/down.png doxy/flake.gif doxy/foot.html \
		doxy/head.html doxy/search.gif doxy/typical.css	documentation.h \
		tests/Makefile tests/test.py tests/Readme tests/nist_tests.c tests/numacc4.dat \
		tests/pontius.dat tests/test.c tests/test_data_mixed tests/test_data tests/test_data2 \
		tests/data-mixed tests/test_data_nans tests/wampler1.dat eg/* \
		apop_db_mysql.c apop_db_sqlite.c apop.h apop.i internal.h COPYING2

#The Doxygen documentation, which you'll have to call yourself (via make doc).

doc: 
	cp model/*.c .
	doxygen doxy/doxyconfig
	cp doxy/flake.gif doxy/search.gif doxy/right.png doxy/down.png html/
	sed -i -e 's|Outlineheader \([^ ]*\)\(.*\)<p>|<h2><a class="anchor" name="\1"><div class="trigger" onClick="showBranch('\''\1d'\'');swapFolder('\''\1f'\'')"><img src="right.png" border="0" id="\1f" alt="pip">\2</div></a></h2><div class="branch" id="\1d">|' \
-e 's|endofdiv<p>|</div>|' \
-e 's|ALLBUTTON|<span class="trigger" onClick="showAll();"<a>Expand all </a></span> \| <span class="trigger" onClick="closeAll();"<a>Collapse all </a></span>|'  html/index.html html/outline.html
	sed -i -e '/index_a/s/- a -/ /' -e '15,27d' -e '/Here is a list of all/d' html/globals.html
	sed -i 's|<td width="100%"></td>||g' html/*html
	cp doxy/*js html/
	sudo cp man/man3/* /usr/share/man/man3/

test:
	cd tests; rm test.o; make; ./test;

install-data-local:
	test -z "$(includedir)" || $(mkdir_p) "$(includedir)"
	cp apop.h $(includedir)/
	@echo 
	@echo "OK. If you'd like to generate documentation, run make doc; to test, run make test."
	@echo 
