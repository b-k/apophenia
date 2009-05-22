
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
			apop_tests.c apop_smoothing.c					\
			apop_settings.c									\
			vasprintf/asnprintf.c vasprintf/printf-args.c 	\
			vasprintf/vasnprintf.c vasprintf/asprintf.c 	\
			vasprintf/printf-parse.c vasprintf/vasprintf.c  \
			model/apop_multivariate_normal.c				\
			model/apop_multinomial_probit.c					\
			model/apop_normal.c model/apop_OLS.c 		 	\
			model/apop_WLS.c model/apop_exponential.c 		\
			model/apop_gamma.c model/apop_zipf.c 		    \
			model/apop_poisson.c	                    	\
			model/apop_uniform.c model/apop_beta.c			\
			model/apop_waring.c model/apop_yule.c 	 		\
			model/apop_bernoulli.c model/apop_binomial.c 	\
			model/apop_histogram.c model/apop_kernel.c 

pkginclude_HEADERS = 	conversions.h db.h likelihoods.h linear_algebra.h \
				bootstrap.h regression.h stats.h headers.h types.h output.h \
				histogram.h	mapply.h missing_data.h			\
		     	model/model.h asst.h settings.h				\
		     	variadic.h									\
			vasprintf/alloca_.h vasprintf/printf-parse.h  	\
			vasprintf/vasprintf.h 		vasprintf/xsize.h 	\
			vasprintf/printf-args.h vasprintf/vasnprintf.h 	

#SWIG

python_PYTHON = apop.py

apop.py: apop.i
	@if test -e /usr/include/python$(PYTHON_VERSION)/Python.h ; then \
	swig -python -globals avar apop.i ; \
	sed -i "s/:print/:_print/" apop.py ; \
	sed -i "s/def print/def _print/" apop.py ; \
	$(CC) -c -fpic -std=gnu99  apop_wrap.c -I/usr/local/include -I/usr/include/python$(PYTHON_VERSION)/ ; \
	$(CC) -shared apop_wrap.o -o _apop.so -lapophenia -lgsl -lgslcblas -lsqlite3 ; \
	else \
	echo "NOT setting up the Python interface, because I couldn't find /usr/include/python$(PYTHON_VERSION)/Python.h. Maybe install a Python$(PYTHON_VERSION)-dev package?" ; \
	touch apop.py ; \
	fi

install-exec-hook:
	@if test -e _apop.so ; then \
	cp _apop.so ${pythondir} ; \
	echo "Python interface NOT set up, because I couldn't find /usr/include/python$(PYTHON_VERSION)/Python.h. Maybe install a Python$(PYTHON_VERSION)-dev package?" ; \
	fi

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

EXTRA_DIST = doxy/doxyconfig doxy/right.png doxy/down.png doxy/flake.gif doxy/foot.html \
		doxy/head.html doxy/search.gif doxy/typical.css	documentation.h \
		tests/Makefile tests/Readme tests/slower.c tests/nist_tests.c tests/numacc4.dat \
		tests/pontius.dat tests/test.c tests/test_data_mixed tests/test_data tests/test_data2 \
		tests/data-mixed tests/test_data_nans tests/wampler1.dat \
		apop_db_mysql.c apop_db_sqlite.c	apop_findzeros.c	\
		apop.h	apop.i COPYING2

#The documentation

#the texinfo version, which I didn't like in the end. 
#info_TEXINFOS = apophenia.texi
#AM_MAKEINFOHTMLFLAGS = --cssinclude = typical.css

#the doxygen version. You've gotta call it yourself.
#the doxygen version. If you're Ben K., then it will post to Sourceforge.

doc: 
	cp model/*.c .
	doxygen doxy/doxyconfig
	cp doxy/flake.gif doxy/search.gif doxy/right.png doxy/down.png html/
	sed -i -e 's|Outlineheader \([^ ]*\)\(.*\)<p>|<h2><a class="anchor" name="\1"><div class="trigger" onClick="showBranch('\''\1d'\'');swapFolder('\''\1f'\'')"><img src="right.png" border="0" id="\1f" alt="pip">\2</div></a></h2><div class="branch" id="\1d">|' \
-e 's|endofdiv<p>|</div>|' \
-e 's|ALLBUTTON|<span class="trigger" onClick="showAll();"<a>Expand all </a></span> \| <span class="trigger" onClick="closeAll();"<a>Collapse all </a></span>|'  html/outline.html
	sed -i -e '/index_a/s/- a -/ /' -e '15,27d' html/globals.html
	sed -i 's|<td width="100%"></td>||g' html/*html
	cp doxy/*js html/
	rsync -avPe ssh html/ ben_k,apophenia@web.sourceforge.net:htdocs/doc/
	sudo cp man/man3/* /usr/share/man/man3/

test:
	cd tests; rm test.o; make; ./test;

install-data-local:
	test -z "$(includedir)" || $(mkdir_p) "$(includedir)"
	cp apop.h $(includedir)/
