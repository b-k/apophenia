AM_CFLAGS = -g -Wall -O3 $(SQLITE_CFLAGS) $(MYSQL_CFLAGS) $(OPENMP_CFLAGS)
AM_LDFLAGS = $(SQLITE_LDFLAGS) $(OPENMP_CFLAGS) $(MYSQL_LDFLAGS)
CC = $(PTHREAD_CC)
ACLOCAL_AMFLAGS=-I m4

lib_LTLIBRARIES = libapophenia.la
libapophenia_la_SOURCES = \
            apop_arms.c apop_asst.c apop_bootstrap.c apop_conversions.c \
            apop_data.c apop_db.c apop_fexact.c apop_hist.c 	        \
			apop_linear_algebra.c apop_linear_constraint.c apop_mapply.c \
			apop_missing_data.c apop_mle.c apop_model.c   \
			apop_name.c apop_output.c apop_rake.c     \
            apop_regression.c apop_settings.c apop_smoothing.c apop_sort.c \
            apop_stats.c apop_tests.c apop_update.c	 apop_vtables.c            \
			asprintf.c 					\
			transform/apop_dcompose.c   \
			transform/apop_dconstrain.c transform/apop_fix_params.c \
			transform/apop_coordinate_transform.c 	\
			transform/apop_mixture.c 	transform/apop_stack.c    \
			model/apop_bernoulli.c      model/apop_beta.c			\
			model/apop_dirichlet.c      \
			model/apop_exponential.c    model/apop_gamma.c          \
            model/apop_histogram.c      model/apop_loess.c          \
			model/apop_multinomial.c    model/apop_multivariate_normal.c  \
			model/apop_normal.c         model/apop_ols.c 		 	\
            model/apop_pmf.c            model/apop_poisson.c        \
			model/apop_probit.c	        model/apop_t.c				\
			model/apop_uniform.c        model/apop_wishart.c        \
            model/apop_yule.c 	 		model/apop_zipf.c 

include_HEADERS = apop.h

pkgconfig_DATA= apophenia.pc
pkgconfigdir = `pkg-config --variable pc_path pkg-config | cut -f1 -d:`

SUBDIRS=. tests

#The programs
bin_PROGRAMS = apop_text_to_db apop_db_to_crosstab apop_plot_query 
apop_text_to_db_SOURCES = cmd/apop_text_to_db.c
apop_db_to_crosstab_SOURCES = cmd/apop_db_to_crosstab.c
apop_plot_query_SOURCES = cmd/apop_plot_query.c

apop_text_to_db_LDADD = libapophenia.la
apop_db_to_crosstab_LDADD = libapophenia.la
apop_plot_query_LDADD = libapophenia.la

VPATH= @srcdir@

EXTRA_DIST = docs/doxyconfig $(srcdir)/docs/*.png $(srcdir)/docs/*.gif \
		$(srcdir)/docs/*.js 	\
		docs/edit_group $(srcdir)/docs/*.html docs/typical.css \
		docs/make_model_doc.awk docs/documentation.h rpm.spec \
		eg apop_db_mysql.c apop_db_sqlite.c \
		apop_internal.h COPYING2 apophenia.pc

DISTCLEANFILES = apophenia.pc* apophenia-uninstalled.pc apophenia-uninstalled.sh tests/*dat* tests/printing_sample 1

%.c: %.m4.c
	m4 -P prep_variadics.m4 $< > $@

%.h: %.m4.h
	m4 -P prep_variadics.m4 $< > $@

#The Doxygen documentation, which you'll have to call yourself (via make doc).

doc:
	-rm tests/apophenia/*h
	cat model/*.c transform/*.c | awk -f docs/make_model_doc.awk > docs/model_doc.h
	doxygen docs/doxyconfig
	cp docs/flake.gif docs/search.gif docs/right.png docs/down.png html/
	@sed -i -e 's|Outlineheader \([^ ]*\)\(.*\)</p>|<h2><a class="anchor" name="\1"><div class="trigger" onClick="showBranch('\''\1d'\'');swapFolder('\''\1f'\'')"><img src="right.png" border="0" id="\1f" alt="pip">\2</div></a></h2><div class="branch" id="\1d">|' \
-e 's|endofdiv</p>|</div>|' \
-e 's|ALLBUTTON|<span class="trigger" onClick="showAll();"<a>Expand all </a></span> \| <span class="trigger" onClick="closeAll();"<a>Collapse all </a></span>|'  html/index.html html/outline.html
	@sed -i -e '/index_x/s/- x -/ /' -e '/index_a/s/- a -/ /' -e '/Here is a list of all/d' -e '/<span>[a-z]<\/span><\/a><\/li>/d' html/globals.html
	@sed -i -e '/div class="contents"/aSee also the <a class="el" href="group__models.html">models</a> and <a class="el" href="group__settings.html">settings</a> pages.' html/globals.html
	@sed -i 's|<td width="100%"></td>||g' html/*html
	@sed -i -f docs/edit_group html/group__models.html 
	cp docs/*js html/
	sudo cp man/man3/* /usr/share/man/man3/

install-data-local:
	@echo 
	@echo "OK. If you'd like to generate documentation via Doxygen, run make doc; to test, run make check."
	@echo 
