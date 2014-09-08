#The programs
bin_PROGRAMS = \
	apop_text_to_db \
	apop_db_to_crosstab \
	apop_plot_query

AM_CFLAGS = \
	-I $(top_srcdir) \
	$(GSL_CFLAGS)

LDADD = \
	$(top_builddir)/libapophenia.la \
	$(GSL_LIBS)

