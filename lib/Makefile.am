lib_LTLIBRARIES = libapophenia.la
libapophenia_la_SOURCES = conversions.c db.c likelihoods.c linear_algebra.c regression.c stats.c estimate.c name.c
pkginclude_HEADERS = conversions.h db.h likelihoods.h linear_algebra.h regression.h stats.h headers.h estimate.h name.h
AM_CFLAGS = -g -Wall -I./
SUBDIRS = gnulib
ACLOCAL_AMFLAGS = -I gnulib/m4
LIBADD = gnulib/libgnu.la -lgsl -lgslcblas -lsqlite3
#AUTOMAKE_OPTIONS = no-dependencies
