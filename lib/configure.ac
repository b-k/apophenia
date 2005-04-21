dnl Process this file with autoconf to produce a configure script.
AC_INIT(db.c)
AM_INIT_AUTOMAKE(apophenia, .05)

dnl Checks for library functions.
AC_CONFIG_FILES(gnulib/Makefile)
dnl AC_CONFIG_SUBDIRS(gnulib)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

dnl Checks for libraries.
dnl Replace `main' with a function in -lgsl:
AC_CHECK_LIB(gsl, main)
dnl Replace `main' with a function in -lgslcblas:
AC_CHECK_LIB(gslcblas, main)
dnl Replace `main' with a function in -lsqlite3:
AC_CHECK_LIB(sqlite3, main)

dnl Checks for header files.
AC_HEADER_STDC

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST


AC_OUTPUT(Makefile)
