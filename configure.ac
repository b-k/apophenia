#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AM_INIT_AUTOMAKE(apophenia, 0.18)
AC_CONFIG_SRCDIR([db.h])
AC_CONFIG_HEADER([config.h])

# Windows:
AC_LIBTOOL_WIN32_DLL

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgsl':
AC_CHECK_LIB([gsl], [gsl_vector_alloc], , , [-lgslcblas])
# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [cblas_dgemm])
# FIXME: Replace `main' with a function in `-lmysqlclient':
AC_CHECK_LIB([mysqlclient], [mysql_query], [libs="$libs -L/usr/lib/mysql -lmysqlclient"], , [-lz])
# FIXME: Replace `main' with a function in `-lsqlite3':
AC_CHECK_LIB([sqlite3], [sqlite3_exec])
# FIXME: Replace `main' with a function in `-lz':
#AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h inttypes.h limits.h stddef.h stdint.h stdlib.h string.h unistd.h wchar.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset pow regcomp sqrt strcasecmp])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
