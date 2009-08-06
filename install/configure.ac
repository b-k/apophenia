# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AM_INIT_AUTOMAKE(apophenia, X.XX)
AC_CONFIG_SRCDIR([db.h])
AC_CONFIG_HEADER([config.h])

# The normal /usr/local default just confused too many people
AC_PREFIX_DEFAULT([/usr])

# Windows:
AC_LIBTOOL_WIN32_DLL

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AC_CHECK_LIB([mysqlclient], [mysql_query], , , [-L/usr/lib/mysql -lz])
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

# For SWIG.

# Does the user want to build Python?
AC_MSG_CHECKING([whether to build the Python extension])
AC_ARG_ENABLE([python],
              [AS_HELP_STRING([--enable-python],
                              [set flags to build Python extension (ignored if
--disable-flag-setting is used)])])
test "x$enable_python" = xyes || enable_python=no
AC_MSG_RESULT($enable_python)
AM_CONDITIONAL([USE_PYTHON], [test "x$enable_python" = xyes ])

if test  "x$enable_python" = xyes ; then
AM_PATH_PYTHON
AC_PYTHON_DEVEL
fi

#Check for database systems, then export results to automake
AX_LIB_MYSQL
AX_LIB_SQLITE3
AM_CONDITIONAL([HAVE_MYSQL], [test x$HAVE_MYSQL = xtrue])
AM_CONDITIONAL([HAVE_PYTHON], [test $PYTHON_VERSION > 1 ])


PACKAGE_NAME=apophenia
AX_CREATE_PKGCONFIG_INFO

ACX_PTHREAD

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
