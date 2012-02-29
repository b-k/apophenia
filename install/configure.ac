# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT([Apophenia], [X.XX], [fluffmail@f-m.fm])
AM_INIT_AUTOMAKE(apophenia, X.XX)
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([db.h])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# The normal /usr/local default confused too many people
AC_PREFIX_DEFAULT([/usr])

# Windows:
LT_INIT

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99

# Checks for libraries.
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AC_CHECK_LIB([mysqlclient], [mysql_query], , , [-L/usr/lib/mysql -lz])
AC_CHECK_LIB([sqlite3], [sqlite3_exec])

AC_CHECK_HEADER([gsl/gsl_matrix.h], , [AC_MSG_ERROR(
    [Couldn't find the GSL header files (I searched for <gsl/gsl_matrix.h> on the include path). \
If you are using a package manager, don't forget to install the libgsl-devel package, as well as \
libgsl itself (and the same for SQLite3).])])

AC_CHECK_HEADER([sqlite3.h], , [AC_MSG_ERROR(
    [Couldn't find the SQLite3 header file, sqlite3.h, on the include path. \
If you are using a package manager, don't forget to install the libsqlite3-devel \
package, as well as libsqlite3 itself.])])

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

echo \
"-----------------------------------------------------------

Thanks for trying ${PACKAGE_NAME} v.${PACKAGE_VERSION}.

 Root for the install: '${prefix}'.
 Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'

 Package features:
   Python: ${enable_python}

 Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all                - build all binaries
     install            - install everything
     test               - test the installation
     doc                - generate documentation via doxygen 

------------------------------------------------------------"

#mySQL claims these vars as its own.
sed '$a\
\#undef PACKAGE     \
\#undef PACKAGE_BUGREPORT   \
\#undef PACKAGE_NAME        \
\#undef PACKAGE_STRING      \
\#undef PACKAGE_TARNAME     \
\#undef PACKAGE_VERSION     \
\#undef VERSION' < config.h > conftmp
mv conftmp config.h
