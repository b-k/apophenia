# Process this file with autoconf to produce a configure script.

m4_define([m4_apop_version], [m4_esyscmd_s(date +%Y%m%d)]) #will switch to this soon.

AC_PREREQ(2.60)
AC_INIT([Apophenia], [0.999], [fluffmail@f-m.fm])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([apop_arms.c])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

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
#AC_CHECK_LIB([mysqlclient], [mysql_query], , , [-L/usr/lib/mysql -lz])
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
AX_C___ATTRIBUTE__

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset pow regcomp sqrt strcasecmp asprintf])

#Check for database systems, then export results to automake
AX_LIB_MYSQL
AX_LIB_SQLITE3
AM_CONDITIONAL([HAVE_MYSQL], [test x$HAVE_MYSQL = xtrue])


PACKAGE_NAME=apophenia

ACX_PTHREAD
AC_OPENMP

AC_CONFIG_FILES([Makefile tests/Makefile apophenia.pc])
AC_OUTPUT

echo \
"-----------------------------------------------------------

Thanks for trying ${PACKAGE_NAME} v.${PACKAGE_VERSION}.

Root for the install: '${prefix}'.
 Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'
"
echo $mysql_message

echo \
"
Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all      - build all binaries
     check    - test the installation
     install  - install everything [e.g., sudo make install]
     doc      - generate documentation via doxygen 

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
