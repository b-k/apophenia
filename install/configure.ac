# Process this file with autoconf to produce a configure script.

m4_define([m4_apop_version], [m4_esyscmd_s(date +%Y%m%d)]) #will switch to this soon.

AC_PREREQ(2.60)
AC_INIT([apophenia], [0.999e], [fluffmail@f-m.fm])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([apop_arms.c])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

# The normal /usr/local default confused too many people
##AC_PREFIX_DEFAULT([/usr])

# libtool:
LT_INIT

# check linker script support
gl_LD_VERSION_SCRIPT

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99

ACX_PTHREAD
AC_OPENMP

# Checks for libraries.
## math library
LT_LIB_M
## GNU Scientific Library (GSL)
AX_PATH_GSL([1.12.0],[],[AC_MSG_ERROR(could not find required version of GSL)])
## DataBase system libraries
#### MySQL library
AX_LIB_MYSQL
#### SQLite3 library
AX_LIB_SQLITE3

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

# Checks for tests tools
AC_PATH_PROGS([BC],[bc],[/usr/bin/bc])
AC_PATH_PROGS([SQLITE3],[sqlite3],[/usr/bin/sqlite3])

# Run only the basic tests unless asked to run the full suite.
AC_MSG_CHECKING([whether to run extended tests])
AC_ARG_ENABLE([extended-tests],
      [AS_HELP_STRING([--enable-extended-tests], [run numeric torture tests (may be slow)])],
                  [enable_extended_tests="yes"], [enable_extended_tests="no"])
AC_MSG_RESULT([$enable_extended_tests])
AM_CONDITIONAL([EXTENDED_TESTS], [test "X$enable_extended_tests" != "Xno"])

AC_CONFIG_FILES([
	apophenia.pc
	Makefile
	transform/Makefile
	model/Makefile
	cmd/Makefile
	tests/Makefile
	docs/doxygen.conf
	docs/Makefile
	eg/Makefile
	])

AC_CONFIG_FILES([
	tests/utilities_test
	],
	[
	chmod a+x tests/utilities_test
	])
AC_OUTPUT

##
## eof
