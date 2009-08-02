# Part one gets the Python development info
# Part two decides whether to use -pthread or -lpthread


# ===========================================================================
#         http://www.nongnu.org/autoconf-archive/ac_python_devel.html
# ===========================================================================
#
# SYNOPSIS
#
#   AC_PYTHON_DEVEL([version])
#
# DESCRIPTION
#
#   Note: Defines as a precious variable "PYTHON_VERSION". Don't override it
#   in your configure.ac.
#
#   This macro checks for Python and tries to get the include path to
#   'Python.h'. It provides the $(PYTHON_CPPFLAGS) and $(PYTHON_LDFLAGS)
#   output variables. It also exports $(PYTHON_EXTRA_LIBS) and
#   $(PYTHON_EXTRA_LDFLAGS) for embedding Python in your code.
#
#   You can search for some particular version of Python by passing a
#   parameter to this macro, for example ">= '2.3.1'", or "== '2.4'". Please
#   note that you *have* to pass also an operator along with the version to
#   match, and pay special attention to the single quotes surrounding the
#   version number. Don't use "PYTHON_VERSION" for this: that environment
#   variable is declared as precious and thus reserved for the end-user.
#
#   This macro should work for all versions of Python >= 2.1.0. As an end
#   user, you can disable the check for the python version by setting the
#   PYTHON_NOVERSIONCHECK environment variable to something else than the
#   empty string.
#
#   If you need to use this macro for an older Python version, please
#   contact the authors. We're always open for feedback.
#
# LICENSE
#
#   Copyright (c) 2009 Sebastian Huber <sebastian-huber@web.de>
#   Copyright (c) 2009 Alan W. Irwin <irwin@beluga.phys.uvic.ca>
#   Copyright (c) 2009 Rafael Laboissiere <rafael@laboissiere.net>
#   Copyright (c) 2009 Andrew Collier <colliera@ukzn.ac.za>
#   Copyright (c) 2009 Matteo Settenvini <matteo@member.fsf.org>
#   Copyright (c) 2009 Horst Knorr <hk_classes@knoda.org>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AC_PYTHON_DEVEL],[
	#
	# Allow the use of a (user set) custom python version
	#
	AC_ARG_VAR([PYTHON_VERSION],[The installed Python
		version to use, for example '2.3'. This string
		will be appended to the Python interpreter
		canonical name.])

	AC_PATH_PROG([PYTHON],[python[$PYTHON_VERSION]])
	if test -z "$PYTHON"; then
	   AC_MSG_ERROR([Cannot find python$PYTHON_VERSION in your system path])
	   PYTHON_VERSION=""
	fi

	#
	# Check for a version of Python >= 2.1.0
	#
	AC_MSG_CHECKING([for a version of Python >= '2.1.0'])
	ac_supports_python_ver=`$PYTHON -c "import sys; \
		ver = sys.version.split ()[[0]]; \
		print (ver >= '2.1.0')"`
	if test "$ac_supports_python_ver" != "True"; then
		if test -z "$PYTHON_NOVERSIONCHECK"; then
			AC_MSG_RESULT([no])
			AC_MSG_FAILURE([
This version of the AC@&t@_PYTHON_DEVEL macro
doesn't work properly with versions of Python before
2.1.0. You may need to re-run configure, setting the
variables PYTHON_CPPFLAGS, PYTHON_LDFLAGS, PYTHON_SITE_PKG,
PYTHON_EXTRA_LIBS and PYTHON_EXTRA_LDFLAGS by hand.
Moreover, to disable this check, set PYTHON_NOVERSIONCHECK
to something else than an empty string.
])
		else
			AC_MSG_RESULT([skip at user request])
		fi
	else
		AC_MSG_RESULT([yes])
	fi

	#
	# if the macro parameter ``version'' is set, honour it
	#
	if test -n "$1"; then
		AC_MSG_CHECKING([for a version of Python $1])
		ac_supports_python_ver=`$PYTHON -c "import sys; \
			ver = sys.version.split ()[[0]]; \
			print (ver $1)"`
		if test "$ac_supports_python_ver" = "True"; then
	   	   AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
			AC_MSG_ERROR([this package requires Python $1.
If you have it installed, but it isn't the default Python
interpreter in your system path, please pass the PYTHON_VERSION
variable to configure. See ``configure --help'' for reference.
])
			PYTHON_VERSION=""
		fi
	fi

	#
	# Check if you have distutils, else fail
	#
	AC_MSG_CHECKING([for the distutils Python package])
	ac_distutils_result=`$PYTHON -c "import distutils" 2>&1`
	if test -z "$ac_distutils_result"; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		AC_MSG_ERROR([cannot import Python module "distutils".
Please check your Python installation. The error was:
$ac_distutils_result])
		PYTHON_VERSION=""
	fi

	#
	# Check for Python include path
	#
	AC_MSG_CHECKING([for Python include path])
	if test -z "$PYTHON_CPPFLAGS"; then
		python_path=`$PYTHON -c "import distutils.sysconfig; \
           		print (distutils.sysconfig.get_python_inc ());"`
		if test -n "${python_path}"; then
		   	python_path="-I$python_path"
		fi
		PYTHON_CPPFLAGS=$python_path
	fi
	AC_MSG_RESULT([$PYTHON_CPPFLAGS])
	AC_SUBST([PYTHON_CPPFLAGS])

	#
	# Check for Python library path
	#
	AC_MSG_CHECKING([for Python library path])
	if test -z "$PYTHON_LDFLAGS"; then
		# (makes two attempts to ensure we've got a version number
		# from the interpreter)
		ac_python_version=`cat<<EOD | $PYTHON -

# join all versioning strings, on some systems
# major/minor numbers could be in different list elements
from distutils.sysconfig import *
ret = ''
for e in get_config_vars ('VERSION'):
	if (e != None):
		ret += e
print (ret)
EOD`

		if test -z "$ac_python_version"; then
			if test -n "$PYTHON_VERSION"; then
				ac_python_version=$PYTHON_VERSION
			else
				ac_python_version=`$PYTHON -c "import sys; \
					print (sys.version[[:3]])"`
			fi
		fi

		# Make the versioning information available to the compiler
		AC_DEFINE_UNQUOTED([HAVE_PYTHON], ["$ac_python_version"],
                                   [If available, contains the Python version number currently in use.])

		# First, the library directory:
		ac_python_libdir=`cat<<EOD | $PYTHON -

# There should be only one
import distutils.sysconfig
for e in distutils.sysconfig.get_config_vars ('LIBDIR'):
	if e != None:
		print (e)
		break
EOD`

		# Before checking for libpythonX.Y, we need to know
		# the extension the OS we're on uses for libraries
		# (we take the first one, if there's more than one fix me!):
		ac_python_soext=`$PYTHON -c \
		  "import distutils.sysconfig; \
		  print (distutils.sysconfig.get_config_vars('SO')[[0]])"`

		# Now, for the library:
		ac_python_soname=`$PYTHON -c \
		  "import distutils.sysconfig; \
		  print (distutils.sysconfig.get_config_vars('LDLIBRARY')[[0]])"`

		# Strip away extension from the end to canonicalize its name:
		ac_python_library=`echo "$ac_python_soname" | sed "s/${ac_python_soext}$//"`

		# This small piece shamelessly adapted from PostgreSQL python macro;
		# credits goes to momjian, I think. I'd like to put the right name
		# in the credits, if someone can point me in the right direction... ?
		#
		if test -n "$ac_python_libdir" -a -n "$ac_python_library" \
			-a x"$ac_python_library" != x"$ac_python_soname"
		then
			# use the official shared library
			ac_python_library=`echo "$ac_python_library" | sed "s/^lib//"`
			PYTHON_LDFLAGS="-L$ac_python_libdir -l$ac_python_library"
		else
			# old way: use libpython from python_configdir
			ac_python_libdir=`$PYTHON -c \
			  "from distutils.sysconfig import get_python_lib as f; \
			  import os; \
			  print (os.path.join(f(plat_specific=1, standard_lib=1), 'config'));"`
			PYTHON_LDFLAGS="-L$ac_python_libdir -lpython$ac_python_version"
		fi

		if test -z "PYTHON_LDFLAGS"; then
			AC_MSG_ERROR([
  Cannot determine location of your Python DSO. Please check it was installed with
  dynamic libraries enabled, or try setting PYTHON_LDFLAGS by hand.
			])
		fi
	fi
	AC_MSG_RESULT([$PYTHON_LDFLAGS])
	AC_SUBST([PYTHON_LDFLAGS])

	#
	# Check for site packages
	#
	AC_MSG_CHECKING([for Python site-packages path])
	if test -z "$PYTHON_SITE_PKG"; then
		PYTHON_SITE_PKG=`$PYTHON -c "import distutils.sysconfig; \
		        print (distutils.sysconfig.get_python_lib(0,0));"`
	fi
	AC_MSG_RESULT([$PYTHON_SITE_PKG])
	AC_SUBST([PYTHON_SITE_PKG])

	#
	# libraries which must be linked in when embedding
	#
	AC_MSG_CHECKING(python extra libraries)
	if test -z "$PYTHON_EXTRA_LIBS"; then
	   PYTHON_EXTRA_LIBS=`$PYTHON -c "import distutils.sysconfig; \
                conf = distutils.sysconfig.get_config_var; \
                print (conf('LOCALMODLIBS') + ' ' + conf('LIBS'))"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LIBS])
	AC_SUBST(PYTHON_EXTRA_LIBS)

	#
	# linking flags needed when embedding
	#
	AC_MSG_CHECKING(python extra linking flags)
	if test -z "$PYTHON_EXTRA_LDFLAGS"; then
		PYTHON_EXTRA_LDFLAGS=`$PYTHON -c "import distutils.sysconfig; \
			conf = distutils.sysconfig.get_config_var; \
			print (conf('LINKFORSHARED'))"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LDFLAGS])
	AC_SUBST(PYTHON_EXTRA_LDFLAGS)

	#
	# final check to see if everything compiles alright
	#
	AC_MSG_CHECKING([consistency of all components of python development environment])
	# save current global flags
	ac_save_LIBS="$LIBS"
	ac_save_CPPFLAGS="$CPPFLAGS"
	LIBS="$ac_save_LIBS $PYTHON_LDFLAGS $PYTHON_EXTRA_LDFLAGS $PYTHON_EXTRA_LIBS"
	CPPFLAGS="$ac_save_CPPFLAGS $PYTHON_CPPFLAGS"
	AC_LANG_PUSH([C])
	AC_LINK_IFELSE([
		AC_LANG_PROGRAM([[#include <Python.h>]],
				[[Py_Initialize();]])
		],[pythonexists=yes],[pythonexists=no])
	AC_LANG_POP([C])
	# turn back to default flags
	CPPFLAGS="$ac_save_CPPFLAGS"
	LIBS="$ac_save_LIBS"

	AC_MSG_RESULT([$pythonexists])

        if test ! "x$pythonexists" = "xyes"; then
	   AC_MSG_FAILURE([
  Could not link test program to Python. Maybe the main Python library has been
  installed in some non-standard library path. If so, pass it to configure,
  via the LDFLAGS environment variable.
  Example: ./configure LDFLAGS="-L/usr/non-standard-path/python/lib"
  ============================================================================
   ERROR!
   You probably have to install the development version of the Python package
   for your distribution.  The exact name of this package varies among them.
  ============================================================================
	   ])
	  PYTHON_VERSION=""
	fi

	#
	# all done!
	#
])
#,  
#AC_DEFUN([AC_PYTHON_DEVEL], )

# ===========================================================================
#           http://www.nongnu.org/autoconf-archive/acx_pthread.html
# ===========================================================================
#
# SYNOPSIS
#
#   ACX_PTHREAD([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro figures out how to build C programs using POSIX threads. It
#   sets the PTHREAD_LIBS output variable to the threads library and linker
#   flags, and the PTHREAD_CFLAGS output variable to any special C compiler
#   flags that are needed. (The user can also force certain compiler
#   flags/libs to be tested by setting these environment variables.)
#
#   Also sets PTHREAD_CC to any special C compiler that is needed for
#   multi-threaded programs (defaults to the value of CC otherwise). (This
#   is necessary on AIX to use the special cc_r compiler alias.)
#
#   NOTE: You are assumed to not only compile your program with these flags,
#   but also link it with them as well. e.g. you should link with
#   $PTHREAD_CC $CFLAGS $PTHREAD_CFLAGS $LDFLAGS ... $PTHREAD_LIBS $LIBS
#
#   If you are only building threads programs, you may wish to use these
#   variables in your default LIBS, CFLAGS, and CC:
#
#          LIBS="$PTHREAD_LIBS $LIBS"
#          CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
#          CC="$PTHREAD_CC"
#
#   In addition, if the PTHREAD_CREATE_JOINABLE thread-attribute constant
#   has a nonstandard name, defines PTHREAD_CREATE_JOINABLE to that name
#   (e.g. PTHREAD_CREATE_UNDETACHED on AIX).
#
#   ACTION-IF-FOUND is a list of shell commands to run if a threads library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_PTHREAD.
#
#   Please let the authors know if this macro fails on any platform, or if
#   you have any other suggestions or comments. This macro was based on work
#   by SGJ on autoconf scripts for FFTW (http://www.fftw.org/) (with help
#   from M. Frigo), as well as ac_pthread and hb_pthread macros posted by
#   Alejandro Forero Cuervo to the autoconf macro repository. We are also
#   grateful for the helpful feedback of numerous users.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([ACX_PTHREAD], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_SAVE
AC_LANG_C
acx_pthread_ok=no

# We used to check for pthread.h first, but this fails if pthread.h
# requires special compiler flags (e.g. on True64 or Sequent).
# It gets checked for in the link test anyway.

# First of all, check if the user has set any of the PTHREAD_LIBS,
# etcetera environment variables, and if threads linking works using
# them:
if test x"$PTHREAD_LIBS$PTHREAD_CFLAGS" != x; then
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        AC_MSG_CHECKING([for pthread_join in LIBS=$PTHREAD_LIBS with CFLAGS=$PTHREAD_CFLAGS])
        AC_TRY_LINK_FUNC(pthread_join, acx_pthread_ok=yes)
        AC_MSG_RESULT($acx_pthread_ok)
        if test x"$acx_pthread_ok" = xno; then
                PTHREAD_LIBS=""
                PTHREAD_CFLAGS=""
        fi
        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"
fi

# We must check for the threads library under a number of different
# names; the ordering is very important because some systems
# (e.g. DEC) have both -lpthread and -lpthreads, where one of the
# libraries is broken (non-POSIX).

# Create a list of thread flags to try.  Items starting with a "-" are
# C compiler flags, and other items are library names, except for "none"
# which indicates that we try without any flags at all, and "pthread-config"
# which is a program returning the flags for the Pth emulation library.

acx_pthread_flags="pthreads none -Kthread -kthread lthread -pthread -pthreads -mthreads pthread --thread-safe -mt pthread-config"

# The ordering *is* (sometimes) important.  Some notes on the
# individual items follow:

# pthreads: AIX (must check this before -lpthread)
# none: in case threads are in libc; should be tried before -Kthread and
#       other compiler flags to prevent continual compiler warnings
# -Kthread: Sequent (threads in libc, but -Kthread needed for pthread.h)
# -kthread: FreeBSD kernel threads (preferred to -pthread since SMP-able)
# lthread: LinuxThreads port on FreeBSD (also preferred to -pthread)
# -pthread: Linux/gcc (kernel threads), BSD/gcc (userland threads)
# -pthreads: Solaris/gcc
# -mthreads: Mingw32/gcc, Lynx/gcc
# -mt: Sun Workshop C (may only link SunOS threads [-lthread], but it
#      doesn't hurt to check since this sometimes defines pthreads too;
#      also defines -D_REENTRANT)
#      ... -mt is also the pthreads flag for HP/aCC
# pthread: Linux, etcetera
# --thread-safe: KAI C++
# pthread-config: use pthread-config program (for GNU Pth library)

case "${host_cpu}-${host_os}" in
        *solaris*)

        # On Solaris (at least, for some versions), libc contains stubbed
        # (non-functional) versions of the pthreads routines, so link-based
        # tests will erroneously succeed.  (We need to link with -pthreads/-mt/
        # -lpthread.)  (The stubs are missing pthread_cleanup_push, or rather
        # a function called by this macro, so we could check for that, but
        # who knows whether they'll stub that too in a future libc.)  So,
        # we'll just look for -pthreads and -lpthread first:

        acx_pthread_flags="-pthreads pthread -mt -pthread $acx_pthread_flags"
        ;;
esac

if test x"$acx_pthread_ok" = xno; then
for flag in $acx_pthread_flags; do

        case $flag in
                none)
                AC_MSG_CHECKING([whether pthreads work without any flags])
                ;;

                -*)
                AC_MSG_CHECKING([whether pthreads work with $flag])
                PTHREAD_CFLAGS="$flag"
                ;;

		pthread-config)
		AC_CHECK_PROG(acx_pthread_config, pthread-config, yes, no)
		if test x"$acx_pthread_config" = xno; then continue; fi
		PTHREAD_CFLAGS="`pthread-config --cflags`"
		PTHREAD_LIBS="`pthread-config --ldflags` `pthread-config --libs`"
		;;

                *)
                AC_MSG_CHECKING([for the pthreads library -l$flag])
                PTHREAD_LIBS="-l$flag"
                ;;
        esac

        save_LIBS="$LIBS"
        save_CFLAGS="$CFLAGS"
        LIBS="$PTHREAD_LIBS $LIBS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Check for various functions.  We must include pthread.h,
        # since some functions may be macros.  (On the Sequent, we
        # need a special flag -Kthread to make this header compile.)
        # We check for pthread_join because it is in -lpthread on IRIX
        # while pthread_create is in libc.  We check for pthread_attr_init
        # due to DEC craziness with -lpthreads.  We check for
        # pthread_cleanup_push because it is one of the few pthread
        # functions on Solaris that doesn't have a non-functional libc stub.
        # We try pthread_create on general principles.
        AC_TRY_LINK([#include <pthread.h>],
                    [pthread_t th; pthread_join(th, 0);
                     pthread_attr_init(0); pthread_cleanup_push(0, 0);
                     pthread_create(0,0,0,0); pthread_cleanup_pop(0); ],
                    [acx_pthread_ok=yes])

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        AC_MSG_RESULT($acx_pthread_ok)
        if test "x$acx_pthread_ok" = xyes; then
                break;
        fi

        PTHREAD_LIBS=""
        PTHREAD_CFLAGS=""
done
fi

# Various other checks:
if test "x$acx_pthread_ok" = xyes; then
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Detect AIX lossage: JOINABLE attribute is called UNDETACHED.
	AC_MSG_CHECKING([for joinable pthread attribute])
	attr_name=unknown
	for attr in PTHREAD_CREATE_JOINABLE PTHREAD_CREATE_UNDETACHED; do
	    AC_TRY_LINK([#include <pthread.h>], [int attr=$attr; return attr;],
                        [attr_name=$attr; break])
	done
        AC_MSG_RESULT($attr_name)
        if test "$attr_name" != PTHREAD_CREATE_JOINABLE; then
            AC_DEFINE_UNQUOTED(PTHREAD_CREATE_JOINABLE, $attr_name,
                               [Define to necessary symbol if this constant
                                uses a non-standard name on your system.])
        fi

        AC_MSG_CHECKING([if more special flags are required for pthreads])
        flag=no
        case "${host_cpu}-${host_os}" in
            *-aix* | *-freebsd* | *-darwin*) flag="-D_THREAD_SAFE";;
            *solaris* | *-osf* | *-hpux*) flag="-D_REENTRANT";;
        esac
        AC_MSG_RESULT(${flag})
        if test "x$flag" != xno; then
            PTHREAD_CFLAGS="$flag $PTHREAD_CFLAGS"
        fi

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        # More AIX lossage: must compile with xlc_r or cc_r
	if test x"$GCC" != xyes; then
          AC_CHECK_PROGS(PTHREAD_CC, xlc_r cc_r, ${CC})
        else
          PTHREAD_CC=$CC
	fi
else
        PTHREAD_CC="$CC"
fi

AC_SUBST(PTHREAD_LIBS)
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_CC)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pthread_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_PTHREAD,1,[Define if you have POSIX threads libraries and header files.]),[$1])
        :
else
        acx_pthread_ok=no
        $2
fi
AC_LANG_RESTORE
])dnl ACX_PTHREAD
# ===========================================================================
#          http://www.nongnu.org/autoconf-archive/ax_lib_mysql.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_MYSQL([MINIMUM-VERSION])
#
# DESCRIPTION
#
#   This macro provides tests of availability of MySQL client library of
#   particular version or newer.
#
#   AX_LIB_MYSQL macro takes only one argument which is optional. If there
#   is no required version passed, then macro does not run version test.
#
#   The --with-mysql option takes one of three possible values:
#
#   no - do not check for MySQL client library
#
#   yes - do check for MySQL library in standard locations (mysql_config
#   should be in the PATH)
#
#   path - complete path to mysql_config utility, use this option if
#   mysql_config can't be found in the PATH
#
#   This macro calls:
#
#     AC_SUBST(MYSQL_CFLAGS)
#     AC_SUBST(MYSQL_LDFLAGS)
#     AC_SUBST(MYSQL_VERSION)
#
#   And sets:
#
#     HAVE_MYSQL
#
# LICENSE
#
#   Copyright (c) 2008 Mateusz Loskot <mateusz@loskot.net>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_LIB_MYSQL],
[
    AC_ARG_WITH([mysql],
        AC_HELP_STRING([--with-mysql=@<:@ARG@:>@],
            [use MySQL client library @<:@default=yes@:>@, optionally specify path to mysql_config]
        ),
        [
        if test "$withval" = "no"; then
            want_mysql="no"
        elif test "$withval" = "yes"; then
            want_mysql="yes"
        else
            want_mysql="yes"
            MYSQL_CONFIG="$withval"
        fi
        ],
        [want_mysql="yes"]
    )

    MYSQL_CFLAGS=""
    MYSQL_LDFLAGS=""
    MYSQL_VERSION=""

    dnl
    dnl Check MySQL libraries (libpq)
    dnl

    if test "$want_mysql" = "yes"; then

        if test -z "$MYSQL_CONFIG" -o test; then
            AC_PATH_PROG([MYSQL_CONFIG], [mysql_config], [no])
        fi

        if test "$MYSQL_CONFIG" != "no"; then
            AC_MSG_CHECKING([for MySQL libraries])

            MYSQL_CFLAGS="`$MYSQL_CONFIG --cflags`"
            MYSQL_LDFLAGS="`$MYSQL_CONFIG --libs`"

            MYSQL_VERSION=`$MYSQL_CONFIG --version`

            AC_DEFINE([HAVE_MYSQL], [1],
                [Define to 1 if MySQL libraries are available])

            found_mysql="yes"
            AC_MSG_RESULT([yes])
        else
            found_mysql="no"
            AC_MSG_RESULT([no])
        fi
    fi

    dnl
    dnl Check if required version of MySQL is available
    dnl


    mysql_version_req=ifelse([$1], [], [], [$1])

    if test "$found_mysql" = "yes" -a -n "$mysql_version_req"; then

        AC_MSG_CHECKING([if MySQL version is >= $mysql_version_req])

        dnl Decompose required version string of MySQL
        dnl and calculate its number representation
        mysql_version_req_major=`expr $mysql_version_req : '\([[0-9]]*\)'`
        mysql_version_req_minor=`expr $mysql_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
        mysql_version_req_micro=`expr $mysql_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        if test "x$mysql_version_req_micro" = "x"; then
            mysql_version_req_micro="0"
        fi

        mysql_version_req_number=`expr $mysql_version_req_major \* 1000000 \
                                   \+ $mysql_version_req_minor \* 1000 \
                                   \+ $mysql_version_req_micro`

        dnl Decompose version string of installed MySQL
        dnl and calculate its number representation
        mysql_version_major=`expr $MYSQL_VERSION : '\([[0-9]]*\)'`
        mysql_version_minor=`expr $MYSQL_VERSION : '[[0-9]]*\.\([[0-9]]*\)'`
        mysql_version_micro=`expr $MYSQL_VERSION : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        if test "x$mysql_version_micro" = "x"; then
            mysql_version_micro="0"
        fi

        mysql_version_number=`expr $mysql_version_major \* 1000000 \
                                   \+ $mysql_version_minor \* 1000 \
                                   \+ $mysql_version_micro`

        mysql_version_check=`expr $mysql_version_number \>\= $mysql_version_req_number`
        if test "$mysql_version_check" = "1"; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no])
        fi
    fi

    AC_SUBST([MYSQL_VERSION])
    AC_SUBST([MYSQL_CFLAGS])
    AC_SUBST([MYSQL_LDFLAGS])
])
# ===========================================================================
#         http://www.nongnu.org/autoconf-archive/ax_lib_sqlite3.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_SQLITE3([MINIMUM-VERSION])
#
# DESCRIPTION
#
#   Test for the SQLite 3 library of a particular version (or newer)
#
#   This macro takes only one optional argument, required version of SQLite
#   3 library. If required version is not passed, 3.0.0 is used in the test
#   of existance of SQLite 3.
#
#   If no intallation prefix to the installed SQLite library is given the
#   macro searches under /usr, /usr/local, and /opt.
#
#   This macro calls:
#
#     AC_SUBST(SQLITE3_CFLAGS)
#     AC_SUBST(SQLITE3_LDFLAGS)
#     AC_SUBST(SQLITE3_VERSION)
#
#   And sets:
#
#     HAVE_SQLITE3
#
# LICENSE
#
#   Copyright (c) 2008 Mateusz Loskot <mateusz@loskot.net>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_LIB_SQLITE3],
[
    AC_ARG_WITH([sqlite3],
        AC_HELP_STRING(
            [--with-sqlite3=@<:@ARG@:>@],
            [use SQLite 3 library @<:@default=yes@:>@, optionally specify the prefix for sqlite3 library]
        ),
        [
        if test "$withval" = "no"; then
            WANT_SQLITE3="no"
        elif test "$withval" = "yes"; then
            WANT_SQLITE3="yes"
            ac_sqlite3_path=""
        else
            WANT_SQLITE3="yes"
            ac_sqlite3_path="$withval"
        fi
        ],
        [WANT_SQLITE3="yes"]
    )

    SQLITE3_CFLAGS=""
    SQLITE3_LDFLAGS=""
    SQLITE3_VERSION=""

    if test "x$WANT_SQLITE3" = "xyes"; then

        ac_sqlite3_header="sqlite3.h"

        sqlite3_version_req=ifelse([$1], [], [3.0.0], [$1])
        sqlite3_version_req_shorten=`expr $sqlite3_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
        sqlite3_version_req_major=`expr $sqlite3_version_req : '\([[0-9]]*\)'`
        sqlite3_version_req_minor=`expr $sqlite3_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
        sqlite3_version_req_micro=`expr $sqlite3_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        if test "x$sqlite3_version_req_micro" = "x" ; then
            sqlite3_version_req_micro="0"
        fi

        sqlite3_version_req_number=`expr $sqlite3_version_req_major \* 1000000 \
                                   \+ $sqlite3_version_req_minor \* 1000 \
                                   \+ $sqlite3_version_req_micro`

        AC_MSG_CHECKING([for SQLite3 library >= $sqlite3_version_req])

        if test "$ac_sqlite3_path" != ""; then
            ac_sqlite3_ldflags="-L$ac_sqlite3_path/lib"
            ac_sqlite3_cppflags="-I$ac_sqlite3_path/include"
        else
            for ac_sqlite3_path_tmp in /usr /usr/local /opt ; do
                if test -f "$ac_sqlite3_path_tmp/include/$ac_sqlite3_header" \
                    && test -r "$ac_sqlite3_path_tmp/include/$ac_sqlite3_header"; then
                    ac_sqlite3_path=$ac_sqlite3_path_tmp
                    ac_sqlite3_cppflags="-I$ac_sqlite3_path_tmp/include"
                    ac_sqlite3_ldflags="-L$ac_sqlite3_path_tmp/lib"
                    break;
                fi
            done
        fi

        ac_sqlite3_ldflags="$ac_sqlite3_ldflags -lsqlite3"

        saved_CPPFLAGS="$CPPFLAGS"
        CPPFLAGS="$CPPFLAGS $ac_sqlite3_cppflags"

        AC_LANG_PUSH(C++)
        AC_COMPILE_IFELSE(
            [
            AC_LANG_PROGRAM([[@%:@include <sqlite3.h>]],
                [[
#if (SQLITE_VERSION_NUMBER >= $sqlite3_version_req_number)
// Everything is okay
#else
#  error SQLite version is too old
#endif
                ]]
            )
            ],
            [
            AC_MSG_RESULT([yes])
            success="yes"
            ],
            [
            AC_MSG_RESULT([not found])
            succees="no"
            ]
        )
        AC_LANG_POP([C++])

        CPPFLAGS="$saved_CPPFLAGS"

        if test "$success" = "yes"; then

            SQLITE3_CFLAGS="$ac_sqlite3_cppflags"
            SQLITE3_LDFLAGS="$ac_sqlite3_ldflags"

            ac_sqlite3_header_path="$ac_sqlite3_path/include/$ac_sqlite3_header"

            dnl Retrieve SQLite release version
            if test "x$ac_sqlite3_header_path" != "x"; then
                ac_sqlite3_version=`cat $ac_sqlite3_header_path \
                    | grep '#define.*SQLITE_VERSION.*\"' | sed -e 's/.* "//' \
                        | sed -e 's/"//'`
                if test $ac_sqlite3_version != ""; then
                    SQLITE3_VERSION=$ac_sqlite3_version
                else
                    AC_MSG_WARN([Can not find SQLITE_VERSION macro in sqlite3.h header to retrieve SQLite version!])
                fi
            fi

            AC_SUBST(SQLITE3_CFLAGS)
            AC_SUBST(SQLITE3_LDFLAGS)
            AC_SUBST(SQLITE3_VERSION)
            AC_DEFINE([HAVE_SQLITE3], [], [Have the SQLITE3 library])
        fi
    fi
])
# ===========================================================================
#    http://www.nongnu.org/autoconf-archive/ax_create_pkgconfig_info.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CREATE_PKGCONFIG_INFO [(outputfile, [requires [,libs [,summary [,cflags [, ldflags]]]]])]
#
# DESCRIPTION
#
#   Defaults:
#
#     $1 = $PACKAGE_NAME.pc
#     $2 = (empty)
#     $3 = $PACKAGE_LIBS $LIBS (as set at that point in configure.ac)
#     $4 = $PACKAGE_SUMMARY (or $1 Library)
#     $5 = $CPPFLAGS $PACKAGE_CFLAGS (as set at the point in configure.ac)
#     $6 = $LDFLAGS $PACKAGE_LDFLAGS (as set at the point in configure.ac)
#
#     PACKAGE_NAME defaults to $PACKAGE if not set.
#     PACKAGE_LIBS defaults to -l$PACKAGE_NAME if not set.
#
#   The resulting file is called $PACKAGE.pc.in / $PACKAGE.pc
#
#   You will find this macro most useful in conjunction with
#   ax_spec_defaults that can read good initializers from the .spec file. In
#   consequencd, most of the generatable installable stuff can be made from
#   information being updated in a single place for the whole project.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2008 Sven Verdoolaege <skimo@kotnet.org>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_CREATE_PKGCONFIG_INFO],[dnl
AS_VAR_PUSHDEF([PKGCONFIG_suffix],[ax_create_pkgconfig_suffix])dnl
AS_VAR_PUSHDEF([PKGCONFIG_libdir],[ax_create_pkgconfig_libdir])dnl
AS_VAR_PUSHDEF([PKGCONFIG_libfile],[ax_create_pkgconfig_libfile])dnl
AS_VAR_PUSHDEF([PKGCONFIG_libname],[ax_create_pkgconfig_libname])dnl
AS_VAR_PUSHDEF([PKGCONFIG_version],[ax_create_pkgconfig_version])dnl
AS_VAR_PUSHDEF([PKGCONFIG_description],[ax_create_pkgconfig_description])dnl
AS_VAR_PUSHDEF([PKGCONFIG_requires],[ax_create_pkgconfig_requires])dnl
AS_VAR_PUSHDEF([PKGCONFIG_pkglibs],[ax_create_pkgconfig_pkglibs])dnl
AS_VAR_PUSHDEF([PKGCONFIG_libs],[ax_create_pkgconfig_libs])dnl
AS_VAR_PUSHDEF([PKGCONFIG_ldflags],[ax_create_pkgconfig_ldflags])dnl
AS_VAR_PUSHDEF([PKGCONFIG_cppflags],[ax_create_pkgconfig_cppflags])dnl
AS_VAR_PUSHDEF([PKGCONFIG_generate],[ax_create_pkgconfig_generate])dnl
AS_VAR_PUSHDEF([PKGCONFIG_src_libdir],[ax_create_pkgconfig_src_libdir])dnl
AS_VAR_PUSHDEF([PKGCONFIG_src_headers],[ax_create_pkgconfig_src_headers])dnl

# we need the expanded forms...
test "x$prefix" = xNONE && prefix=$ac_default_prefix
test "x$exec_prefix" = xNONE && exec_prefix='${prefix}'

AC_MSG_CHECKING(our pkgconfig libname)
test ".$PKGCONFIG_libname" != "." || \
PKGCONFIG_libname="ifelse($1,,${PACKAGE_NAME},`basename $1 .pc`)"
test ".$PKGCONFIG_libname" != "." || \
PKGCONFIG_libname="$PACKAGE"
PKGCONFIG_libname=`eval echo "$PKGCONFIG_libname"`
PKGCONFIG_libname=`eval echo "$PKGCONFIG_libname"`
AC_MSG_RESULT($PKGCONFIG_libname)

AC_MSG_CHECKING(our pkgconfig version)
test ".$PKGCONFIG_version" != "." || \
PKGCONFIG_version="${PACKAGE_VERSION}"
test ".$PKGCONFIG_version" != "." || \
PKGCONFIG_version="$VERSION"
PKGCONFIG_version=`eval echo "$PKGCONFIG_version"`
PKGCONFIG_version=`eval echo "$PKGCONFIG_version"`
AC_MSG_RESULT($PKGCONFIG_version)

AC_MSG_CHECKING(our pkgconfig_libdir)
test ".$pkgconfig_libdir" = "." && \
pkgconfig_libdir='${libdir}/pkgconfig'
PKGCONFIG_libdir=`eval echo "$pkgconfig_libdir"`
PKGCONFIG_libdir=`eval echo "$PKGCONFIG_libdir"`
PKGCONFIG_libdir=`eval echo "$PKGCONFIG_libdir"`
AC_MSG_RESULT($pkgconfig_libdir)
test "$pkgconfig_libdir" != "$PKGCONFIG_libdir" && (
AC_MSG_RESULT(expanded our pkgconfig_libdir... $PKGCONFIG_libdir))
AC_SUBST([pkgconfig_libdir])

AC_MSG_CHECKING(our pkgconfig_libfile)
test ".$pkgconfig_libfile" != "." || \
pkgconfig_libfile="ifelse($1,,$PKGCONFIG_libname.pc,`basename $1`)"
PKGCONFIG_libfile=`eval echo "$pkgconfig_libfile"`
PKGCONFIG_libfile=`eval echo "$PKGCONFIG_libfile"`
AC_MSG_RESULT($pkgconfig_libfile)
test "$pkgconfig_libfile" != "$PKGCONFIG_libfile" && (
AC_MSG_RESULT(expanded our pkgconfig_libfile... $PKGCONFIG_libfile))
AC_SUBST([pkgconfig_libfile])

AC_MSG_CHECKING(our package / suffix)
PKGCONFIG_suffix="$program_suffix"
test ".$PKGCONFIG_suffix" != .NONE || PKGCONFIG_suffix=""
AC_MSG_RESULT(${PACKAGE_NAME} / ${PKGCONFIG_suffix})

AC_MSG_CHECKING(our pkgconfig description)
PKGCONFIG_description="ifelse($4,,$PACKAGE_SUMMARY,$4)"
test ".$PKGCONFIG_description" != "." || \
PKGCONFIG_description="$PKGCONFIG_libname Library"
PKGCONFIG_description=`eval echo "$PKGCONFIG_description"`
PKGCONFIG_description=`eval echo "$PKGCONFIG_description"`
AC_MSG_RESULT($PKGCONFIG_description)

AC_MSG_CHECKING(our pkgconfig requires)
PKGCONFIG_requires="ifelse($2,,$PACKAGE_REQUIRES,$2)"
PKGCONFIG_requires=`eval echo "$PKGCONFIG_requires"`
PKGCONFIG_requires=`eval echo "$PKGCONFIG_requires"`
AC_MSG_RESULT($PKGCONFIG_requires)

AC_MSG_CHECKING(our pkgconfig ext libs)
PKGCONFIG_pkglibs="$PACKAGE_LIBS"
test ".$PKGCONFIG_pkglibs" != "." || PKGCONFIG_pkglibs="-l$PKGCONFIG_libname"
PKGCONFIG_libs="ifelse($3,,$PKGCONFIG_pkglibs $LIBS,$3)"
PKGCONFIG_libs=`eval echo "$PKGCONFIG_libs"`
PKGCONFIG_libs=`eval echo "$PKGCONFIG_libs"`
AC_MSG_RESULT($PKGCONFIG_libs)

AC_MSG_CHECKING(our pkgconfig cppflags)
PKGCONFIG_cppflags="ifelse($5,,$CPPFLAGS $PACKAGE_CFLAGS,$5)"
PKGCONFIG_cppflags=`eval echo "$PKGCONFIG_cppflags"`
PKGCONFIG_cppflags=`eval echo "$PKGCONFIG_cppflags"`
AC_MSG_RESULT($PKGCONFIG_cppflags)

AC_MSG_CHECKING(our pkgconfig ldflags)
PKGCONFIG_ldflags="ifelse($6,,$LDFLAGS $PACKAGE_LDFLAGS,$5)"
PKGCONFIG_ldflags=`eval echo "$PKGCONFIG_ldflags"`
PKGCONFIG_ldflags=`eval echo "$PKGCONFIG_ldflags"`
AC_MSG_RESULT($PKGCONFIG_ldflags)

test ".$PKGCONFIG_generate" != "." || \
PKGCONFIG_generate="ifelse($1,,$PKGCONFIG_libname.pc,$1)"
PKGCONFIG_generate=`eval echo "$PKGCONFIG_generate"`
PKGCONFIG_generate=`eval echo "$PKGCONFIG_generate"`
test "$pkgconfig_libfile" != "$PKGCONFIG_generate" && (
AC_MSG_RESULT(generate the pkgconfig later... $PKGCONFIG_generate))

if test ".$PKGCONFIG_src_libdir" = "." ; then
PKGCONFIG_src_libdir=`pwd`
PKGCONFIG_src_libdir=`AS_DIRNAME("$PKGCONFIG_src_libdir/$PKGCONFIG_generate")`
test ! -d $PKGCONFIG_src_libdir/src || \
PKGCONFIG_src_libdir="$PKGCONFIG_src_libdir/src"
case ".$objdir" in
*libs) PKGCONFIG_src_libdir="$PKGCONFIG_src_libdir/$objdir" ;; esac
AC_MSG_RESULT(noninstalled pkgconfig -L $PKGCONFIG_src_libdir)
fi

if test ".$PKGCONFIG_src_headers" = "." ; then
PKGCONFIG_src_headers=`pwd`
v="$ac_top_srcdir" ;
test ".$v" != "." || v="$ax_spec_dir"
test ".$v" != "." || v="$srcdir"
case "$v" in /*) PKGCONFIG_src_headers="" ;; esac
PKGCONFIG_src_headers=`AS_DIRNAME("$PKGCONFIG_src_headers/$v/x")`
test ! -d $PKGCONFIG_src_headers/incl[]ude || \
PKGCONFIG_src_headers="$PKGCONFIG_src_headers/incl[]ude"
AC_MSG_RESULT(noninstalled pkgconfig -I $PKGCONFIG_src_headers)
fi


dnl AC_CONFIG_COMMANDS crap disallows to use $PKGCONFIG_libfile here...
AC_CONFIG_COMMANDS([$ax_create_pkgconfig_generate],[
pkgconfig_generate="$ax_create_pkgconfig_generate"
if test ! -f "$pkgconfig_generate.in"
then generate="true"
elif grep ' generated by configure ' $pkgconfig_generate.in >/dev/null
then generate="true"
else generate="false";
fi
if $generate ; then
AC_MSG_NOTICE(creating $pkgconfig_generate.in)
cat > $pkgconfig_generate.in <<AXEOF
# generated by configure / remove this line to disable regeneration
prefix=@prefix@
exec_prefix=@exec_prefix@
bindir=@bindir@
libdir=@libdir@
datarootdir=@datarootdir@
datadir=@datadir@
sysconfdir=@sysconfdir@
includedir=@includedir@
package=@PACKAGE@
suffix=@suffix@

Name: @PACKAGE_NAME@
Description: @PACKAGE_DESCRIPTION@
Version: @PACKAGE_VERSION@
Requires: @PACKAGE_REQUIRES@
Libs: -L\${libdir} @LDFLAGS@ @LIBS@
Cflags: -I\${includedir} @CPPFLAGS@
AXEOF
fi # DONE generate $pkgconfig_generate.in
AC_MSG_NOTICE(creating $pkgconfig_generate)
cat >conftest.sed <<AXEOF
s|@prefix@|${pkgconfig_prefix}|
s|@exec_prefix@|${pkgconfig_execprefix}|
s|@bindir@|${pkgconfig_bindir}|
s|@libdir@|${pkgconfig_libdir}|
s|@datarootdir@|${pkgconfig_datarootdir}|
s|@datadir@|${pkgconfig_datadir}|
s|@sysconfdir@|${pkgconfig_sysconfdir}|
s|@includedir@|${pkgconfig_includedir}|
s|@suffix@|${pkgconfig_suffix}|
s|@PACKAGE@|${pkgconfig_package}|
s|@PACKAGE_NAME@|${pkgconfig_libname}|
s|@PACKAGE_DESCRIPTION@|${pkgconfig_description}|
s|@PACKAGE_VERSION@|${pkgconfig_version}|
s|@PACKAGE_REQUIRES@|${pkgconfig_requires}|
s|@LIBS@|${pkgconfig_libs}|
s|@LDFLAGS@|${pkgconfig_ldflags}|
s|@CPPFLAGS@|${pkgconfig_cppflags}|
AXEOF
sed -f conftest.sed  $pkgconfig_generate.in > $pkgconfig_generate
if test ! -s $pkgconfig_generate ; then
    AC_MSG_ERROR([$pkgconfig_generate is empty])
fi ; rm conftest.sed # DONE generate $pkgconfig_generate
pkgconfig_uninstalled=`echo $pkgconfig_generate |sed 's/.pc$/-uninstalled.pc/'`
AC_MSG_NOTICE(creating $pkgconfig_uninstalled)
cat >conftest.sed <<AXEOF
s|@prefix@|${pkgconfig_prefix}|
s|@exec_prefix@|${pkgconfig_execprefix}|
s|@bindir@|${pkgconfig_bindir}|
s|@libdir@|${pkgconfig_src_libdir}|
s|@datarootdir@|${pkgconfig_datarootdir}|
s|@datadir@|${pkgconfig_datadir}|
s|@sysconfdir@|${pkgconfig_sysconfdir}|
s|@includedir@|${pkgconfig_src_headers}|
s|@suffix@|${pkgconfig_suffix}|
s|@PACKAGE@|${pkgconfig_package}|
s|@PACKAGE_NAME@|${pkgconfig_libname}|
s|@PACKAGE_DESCRIPTION@|${pkgconfig_description}|
s|@PACKAGE_VERSION@|${pkgconfig_version}|
s|@PACKAGE_REQUIRES@|${pkgconfig_requires}|
s|@LIBS@|${pkgconfig_libs}|
s|@LDFLAGS@|${pkgconfig_ldflags}|
s|@CPPFLAGS@|${pkgconfig_cppflags}|
AXEOF
sed -f conftest.sed $pkgconfig_generate.in > $pkgconfig_uninstalled
if test ! -s $pkgconfig_uninstalled ; then
    AC_MSG_ERROR([$pkgconfig_uninstalled is empty])
fi ; rm conftest.sed # DONE generate $pkgconfig_uninstalled
           pkgconfig_requires_add=`echo ${pkgconfig_requires}`
if test ".$pkgconfig_requires_add" != "." ; then
           pkgconfig_requires_add="pkg-config $pkgconfig_requires_add"
    else   pkgconfig_requires_add=":" ; fi
pkgconfig_uninstalled=`echo $pkgconfig_generate |sed 's/.pc$/-uninstalled.sh/'`
AC_MSG_NOTICE(creating $pkgconfig_uninstalled)
cat >conftest.sed <<AXEOF
s|@prefix@|\"${pkgconfig_prefix}\"|
s|@exec_prefix@|\"${pkgconfig_execprefix}\"|
s|@bindir@|\"${pkgconfig_bindir}\"|
s|@libdir@|\"${pkgconfig_src_libdir}\"|
s|@datarootdir@|\"${pkgconfig_datarootdir}\"|
s|@datadir@|\"${pkgconfig_datadir}\"|
s|@sysconfdir@|\"${pkgconfig_sysconfdir}\"|
s|@includedir@|\"${pkgconfig_src_headers}\"|
s|@suffix@|\"${pkgconfig_suffix}\"|
s|@PACKAGE@|\"${pkgconfig_package}\"|
s|@PACKAGE_NAME@|\"${pkgconfig_libname}\"|
s|@PACKAGE_DESCRIPTION@|\"${pkgconfig_description}\"|
s|@PACKAGE_VERSION@|\"${pkgconfig_version}\"|
s|@PACKAGE_REQUIRES@|\"${pkgconfig_requires}\"|
s|@LIBS@|\"${pkgconfig_libs}\"|
s|@LDFLAGS@|\"${pkgconfig_ldflags}\"|
s|@CPPFLAGS@|\"${pkgconfig_cppflags}\"|
s>Name:>for option\\; do case \"\$option\" in --list-all|--name) echo >
s>Description: *>\\;\\; --help) pkg-config --help \\; echo Buildscript Of >
s>Version: *>\\;\\; --modversion|--version) echo >
s>Requires:>\\;\\; --requires) echo $pkgconfig_requires_add>
s>Libs: *>\\;\\; --libs) echo >
s>Cflags: *>\\;\\; --cflags) echo >
/--libs)/a\\
       $pkgconfig_requires_add
/--cflags)/a\\
       $pkgconfig_requires_add\\
;; --variable=*) eval echo '\$'\`echo \$option | sed -e 's/.*=//'\`\\
;; --uninstalled) exit 0 \\
;; *) ;; esac done
AXEOF
sed -f conftest.sed  $pkgconfig_generate.in > $pkgconfig_uninstalled
if test ! -s $pkgconfig_uninstalled ; then
    AC_MSG_ERROR([$pkgconfig_uninstalled is empty])
fi ; rm conftest.sed # DONE generate $pkgconfig_uninstalled
],[
dnl AC_CONFIG_COMMANDS crap, the AS_PUSHVAR defines are invalid here...
ax_create_pkgconfig_generate="$ax_create_pkgconfig_generate"
pkgconfig_prefix='$prefix'
pkgconfig_execprefix='$exec_prefix'
pkgconfig_bindir='$bindir'
pkgconfig_libdir='$libdir'
pkgconfig_includedir='$includedir'
pkgconfig_datarootdir='$datarootdir'
pkgconfig_datadir='$datadir'
pkgconfig_sysconfdir='$sysconfdir'
pkgconfig_suffix='$ax_create_pkgconfig_suffix'
pkgconfig_package='$PACKAGE_NAME'
pkgconfig_libname='$ax_create_pkgconfig_libname'
pkgconfig_description='$ax_create_pkgconfig_description'
pkgconfig_version='$ax_create_pkgconfig_version'
pkgconfig_requires='$ax_create_pkgconfig_requires'
pkgconfig_libs='$ax_create_pkgconfig_libs'
pkgconfig_ldflags='$ax_create_pkgconfig_ldflags'
pkgconfig_cppflags='$ax_create_pkgconfig_cppflags'
pkgconfig_src_libdir='$ax_create_pkgconfig_src_libdir'
pkgconfig_src_headers='$ax_create_pkgconfig_src_headers'
])dnl
AS_VAR_POPDEF([PKGCONFIG_suffix])dnl
AS_VAR_POPDEF([PKGCONFIG_libdir])dnl
AS_VAR_POPDEF([PKGCONFIG_libfile])dnl
AS_VAR_POPDEF([PKGCONFIG_libname])dnl
AS_VAR_POPDEF([PKGCONFIG_version])dnl
AS_VAR_POPDEF([PKGCONFIG_description])dnl
AS_VAR_POPDEF([PKGCONFIG_requires])dnl
AS_VAR_POPDEF([PKGCONFIG_pkglibs])dnl
AS_VAR_POPDEF([PKGCONFIG_libs])dnl
AS_VAR_POPDEF([PKGCONFIG_ldflags])dnl
AS_VAR_POPDEF([PKGCONFIG_cppflags])dnl
AS_VAR_POPDEF([PKGCONFIG_generate])dnl
AS_VAR_POPDEF([PKGCONFIG_src_libdir])dnl
AS_VAR_POPDEF([PKGCONFIG_src_headers])dnl
])
