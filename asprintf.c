    /** \cond doxy_ignore */
/* Formatted output to strings.
   Copyright (C) 1999, 2002 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */

/* Mashed into one file by BK */

/* Tell glibc's <stdio.h> to provide a prototype for snprintf().
   This must come before <config.h> because <config.h> may include
   <features.h>, and once <features.h> has been included, it's too late.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE    1
#endif

/* vsprintf with automatic memory allocation. */

#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#if HAVE_ASPRINTF
#else

#ifndef __attribute__
/* This feature is available in gcc versions 2.5 and later.  */
# if __GNUC__ < 2 || (__GNUC__ == 2 && __GNUC_MINOR__ < 5) || __STRICT_ANSI__
#  define __attribute__(Spec) /* empty */
# endif
/* The __-protected variants of `format' and `printf' attributes
   are accepted by gcc versions 2.6.4 (effectively 2.7) and later.  */
# if __GNUC__ < 2 || (__GNUC__ == 2 && __GNUC_MINOR__ < 7)
#  define __format__ format
#  define __printf__ printf
# endif
#endif

#ifdef	__cplusplus
extern "C" {
#endif

/* Write formatted output to a string dynamically allocated with malloc().
   If the memory allocation succeeds, store the address of the string in
   *RESULT and return the number of resulting bytes, excluding the trailing
   NUL.  Upon memory allocation error, or some other error, return -1.  */
extern int asprintf (char **result, const char *format, ...)
       __attribute__ ((__format__ (__printf__, 2, 3)));
extern int vasprintf (char **result, const char *format, va_list args)
       __attribute__ ((__format__ (__printf__, 2, 0)));

//from vasnprintf
extern char * asnprintf (char *resultbuf, size_t *lengthp, const char *format, ...)
       __attribute__ ((__format__ (__printf__, 3, 4)));
extern char * vasnprintf (char *resultbuf, size_t *lengthp, const char *format, va_list args)
       __attribute__ ((__format__ (__printf__, 3, 0)));

#ifdef	__cplusplus
}
#endif


/* Decomposed printf argument list. */

/* Get wint_t.  */
#ifdef HAVE_WINT_T
# include <wchar.h>
#endif

/* Argument types */
typedef enum {
  TYPE_NONE,
  TYPE_SCHAR,
  TYPE_UCHAR,
  TYPE_SHORT,
  TYPE_USHORT,
  TYPE_INT,
  TYPE_UINT,
  TYPE_LONGINT,
  TYPE_ULONGINT,
#ifdef HAVE_LONG_LONG
  TYPE_LONGLONGINT,
  TYPE_ULONGLONGINT,
#endif
  TYPE_DOUBLE,
#ifdef HAVE_LONG_DOUBLE
  TYPE_LONGDOUBLE,
#endif
  TYPE_CHAR,
#ifdef HAVE_WINT_T
  TYPE_WIDE_CHAR,
#endif
  TYPE_STRING,
#ifdef HAVE_WCHAR_T
  TYPE_WIDE_STRING,
#endif
  TYPE_POINTER,
  TYPE_COUNT_SCHAR_POINTER,
  TYPE_COUNT_SHORT_POINTER,
  TYPE_COUNT_INT_POINTER,
  TYPE_COUNT_LONGINT_POINTER
#ifdef HAVE_LONG_LONG
, TYPE_COUNT_LONGLONGINT_POINTER
#endif
} arg_type;

/* Polymorphic argument */
typedef struct {
  arg_type type;
  union
  {
    signed char			a_schar;
    unsigned char		a_uchar;
    short			a_short;
    unsigned short		a_ushort;
    int				a_int;
    unsigned int		a_uint;
    long int			a_longint;
    unsigned long int		a_ulongint;
#ifdef HAVE_LONG_LONG
    long long int		a_longlongint;
    unsigned long long int	a_ulonglongint;
#endif
    float			a_float;
    double			a_double;
#ifdef HAVE_LONG_DOUBLE
    long double			a_longdouble;
#endif
    int				a_char;
#ifdef HAVE_WINT_T
    wint_t			a_wide_char;
#endif
    const char*			a_string;
#ifdef HAVE_WCHAR_T
    const wchar_t*		a_wide_string;
#endif
    void*			a_pointer;
    signed char *		a_count_schar_pointer;
    short *			a_count_short_pointer;
    int *			a_count_int_pointer;
    long int *			a_count_longint_pointer;
#ifdef HAVE_LONG_LONG
    long long int *		a_count_longlongint_pointer;
#endif
  }
  a;
}
argument;

typedef struct {
  size_t count;
  argument *arg;
}
arguments;

/* Fetch the arguments, putting them into a. */
#ifdef STATIC
STATIC
#else
extern
#endif
int printf_fetchargs (va_list args, arguments *a);

/* Parse printf format string. */

/* Flags */
#define FLAG_GROUP	 1	/* ' flag */
#define FLAG_LEFT	 2	/* - flag */
#define FLAG_SHOWSIGN	 4	/* + flag */
#define FLAG_SPACE	 8	/* space flag */
#define FLAG_ALT	16	/* # flag */
#define FLAG_ZERO	32

/* arg_index value indicating that no argument is consumed.  */
#define ARG_NONE	(~(size_t)0)

/* A parsed directive.  */
typedef struct {
  const char* dir_start;
  const char* dir_end;
  int flags;
  const char* width_start;
  const char* width_end;
  size_t width_arg_index;
  const char* precision_start;
  const char* precision_end;
  size_t precision_arg_index;
  char conversion; /* d i o u x X f e E g G c s p n U % but not C S */
  size_t arg_index;
}
char_directive;

/* A parsed format string.  */
typedef struct {
  size_t count;
  char_directive *dir;
  size_t max_width_length;
  size_t max_precision_length;
}
char_directives;


/* Parses the format string.  Fills in the number N of directives, and fills
   in directives[0], ..., directives[N-1], and sets directives[N].dir_start
   to the end of the format string.  Also fills in the arg_type fields of the
   arguments and the needed count of arguments.  */
#ifdef STATIC
STATIC
#else
extern
#endif
int printf_parse (const char *format, char_directives *d, arguments *a);

/*end headers */

char * asnprintf (char *resultbuf, size_t *lengthp, const char *format, ...) {
  va_list args;
  char *result;

  va_start (args, format);
  result = vasnprintf (resultbuf, lengthp, format, args);
  va_end (args);
  return result;
}

int asprintf (char **resultp, const char *format, ...) {
  va_list args;
  int result;

  va_start (args, format);
  result = vasprintf (resultp, format, args);
  va_end (args);
  return result;
}

#ifdef STATIC
STATIC
#endif
int printf_fetchargs (va_list args, arguments *a) {
  size_t i;
  argument *ap;

  for (i = 0, ap = &a->arg[0]; i < a->count; i++, ap++)
    switch (ap->type)
      {
      case TYPE_SCHAR:
	ap->a.a_schar = va_arg (args, /*signed char*/ int);
	break;
      case TYPE_UCHAR:
	ap->a.a_uchar = va_arg (args, /*unsigned char*/ int);
	break;
      case TYPE_SHORT:
	ap->a.a_short = va_arg (args, /*short*/ int);
	break;
      case TYPE_USHORT:
	ap->a.a_ushort = va_arg (args, /*unsigned short*/ int);
	break;
      case TYPE_INT:
	ap->a.a_int = va_arg (args, int);
	break;
      case TYPE_UINT:
	ap->a.a_uint = va_arg (args, unsigned int);
	break;
      case TYPE_LONGINT:
	ap->a.a_longint = va_arg (args, long int);
	break;
      case TYPE_ULONGINT:
	ap->a.a_ulongint = va_arg (args, unsigned long int);
	break;
#ifdef HAVE_LONG_LONG
      case TYPE_LONGLONGINT:
	ap->a.a_longlongint = va_arg (args, long long int);
	break;
      case TYPE_ULONGLONGINT:
	ap->a.a_ulonglongint = va_arg (args, unsigned long long int);
	break;
#endif
      case TYPE_DOUBLE:
	ap->a.a_double = va_arg (args, double);
	break;
#ifdef HAVE_LONG_DOUBLE
      case TYPE_LONGDOUBLE:
	ap->a.a_longdouble = va_arg (args, long double);
	break;
#endif
      case TYPE_CHAR:
	ap->a.a_char = va_arg (args, int);
	break;
#ifdef HAVE_WINT_T
      case TYPE_WIDE_CHAR:
	ap->a.a_wide_char = va_arg (args, wint_t);
	break;
#endif
      case TYPE_STRING:
	ap->a.a_string = va_arg (args, const char *);
	break;
#ifdef HAVE_WCHAR_T
      case TYPE_WIDE_STRING:
	ap->a.a_wide_string = va_arg (args, const wchar_t *);
	break;
#endif
      case TYPE_POINTER:
	ap->a.a_pointer = va_arg (args, void *);
	break;
      case TYPE_COUNT_SCHAR_POINTER:
	ap->a.a_count_schar_pointer = va_arg (args, signed char *);
	break;
      case TYPE_COUNT_SHORT_POINTER:
	ap->a.a_count_short_pointer = va_arg (args, short *);
	break;
      case TYPE_COUNT_INT_POINTER:
	ap->a.a_count_int_pointer = va_arg (args, int *);
	break;
      case TYPE_COUNT_LONGINT_POINTER:
	ap->a.a_count_longint_pointer = va_arg (args, long int *);
	break;
#ifdef HAVE_LONG_LONG
      case TYPE_COUNT_LONGLONGINT_POINTER:
	ap->a.a_count_longlongint_pointer = va_arg (args, long long int *);
	break;
#endif
      default:
	/* Unknown type.  */
	return -1;
      }
  return 0;
}

/* Get intmax_t.  */
#if HAVE_STDINT_H_WITH_UINTMAX
# include <stdint.h>
#endif
#if HAVE_INTTYPES_H_WITH_UINTMAX
# include <inttypes.h>
#endif

/* malloc(), realloc(), free().  */
#include <stdlib.h>

/* xsize.h -- Checked size_t computations. */

#ifndef _XSIZE_H
#define _XSIZE_H

/* Get SIZE_MAX.  */
#include <limits.h>
#if HAVE_STDINT_H
# include <stdint.h>
#endif

/* The size of memory objects is often computed through expressions of
   type size_t. Example:
      void* p = malloc (header_size + n * element_size).
   These computations can lead to overflow.  When this happens, malloc()
   returns a piece of memory that is way too small, and the program then
   crashes while attempting to fill the memory.
   To avoid this, the functions and macros in this file check for overflow.
   The convention is that SIZE_MAX represents overflow.
   malloc (SIZE_MAX) is not guaranteed to fail -- think of a malloc
   implementation that uses mmap --, it's recommended to use size_overflow_p()
   or size_in_bounds_p() before invoking malloc().
   The example thus becomes:
      size_t size = xsum (header_size, xtimes (n, element_size));
      void *p = (size_in_bounds_p (size) ? malloc (size) : NULL);
*/

/* Convert an arbitrary value >= 0 to type size_t.  */
#define xcast_size_t(N) \
  ((N) <= SIZE_MAX ? (size_t) (N) : SIZE_MAX)

/* Sum of two sizes, with overflow check.  */
static inline size_t
#if __GNUC__ >= 3
__attribute__ ((__pure__))
#endif
xsum (size_t size1, size_t size2) {
  size_t sum = size1 + size2;
  return (sum >= size1 ? sum : SIZE_MAX);
}

/* Sum of three sizes, with overflow check.  */
static inline size_t
#if __GNUC__ >= 3
__attribute__ ((__pure__))
#endif
xsum3 (size_t size1, size_t size2, size_t size3) {
  return xsum (xsum (size1, size2), size3);
}

/* Sum of four sizes, with overflow check.  */
static inline size_t
#if __GNUC__ >= 3
__attribute__ ((__pure__))
#endif
xsum4 (size_t size1, size_t size2, size_t size3, size_t size4) {
  return xsum (xsum (xsum (size1, size2), size3), size4);
}

/* Maximum of two sizes, with overflow check.  */
static inline size_t
#if __GNUC__ >= 3
__attribute__ ((__pure__))
#endif
xmax (size_t size1, size_t size2) {
  /* No explicit check is needed here, because for any n:
     max (SIZE_MAX, n) == SIZE_MAX and max (n, SIZE_MAX) == SIZE_MAX.  */
  return (size1 >= size2 ? size1 : size2);
}

/* Multiplication of a count with an element size, with overflow check.
   The count must be >= 0 and the element size must be > 0.
   This is a macro, not an inline function, so that it works correctly even
   when N is of a wider tupe and N > SIZE_MAX.  */
#define xtimes(N, ELSIZE) \
  ((N) <= SIZE_MAX / (ELSIZE) ? (size_t) (N) * (ELSIZE) : SIZE_MAX)

/* Check for overflow.  */
#define size_overflow_p(SIZE) \
  ((SIZE) == SIZE_MAX)
/* Check against overflow.  */
#define size_in_bounds_p(SIZE) \
  ((SIZE) != SIZE_MAX)

#endif /* _XSIZE_H */


#if WIDE_CHAR_VERSION
# define PRINTF_PARSE wprintf_parse
# define CHAR_T wchar_t
# define DIRECTIVE wchar_t_directive
# define DIRECTIVES wchar_t_directives
#else
# define PRINTF_PARSE printf_parse
# define CHAR_T char
# define DIRECTIVE char_directive
# define DIRECTIVES char_directives
#endif

#ifdef STATIC
STATIC
#endif
int
PRINTF_PARSE (const CHAR_T *format, DIRECTIVES *d, arguments *a) {
  const CHAR_T *cp = format;		/* pointer into format */
  size_t arg_posn = 0;		/* number of regular arguments consumed */
  size_t d_allocated;			/* allocated elements of d->dir */
  size_t a_allocated;			/* allocated elements of a->arg */
  size_t max_width_length = 0;
  size_t max_precision_length = 0;

  d->count = 0;
  d_allocated = 1;
  d->dir = malloc (d_allocated * sizeof (DIRECTIVE));
  if (d->dir == NULL)
    /* Out of memory.  */
    return -1;

  a->count = 0;
  a_allocated = 0;
  a->arg = NULL;

#define REGISTER_ARG(_index_,_type_) \
  {									\
    size_t n = (_index_);						\
    if (n >= a_allocated)						\
      {									\
	size_t memory_size;						\
	argument *memory;						\
									\
	a_allocated = xtimes (a_allocated, 2);				\
	if (a_allocated <= n)						\
	  a_allocated = xsum (n, 1);					\
	memory_size = xtimes (a_allocated, sizeof (argument));		\
	if (size_overflow_p (memory_size))				\
	  /* Overflow, would lead to out of memory.  */			\
	  goto error;							\
	memory = (a->arg						\
		  ? realloc (a->arg, memory_size)			\
		  : malloc (memory_size));				\
	if (memory == NULL)						\
	  /* Out of memory.  */						\
	  goto error;							\
	a->arg = memory;						\
      }									\
    while (a->count <= n)						\
      a->arg[a->count++].type = TYPE_NONE;				\
    if (a->arg[n].type == TYPE_NONE)					\
      a->arg[n].type = (_type_);					\
    else if (a->arg[n].type != (_type_))				\
      /* Ambiguous type for positional argument.  */			\
      goto error;							\
  }

  while (*cp != '\0')
    {
      CHAR_T c = *cp++;
      if (c == '%')
	{
	  size_t arg_index = ARG_NONE;
	  DIRECTIVE *dp = &d->dir[d->count];/* pointer to next directive */

	  /* Initialize the next directive.  */
	  dp->dir_start = cp - 1;
	  dp->flags = 0;
	  dp->width_start = NULL;
	  dp->width_end = NULL;
	  dp->width_arg_index = ARG_NONE;
	  dp->precision_start = NULL;
	  dp->precision_end = NULL;
	  dp->precision_arg_index = ARG_NONE;
	  dp->arg_index = ARG_NONE;

	  /* Test for positional argument.  */
	  if (*cp >= '0' && *cp <= '9') {
	      const CHAR_T *np;

	      for (np = cp; *np >= '0' && *np <= '9'; np++)
            ;
	      if (*np == '$') {
              size_t n = 0;

              for (np = cp; *np >= '0' && *np <= '9'; np++)
                n = xsum (xtimes (n, 10), *np - '0');
              if (n == 0)
                /* Positional argument 0.  */
                goto error;
              if (size_overflow_p (n))
                /* n too large, would lead to out of memory later.  */
                goto error;
              arg_index = n - 1;
              cp = np + 1;
          }
	    }

	  /* Read the flags.  */
	  for (;;) {
	      if (*cp == '\'') {
              dp->flags |= FLAG_GROUP;
              cp++;
		  } else if (*cp == '-') {
              dp->flags |= FLAG_LEFT;
              cp++;
		  } else if (*cp == '+') {
              dp->flags |= FLAG_SHOWSIGN;
              cp++;
		  } else if (*cp == ' ') {
              dp->flags |= FLAG_SPACE;
              cp++;
		  } else if (*cp == '#') {
              dp->flags |= FLAG_ALT;
              cp++;
		  } else if (*cp == '0') {
              dp->flags |= FLAG_ZERO;
              cp++;
		  } else
		      break;
      }

	  /* Parse the field width.  */
	  if (*cp == '*')
	    {
	      dp->width_start = cp;
	      cp++;
	      dp->width_end = cp;
	      if (max_width_length < 1)
		max_width_length = 1;

	      /* Test for positional argument.  */
	      if (*cp >= '0' && *cp <= '9')
		{
		  const CHAR_T *np;

		  for (np = cp; *np >= '0' && *np <= '9'; np++)
		    ;
		  if (*np == '$')
		    {
		      size_t n = 0;

		      for (np = cp; *np >= '0' && *np <= '9'; np++)
			n = xsum (xtimes (n, 10), *np - '0');
		      if (n == 0)
			/* Positional argument 0.  */
			goto error;
		      if (size_overflow_p (n))
			/* n too large, would lead to out of memory later.  */
			goto error;
		      dp->width_arg_index = n - 1;
		      cp = np + 1;
		    }
		}
	      if (dp->width_arg_index == ARG_NONE)
		{
		  dp->width_arg_index = arg_posn++;
		  if (dp->width_arg_index == ARG_NONE)
		    /* arg_posn wrapped around.  */
		    goto error;
		}
	      REGISTER_ARG (dp->width_arg_index, TYPE_INT);
	    }
	  else if (*cp >= '0' && *cp <= '9')
	    {
	      size_t width_length;

	      dp->width_start = cp;
	      for (; *cp >= '0' && *cp <= '9'; cp++)
		;
	      dp->width_end = cp;
	      width_length = dp->width_end - dp->width_start;
	      if (max_width_length < width_length)
		max_width_length = width_length;
	    }

	  /* Parse the precision.  */
	  if (*cp == '.')
	    {
	      cp++;
	      if (*cp == '*')
		{
		  dp->precision_start = cp - 1;
		  cp++;
		  dp->precision_end = cp;
		  if (max_precision_length < 2)
		    max_precision_length = 2;

		  /* Test for positional argument.  */
		  if (*cp >= '0' && *cp <= '9')
		    {
		      const CHAR_T *np;

		      for (np = cp; *np >= '0' && *np <= '9'; np++)
			;
		      if (*np == '$')
			{
			  size_t n = 0;

			  for (np = cp; *np >= '0' && *np <= '9'; np++)
			    n = xsum (xtimes (n, 10), *np - '0');
			  if (n == 0)
			    /* Positional argument 0.  */
			    goto error;
			  if (size_overflow_p (n))
			    /* n too large, would lead to out of memory
			       later.  */
			    goto error;
			  dp->precision_arg_index = n - 1;
			  cp = np + 1;
			}
		    }
		  if (dp->precision_arg_index == ARG_NONE)
		    {
		      dp->precision_arg_index = arg_posn++;
		      if (dp->precision_arg_index == ARG_NONE)
			/* arg_posn wrapped around.  */
			goto error;
		    }
		  REGISTER_ARG (dp->precision_arg_index, TYPE_INT);
		}
	      else
		{
		  size_t precision_length;

		  dp->precision_start = cp - 1;
		  for (; *cp >= '0' && *cp <= '9'; cp++)
		    ;
		  dp->precision_end = cp;
		  precision_length = dp->precision_end - dp->precision_start;
		  if (max_precision_length < precision_length)
		    max_precision_length = precision_length;
		}
	    }

	  {
	    arg_type type;

	    /* Parse argument type/size specifiers.  */
	    {
	      int flags = 0;

	      for (;;)
		{
		  if (*cp == 'h')
		    {
		      flags |= (1 << (flags & 1));
		      cp++;
		    }
		  else if (*cp == 'L')
		    {
		      flags |= 4;
		      cp++;
		    }
		  else if (*cp == 'l')
		    {
		      flags += 8;
		      cp++;
		    }
#ifdef HAVE_INTMAX_T
		  else if (*cp == 'j')
		    {
		      if (sizeof (intmax_t) > sizeof (long))
			{
			  /* intmax_t = long long */
			  flags += 16;
			}
		      else if (sizeof (intmax_t) > sizeof (int))
			{
			  /* intmax_t = long */
			  flags += 8;
			}
		      cp++;
		    }
#endif
		  else if (*cp == 'z' || *cp == 'Z')
		    {
		      /* 'z' is standardized in ISO C 99, but glibc uses 'Z'
			 because the warning facility in gcc-2.95.2 understands
			 only 'Z' (see gcc-2.95.2/gcc/c-common.c:1784).  */
		      if (sizeof (size_t) > sizeof (long))
			{
			  /* size_t = long long */
			  flags += 16;
			}
		      else if (sizeof (size_t) > sizeof (int))
			{
			  /* size_t = long */
			  flags += 8;
			}
		      cp++;
		    }
		  else if (*cp == 't')
		    {
		      if (sizeof (ptrdiff_t) > sizeof (long))
			{
			  /* ptrdiff_t = long long */
			  flags += 16;
			}
		      else if (sizeof (ptrdiff_t) > sizeof (int))
			{
			  /* ptrdiff_t = long */
			  flags += 8;
			}
		      cp++;
		    }
		  else
		    break;
		}

	      /* Read the conversion character.  */
	      c = *cp++;
	      switch (c)
		{
		case 'd': case 'i':
#ifdef HAVE_LONG_LONG
		  if (flags >= 16 || (flags & 4))
		    type = TYPE_LONGLONGINT;
		  else
#endif
		  if (flags >= 8)
		    type = TYPE_LONGINT;
		  else if (flags & 2)
		    type = TYPE_SCHAR;
		  else if (flags & 1)
		    type = TYPE_SHORT;
		  else
		    type = TYPE_INT;
		  break;
		case 'o': case 'u': case 'x': case 'X':
#ifdef HAVE_LONG_LONG
		  if (flags >= 16 || (flags & 4))
		    type = TYPE_ULONGLONGINT;
		  else
#endif
		  if (flags >= 8)
		    type = TYPE_ULONGINT;
		  else if (flags & 2)
		    type = TYPE_UCHAR;
		  else if (flags & 1)
		    type = TYPE_USHORT;
		  else
		    type = TYPE_UINT;
		  break;
		case 'f': case 'F': case 'e': case 'E': case 'g': case 'G':
		case 'a': case 'A':
#ifdef HAVE_LONG_DOUBLE
		  if (flags >= 16 || (flags & 4))
		    type = TYPE_LONGDOUBLE;
		  else
#endif
		  type = TYPE_DOUBLE;
		  break;
		case 'c':
		  if (flags >= 8)
#ifdef HAVE_WINT_T
		    type = TYPE_WIDE_CHAR;
#else
		    goto error;
#endif
		  else
		    type = TYPE_CHAR;
		  break;
#ifdef HAVE_WINT_T
		case 'C':
		  type = TYPE_WIDE_CHAR;
		  c = 'c';
		  break;
#endif
		case 's':
		  if (flags >= 8)
#ifdef HAVE_WCHAR_T
		    type = TYPE_WIDE_STRING;
#else
		    goto error;
#endif
		  else
		    type = TYPE_STRING;
		  break;
#ifdef HAVE_WCHAR_T
		case 'S':
		  type = TYPE_WIDE_STRING;
		  c = 's';
		  break;
#endif
		case 'p':
		  type = TYPE_POINTER;
		  break;
		case 'n':
#ifdef HAVE_LONG_LONG
		  if (flags >= 16 || (flags & 4))
		    type = TYPE_COUNT_LONGLONGINT_POINTER;
		  else
#endif
		  if (flags >= 8)
		    type = TYPE_COUNT_LONGINT_POINTER;
		  else if (flags & 2)
		    type = TYPE_COUNT_SCHAR_POINTER;
		  else if (flags & 1)
		    type = TYPE_COUNT_SHORT_POINTER;
		  else
		    type = TYPE_COUNT_INT_POINTER;
		  break;
		case '%':
		  type = TYPE_NONE;
		  break;
		default:
		  /* Unknown conversion character.  */
		  goto error;
		}
	    }

	    if (type != TYPE_NONE)
	      {
		dp->arg_index = arg_index;
		if (dp->arg_index == ARG_NONE)
		  {
		    dp->arg_index = arg_posn++;
		    if (dp->arg_index == ARG_NONE)
		      /* arg_posn wrapped around.  */
		      goto error;
		  }
		REGISTER_ARG (dp->arg_index, type);
	      }
	    dp->conversion = c;
	    dp->dir_end = cp;
	  }

	  d->count++;
	  if (d->count >= d_allocated)
	    {
	      size_t memory_size;
	      DIRECTIVE *memory;

	      d_allocated = xtimes (d_allocated, 2);
	      memory_size = xtimes (d_allocated, sizeof (DIRECTIVE));
	      if (size_overflow_p (memory_size))
		/* Overflow, would lead to out of memory.  */
		goto error;
	      memory = realloc (d->dir, memory_size);
	      if (memory == NULL)
		/* Out of memory.  */
		goto error;
	      d->dir = memory;
	    }
	}
    }
  d->dir[d->count].dir_start = cp;

  d->max_width_length = max_width_length;
  d->max_precision_length = max_precision_length;
  return 0;

error:
  if (a->arg)
    free (a->arg);
  if (d->dir)
    free (d->dir);
  return -1;
}

#undef DIRECTIVES
#undef DIRECTIVE
#undef CHAR_T
#undef PRINTF_PARSE


#ifndef IN_LIBINTL
# include <alloca.h>
#endif


#include <string.h>	/* memcpy(), strlen() */
#include <errno.h>	/* errno */
#include <float.h>	/* DBL_MAX_EXP, LDBL_MAX_EXP */

/* Some systems, like OSF/1 4.0 and Woe32, don't have EOVERFLOW.  */
#ifndef EOVERFLOW
# define EOVERFLOW E2BIG
#endif

#ifdef HAVE_WCHAR_T
# ifdef HAVE_WCSLEN
#  define local_wcslen wcslen
# else
   /* Solaris 2.5.1 has wcslen() in a separate library libw.so. To avoid
      a dependency towards this library, here is a local substitute.
      Define this substitute only once, even if this file is included
      twice in the same compilation unit.  */
#  ifndef local_wcslen_defined
#   define local_wcslen_defined 1
static size_t local_wcslen (const wchar_t *s) {
  const wchar_t *ptr;

  for (ptr = s; *ptr != (wchar_t) 0; ptr++)
    ;
  return ptr - s;
}
#  endif
# endif
#endif

#if WIDE_CHAR_VERSION
# define VASNPRINTF vasnwprintf
# define CHAR_T wchar_t
# define DIRECTIVE wchar_t_directive
# define DIRECTIVES wchar_t_directives
# define PRINTF_PARSE wprintf_parse
# define USE_SNPRINTF 1
# if HAVE_DECL__SNWPRINTF
   /* On Windows, the function swprintf() has a different signature than
      on Unix; we use the _snwprintf() function instead.  */
#  define SNPRINTF _snwprintf
# else
   /* Unix.  */
#  define SNPRINTF swprintf
# endif
#else
# define VASNPRINTF vasnprintf
# define CHAR_T char
# define DIRECTIVE char_directive
# define DIRECTIVES char_directives
# define PRINTF_PARSE printf_parse
# define USE_SNPRINTF (HAVE_DECL__SNPRINTF || HAVE_SNPRINTF)
# if HAVE_DECL__SNPRINTF
   /* Windows.  */
#  define SNPRINTF _snprintf
# else
   /* Unix.  */
#  define SNPRINTF snprintf
# endif
#endif

CHAR_T * VASNPRINTF (CHAR_T *resultbuf, size_t *lengthp, const CHAR_T *format, va_list args) {
  DIRECTIVES d;
  arguments a;

  if (PRINTF_PARSE (format, &d, &a) < 0) {
      errno = EINVAL;
      return NULL;
  }

#define CLEANUP() \
  free (d.dir);								\
  if (a.arg)								\
    free (a.arg);

  if (printf_fetchargs (args, &a) < 0) {
      CLEANUP ();
      errno = EINVAL;
      return NULL;
  }

  {
    size_t buf_neededlength;
    CHAR_T *buf;
    CHAR_T *buf_malloced;
    const CHAR_T *cp;
    size_t i;
    DIRECTIVE *dp;
    /* Output string accumulator.  */
    CHAR_T *result;
    size_t allocated;
    size_t length;

    /* Allocate a small buffer that will hold a directive passed to
       sprintf or snprintf.  */
    buf_neededlength =
      xsum4 (7, d.max_width_length, d.max_precision_length, 6);
#if HAVE_ALLOCA
    if (buf_neededlength < 4000 / sizeof (CHAR_T)) {
        buf = (CHAR_T *) alloca (buf_neededlength * sizeof (CHAR_T));
        buf_malloced = NULL;
    } else
#endif
    {
        size_t buf_memsize = xtimes (buf_neededlength, sizeof (CHAR_T));
        if (size_overflow_p (buf_memsize))
          goto out_of_memory_1;
        buf = (CHAR_T *) malloc (buf_memsize);
        if (buf == NULL)
          goto out_of_memory_1;
        buf_malloced = buf;
    }

    if (resultbuf != NULL) {
        result = resultbuf;
        allocated = *lengthp;
    } else {
        result = NULL;
        allocated = 0;
    }
    length = 0;
    /* Invariants:
       result is either == resultbuf or == NULL or malloc-allocated.
       If length > 0, then result != NULL.  */

    /* Ensures that allocated >= needed.  Aborts through a jump to
       out_of_memory if needed is SIZE_MAX or otherwise too big.  */
#define ENSURE_ALLOCATION(needed) \
    if ((needed) > allocated)						     \
      {									     \
	size_t memory_size;						     \
	CHAR_T *memory;							     \
									     \
	allocated = (allocated > 0 ? xtimes (allocated, 2) : 12);	     \
	if ((needed) > allocated)					     \
	  allocated = (needed);						     \
	memory_size = xtimes (allocated, sizeof (CHAR_T));		     \
	if (size_overflow_p (memory_size))				     \
	  goto out_of_memory;						     \
	if (result == resultbuf || result == NULL)			     \
	  memory = (CHAR_T *) malloc (memory_size);			     \
	else								     \
	  memory = (CHAR_T *) realloc (result, memory_size);		     \
	if (memory == NULL)						     \
	  goto out_of_memory;						     \
	if (result == resultbuf && length > 0)				     \
	  memcpy (memory, result, length * sizeof (CHAR_T));		     \
	result = memory;						     \
      }

    for (cp = format, i = 0, dp = &d.dir[0]; ; cp = dp->dir_end, i++, dp++)
      {
	if (cp != dp->dir_start)
	  {
	    size_t n = dp->dir_start - cp;
	    size_t augmented_length = xsum (length, n);

	    ENSURE_ALLOCATION (augmented_length);
	    memcpy (result + length, cp, n * sizeof (CHAR_T));
	    length = augmented_length;
	  }
	if (i == d.count)
	  break;

	/* Execute a single directive.  */
	if (dp->conversion == '%') {
	    size_t augmented_length;

	    if (!(dp->arg_index == ARG_NONE))
	      abort ();
	    augmented_length = xsum (length, 1);
	    ENSURE_ALLOCATION (augmented_length);
	    result[length] = '%';
	    length = augmented_length;
	} else {
	    if (!(dp->arg_index != ARG_NONE))
	      abort ();

	    if (dp->conversion == 'n') {
		switch (a.arg[dp->arg_index].type) {
		  case TYPE_COUNT_SCHAR_POINTER:
		    *a.arg[dp->arg_index].a.a_count_schar_pointer = length;
		    break;
		  case TYPE_COUNT_SHORT_POINTER:
		    *a.arg[dp->arg_index].a.a_count_short_pointer = length;
		    break;
		  case TYPE_COUNT_INT_POINTER:
		    *a.arg[dp->arg_index].a.a_count_int_pointer = length;
		    break;
		  case TYPE_COUNT_LONGINT_POINTER:
		    *a.arg[dp->arg_index].a.a_count_longint_pointer = length;
		    break;
#ifdef HAVE_LONG_LONG
		  case TYPE_COUNT_LONGLONGINT_POINTER:
		    *a.arg[dp->arg_index].a.a_count_longlongint_pointer = length;
		    break;
#endif
		  default:
		    abort ();
		  }
	      }
	    else {
		arg_type type = a.arg[dp->arg_index].type;
		CHAR_T *p;
		unsigned int prefix_count;
		int prefixes[2];
#if !USE_SNPRINTF
		size_t tmp_length;
		CHAR_T tmpbuf[700];
		CHAR_T *tmp;

		/* Allocate a temporary buffer of sufficient size for calling
		   sprintf.  */
		{
		  size_t width;
		  size_t precision;

		  width = 0;
		  if (dp->width_start != dp->width_end) {
		      if (dp->width_arg_index != ARG_NONE) {
                  int arg;

                  if (!(a.arg[dp->width_arg_index].type == TYPE_INT))
                    abort ();
                  arg = a.arg[dp->width_arg_index].a.a_int;
                  width = (arg < 0 ? (unsigned int) (-arg) : arg);
			  } else {
                  const CHAR_T *digitp = dp->width_start;

                  do
                    width = xsum (xtimes (width, 10), *digitp++ - '0');
                  while (digitp != dp->width_end);
			  }
		  }

		  precision = 6;
		  if (dp->precision_start != dp->precision_end) {
		      if (dp->precision_arg_index != ARG_NONE) {
                  int arg;

                  if (!(a.arg[dp->precision_arg_index].type == TYPE_INT))
                    abort ();
                  arg = a.arg[dp->precision_arg_index].a.a_int;
                  precision = (arg < 0 ? 0 : arg);
		      } else {
                  const CHAR_T *digitp = dp->precision_start + 1;

                  precision = 0;
                  while (digitp != dp->precision_end)
                    precision = xsum (xtimes (precision, 10), *digitp++ - '0');
              }
		  }

		  switch (dp->conversion) {

		    case 'd': case 'i': case 'u':
# ifdef HAVE_LONG_LONG
		      if (type == TYPE_LONGLONGINT || type == TYPE_ULONGLONGINT)
			tmp_length =
			  (unsigned int) (sizeof (unsigned long long) * CHAR_BIT
					  * 0.30103 /* binary -> decimal */
					  * 2 /* estimate for FLAG_GROUP */
					 )
			  + 1 /* turn floor into ceil */
			  + 1; /* account for leading sign */
		      else
# endif
		      if (type == TYPE_LONGINT || type == TYPE_ULONGINT)
			tmp_length =
			  (unsigned int) (sizeof (unsigned long) * CHAR_BIT
					  * 0.30103 /* binary -> decimal */
					  * 2 /* estimate for FLAG_GROUP */
					 )
			  + 1 /* turn floor into ceil */
			  + 1; /* account for leading sign */
		      else
			tmp_length =
			  (unsigned int) (sizeof (unsigned int) * CHAR_BIT
					  * 0.30103 /* binary -> decimal */
					  * 2 /* estimate for FLAG_GROUP */
					 )
			  + 1 /* turn floor into ceil */
			  + 1; /* account for leading sign */
		      break;

		    case 'o':
# ifdef HAVE_LONG_LONG
		      if (type == TYPE_LONGLONGINT || type == TYPE_ULONGLONGINT)
			tmp_length =
			  (unsigned int) (sizeof (unsigned long long) * CHAR_BIT
					  * 0.333334 /* binary -> octal */
					 )
			  + 1 /* turn floor into ceil */
			  + 1; /* account for leading sign */
		      else
# endif
		      if (type == TYPE_LONGINT || type == TYPE_ULONGINT)
			tmp_length =
			  (unsigned int) (sizeof (unsigned long) * CHAR_BIT
					  * 0.333334 /* binary -> octal */
					 )
			  + 1 /* turn floor into ceil */
			  + 1; /* account for leading sign */
		      else
			tmp_length =
			  (unsigned int) (sizeof (unsigned int) * CHAR_BIT
					  * 0.333334 /* binary -> octal */
					 )
			  + 1 /* turn floor into ceil */
			  + 1; /* account for leading sign */
		      break;

		    case 'x': case 'X':
# ifdef HAVE_LONG_LONG
		      if (type == TYPE_LONGLONGINT || type == TYPE_ULONGLONGINT)
			tmp_length =
			  (unsigned int) (sizeof (unsigned long long) * CHAR_BIT
					  * 0.25 /* binary -> hexadecimal */
					 )
			  + 1 /* turn floor into ceil */
			  + 2; /* account for leading sign or alternate form */
		      else
# endif
		      if (type == TYPE_LONGINT || type == TYPE_ULONGINT)
			tmp_length =
			  (unsigned int) (sizeof (unsigned long) * CHAR_BIT
					  * 0.25 /* binary -> hexadecimal */
					 )
			  + 1 /* turn floor into ceil */
			  + 2; /* account for leading sign or alternate form */
		      else
			tmp_length =
			  (unsigned int) (sizeof (unsigned int) * CHAR_BIT
					  * 0.25 /* binary -> hexadecimal */
					 )
			  + 1 /* turn floor into ceil */
			  + 2; /* account for leading sign or alternate form */
		      break;

		    case 'f': case 'F':
# ifdef HAVE_LONG_DOUBLE
		      if (type == TYPE_LONGDOUBLE)
			tmp_length =
			  (unsigned int) (LDBL_MAX_EXP
					  * 0.30103 /* binary -> decimal */
					  * 2 /* estimate for FLAG_GROUP */
					 )
			  + 1 /* turn floor into ceil */
			  + 10; /* sign, decimal point etc. */
		      else
# endif
			tmp_length =
			  (unsigned int) (DBL_MAX_EXP
					  * 0.30103 /* binary -> decimal */
					  * 2 /* estimate for FLAG_GROUP */
					 )
			  + 1 /* turn floor into ceil */
			  + 10; /* sign, decimal point etc. */
		      tmp_length = xsum (tmp_length, precision);
		      break;

		    case 'e': case 'E': case 'g': case 'G':
		    case 'a': case 'A':
		      tmp_length =
			12; /* sign, decimal point, exponent etc. */
		      tmp_length = xsum (tmp_length, precision);
		      break;

		    case 'c':
# if defined HAVE_WINT_T && !WIDE_CHAR_VERSION
		      if (type == TYPE_WIDE_CHAR)
			tmp_length = MB_CUR_MAX;
		      else
# endif
			tmp_length = 1;
		      break;

		    case 's':
# ifdef HAVE_WCHAR_T
		      if (type == TYPE_WIDE_STRING) {
			  tmp_length =
			    local_wcslen (a.arg[dp->arg_index].a.a_wide_string);

#  if !WIDE_CHAR_VERSION
			  tmp_length = xtimes (tmp_length, MB_CUR_MAX);
#  endif
			}
		      else
# endif
			tmp_length = strlen (a.arg[dp->arg_index].a.a_string);
		      break;

		    case 'p':
		      tmp_length =
			(unsigned int) (sizeof (void *) * CHAR_BIT
					* 0.25 /* binary -> hexadecimal */
				       )
			  + 1 /* turn floor into ceil */
			  + 2; /* account for leading 0x */
		      break;

		    default:
		      abort ();
		    }

		  if (tmp_length < width)
		    tmp_length = width;

		  tmp_length = xsum (tmp_length, 1); /* account for trailing NUL */
		}

		if (tmp_length <= sizeof (tmpbuf) / sizeof (CHAR_T))
		  tmp = tmpbuf;
		else {
		    size_t tmp_memsize = xtimes (tmp_length, sizeof (CHAR_T));

		    if (size_overflow_p (tmp_memsize))
		      /* Overflow, would lead to out of memory.  */
		      goto out_of_memory;
		    tmp = (CHAR_T *) malloc (tmp_memsize);
		    if (tmp == NULL)
		      /* Out of memory.  */
		      goto out_of_memory;
		}
#endif

		/* Construct the format string for calling snprintf or
		   sprintf.  */
		p = buf;
		*p++ = '%';
		if (dp->flags & FLAG_GROUP)
		  *p++ = '\'';
		if (dp->flags & FLAG_LEFT)
		  *p++ = '-';
		if (dp->flags & FLAG_SHOWSIGN)
		  *p++ = '+';
		if (dp->flags & FLAG_SPACE)
		  *p++ = ' ';
		if (dp->flags & FLAG_ALT)
		  *p++ = '#';
		if (dp->flags & FLAG_ZERO)
		  *p++ = '0';
		if (dp->width_start != dp->width_end) {
		    size_t n = dp->width_end - dp->width_start;
		    memcpy (p, dp->width_start, n * sizeof (CHAR_T));
		    p += n;
		}
		if (dp->precision_start != dp->precision_end) {
		    size_t n = dp->precision_end - dp->precision_start;
		    memcpy (p, dp->precision_start, n * sizeof (CHAR_T));
		    p += n;
		}

		switch (type) {
#ifdef HAVE_LONG_LONG
		  case TYPE_LONGLONGINT:
		  case TYPE_ULONGLONGINT:
		    *p++ = 'l';
		    /*FALLTHROUGH*/
#endif
		  case TYPE_LONGINT:
		  case TYPE_ULONGINT:
#ifdef HAVE_WINT_T
		  case TYPE_WIDE_CHAR:
#endif
#ifdef HAVE_WCHAR_T
		  case TYPE_WIDE_STRING:
#endif
		    *p++ = 'l';
		    break;
#ifdef HAVE_LONG_DOUBLE
		  case TYPE_LONGDOUBLE:
		    *p++ = 'L';
		    break;
#endif
		  default:
		    break;
		  }
		*p = dp->conversion;
#if USE_SNPRINTF
		p[1] = '%';
		p[2] = 'n';
		p[3] = '\0';
#else
		p[1] = '\0';
#endif

		/* Construct the arguments for calling snprintf or sprintf.  */
		prefix_count = 0;
		if (dp->width_arg_index != ARG_NONE) {
		    if (!(a.arg[dp->width_arg_index].type == TYPE_INT))
		      abort ();
		    prefixes[prefix_count++] = a.arg[dp->width_arg_index].a.a_int;
	    }
		if (dp->precision_arg_index != ARG_NONE) {
		    if (!(a.arg[dp->precision_arg_index].type == TYPE_INT))
		      abort ();
		    prefixes[prefix_count++] = a.arg[dp->precision_arg_index].a.a_int;
	    }

#if USE_SNPRINTF
		/* Prepare checking whether snprintf returns the count
		   via %n.  */
		ENSURE_ALLOCATION (xsum (length, 1));
		result[length] = '\0';
#endif

		for (;;) {
		    size_t maxlen = allocated - length;
		    int count = -1;

#if USE_SNPRINTF
		    int retcount = 0;
# define SNPRINTF_BUF(arg) \
		    switch (prefix_count)				    \
		      {							    \
		      case 0:						    \
			retcount = SNPRINTF (result + length, maxlen, buf,  \
					     arg, &count);		    \
			break;						    \
		      case 1:						    \
			retcount = SNPRINTF (result + length, maxlen, buf,  \
					     prefixes[0], arg, &count);	    \
			break;						    \
		      case 2:						    \
			retcount = SNPRINTF (result + length, maxlen, buf,  \
					     prefixes[0], prefixes[1], arg, \
					     &count);			    \
			break;						    \
		      default:						    \
			abort ();					    \
		      }
#else
# define SNPRINTF_BUF(arg) \
		    switch (prefix_count)				    \
		      {							    \
		      case 0:						    \
			count = sprintf (tmp, buf, arg);		    \
			break;						    \
		      case 1:						    \
			count = sprintf (tmp, buf, prefixes[0], arg);	    \
			break;						    \
		      case 2:						    \
			count = sprintf (tmp, buf, prefixes[0], prefixes[1],\
					 arg);				    \
			break;						    \
		      default:						    \
			abort ();					    \
		      }
#endif

		    switch (type) {
		      case TYPE_SCHAR: {
			  int arg = a.arg[dp->arg_index].a.a_schar;
			  SNPRINTF_BUF (arg);
			}
			break;
		      case TYPE_UCHAR: {
			  unsigned int arg = a.arg[dp->arg_index].a.a_uchar;
			  SNPRINTF_BUF (arg);
			}
			break;
		      case TYPE_SHORT: {
			  int arg = a.arg[dp->arg_index].a.a_short;
			  SNPRINTF_BUF (arg);
			}
			break;
		      case TYPE_USHORT: {
			  unsigned int arg = a.arg[dp->arg_index].a.a_ushort;
			  SNPRINTF_BUF (arg);
			}
			break;
		      case TYPE_INT: {
			  int arg = a.arg[dp->arg_index].a.a_int;
			  SNPRINTF_BUF (arg);
			}
			break;
		      case TYPE_UINT: {
			  unsigned int arg = a.arg[dp->arg_index].a.a_uint;
			  SNPRINTF_BUF (arg);
			}
			break;
		      case TYPE_LONGINT: {
			  long int arg = a.arg[dp->arg_index].a.a_longint;
			  SNPRINTF_BUF (arg);
			}
			break;
		      case TYPE_ULONGINT: {
			  unsigned long int arg = a.arg[dp->arg_index].a.a_ulongint;
			  SNPRINTF_BUF (arg);
			}
			break;
#ifdef HAVE_LONG_LONG
		      case TYPE_LONGLONGINT: {
			  long long int arg = a.arg[dp->arg_index].a.a_longlongint;
			  SNPRINTF_BUF (arg);
			}
			break;
		      case TYPE_ULONGLONGINT: {
			  unsigned long long int arg = a.arg[dp->arg_index].a.a_ulonglongint;
			  SNPRINTF_BUF (arg);
			}
			break;
#endif
		      case TYPE_DOUBLE: {
			  double arg = a.arg[dp->arg_index].a.a_double;
			  SNPRINTF_BUF (arg);
			}
			break;
#ifdef HAVE_LONG_DOUBLE
		      case TYPE_LONGDOUBLE: {
			  long double arg = a.arg[dp->arg_index].a.a_longdouble;
			  SNPRINTF_BUF (arg);
			}
			break;
#endif
		      case TYPE_CHAR: {
			  int arg = a.arg[dp->arg_index].a.a_char;
			  SNPRINTF_BUF (arg);
			}
			break;
#ifdef HAVE_WINT_T
		      case TYPE_WIDE_CHAR: {
			  wint_t arg = a.arg[dp->arg_index].a.a_wide_char;
			  SNPRINTF_BUF (arg);
			}
			break;
#endif
		      case TYPE_STRING: {
			  const char *arg = a.arg[dp->arg_index].a.a_string;
			  SNPRINTF_BUF (arg);
			}
			break;
#ifdef HAVE_WCHAR_T
		      case TYPE_WIDE_STRING: {
			  const wchar_t *arg = a.arg[dp->arg_index].a.a_wide_string;
			  SNPRINTF_BUF (arg);
			}
			break;
#endif
		      case TYPE_POINTER: {
			  void *arg = a.arg[dp->arg_index].a.a_pointer;
			  SNPRINTF_BUF (arg);
			}
			break;
		      default:
			abort ();
	      }

#if USE_SNPRINTF
		    /* Portability: Not all implementations of snprintf()
		       are ISO C 99 compliant.  Determine the number of
		       bytes that snprintf() has produced or would have
		       produced.  */
		    if (count >= 0) {
                /* Verify that snprintf() has NUL-terminated its
                   result.  */
                if (count < maxlen && result[length + count] != '\0')
                  abort ();
                /* Portability hack.  */
                if (retcount > count)
                  count = retcount;
		    } else {
                /* snprintf() doesn't understand the '%n'
                   directive.  */
                if (p[1] != '\0') {
                    /* Don't use the '%n' directive; instead, look
                       at the snprintf() return value.  */
                    p[1] = '\0';
                    continue;
                } else {
                    /* Look at the snprintf() return value.  */
                    if (retcount < 0) {
                    /* HP-UX 10.20 snprintf() is doubly deficient:
                       It doesn't understand the '%n' directive,
                       *and* it returns -1 (rather than the length
                       that would have been required) when the
                       buffer is too small.  */
                    size_t bigger_need =
                      xsum (xtimes (allocated, 2), 12);
                    ENSURE_ALLOCATION (bigger_need);
                    continue;
                    } else
                      count = retcount;
               }
		    }
#endif

		    /* Attempt to handle failure.  */
		    if (count < 0) {
                if (!(result == resultbuf || result == NULL))
                  free (result);
                if (buf_malloced != NULL)
                  free (buf_malloced);
                CLEANUP ();
                errno = EINVAL;
                return NULL;
		    }

#if !USE_SNPRINTF
		    if (count >= tmp_length)
		      /* tmp_length was incorrectly calculated - fix the
			 code above!  */
		      abort ();
#endif

		    /* Make room for the result.  */
		    if (count >= maxlen)
		      {
			/* Need at least count bytes.  But allocate
			   proportionally, to avoid looping eternally if
			   snprintf() reports a too small count.  */
			size_t n =
			  xmax (xsum (length, count), xtimes (allocated, 2));

			ENSURE_ALLOCATION (n);
#if USE_SNPRINTF
			continue;
#endif
		      }

#if USE_SNPRINTF
		    /* The snprintf() result did fit.  */
#else
		    /* Append the sprintf() result.  */
		    memcpy (result + length, tmp, count * sizeof (CHAR_T));
		    if (tmp != tmpbuf)
		      free (tmp);
#endif

		    length += count;
		    break;
		  }
	      }
	  }
      }

    /* Add the final NUL.  */
    ENSURE_ALLOCATION (xsum (length, 1));
    result[length] = '\0';

    if (result != resultbuf && length + 1 < allocated)
      {
	/* Shrink the allocated memory if possible.  */
	CHAR_T *memory;

	memory = (CHAR_T *) realloc (result, (length + 1) * sizeof (CHAR_T));
	if (memory != NULL)
	  result = memory;
      }

    if (buf_malloced != NULL)
      free (buf_malloced);
    CLEANUP ();
    *lengthp = length;
    if (length > INT_MAX)
      goto length_overflow;
    return result;

  length_overflow:
    /* We could produce such a big string, but its length doesn't fit into
       an 'int'.  POSIX says that snprintf() fails with errno = EOVERFLOW in
       this case.  */
    if (result != resultbuf)
      free (result);
    errno = EOVERFLOW;
    return NULL;

  out_of_memory:
    if (!(result == resultbuf || result == NULL))
      free (result);
    if (buf_malloced != NULL)
      free (buf_malloced);
  out_of_memory_1:
    CLEANUP ();
    errno = ENOMEM;
    return NULL;
  }
}

#undef SNPRINTF
#undef USE_SNPRINTF
#undef PRINTF_PARSE
#undef DIRECTIVES
#undef DIRECTIVE
#undef CHAR_T
#undef VASNPRINTF

int vasprintf (char **resultp, const char *format, va_list args) {
  size_t length;
  char *result = vasnprintf (NULL, &length, format, args);
  if (result == NULL)
    return -1;

  *resultp = result;
  /* Return the number of resulting bytes, excluding the trailing NUL.
     If it wouldn't fit in an 'int', vasnprintf() would have returned NULL
     and set errno to EOVERFLOW.  */
  return length;
}
#endif
    /** \endcond */ //Doxygen ignore.
