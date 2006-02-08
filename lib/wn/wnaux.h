/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

AUTHOR:

  Will Naylor

****************************************************************************/
#ifndef wnportH
#define wnportH

#define wn_min(_x,_y)      (  ((_x)<(_y)) ? (_x) : (_y)  )
#define wn_max(_x,_y)      (  ((_x)>(_y)) ? (_x) : (_y)  )


#define WN_IHUGE (1000000000)

#define WN_FHUGE (1.0e+50)
#define WN_FTINY (1.0/WN_FHUGE)


#ifdef _LANGUAGE_C_PLUS_PLUS
#  ifndef __cplusplus
#    define __cplusplus
#  endif
#endif


/*
 * EXTERN is used to allow header files from ANSI C modules to be used by
 * modules compiled under C++
 */
#ifndef EXTERN
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

/*
 * This is used to bracket code that may or may be in C or C++, forcing it
 * to have C linkage. For example, you may have a header file containing a
 * typedef for a pointer to a function. Put EXTERN_BEGIN and EXTERN_END
 * around the typedef to force it to be a pointer to a "C" function, whether
 * it is included in a C file or a C++ file.
 */
#ifdef __cplusplus
#define WN_EXTERN_BEGIN extern "C" {
#define WN_EXTERN_END   };
#else
#define WN_EXTERN_BEGIN
#define WN_EXTERN_END
#endif


#if defined(unix)  || defined(UNIX)      || defined(__unix__)   || \
    defined(linux) || defined(__linux__) || defined(__CYGWIN__) || \
    defined(hpux)  || defined(__hpux__)  || defined(__APPLE__)
# define WN_UNIX
#endif

#if defined(vms)
# define WN_VMS
#endif

#if defined(WIN32)
# define WN_WINDOWS
#endif

#if !defined(WN_UNIX) && !defined(WN_VMS) && !defined(WN_WINDOWS)
# error Unrecognized OS
#endif
#if defined(WN_UNIX) && defined(WN_VMS)
# error Unix and VMS
#endif
#if defined(WN_UNIX) && defined(WN_WINDOWS)
# error Unix and Windows
#endif
#if defined(WN_VMS) && defined(WN_WINDOWS)
# error VMS and Windows
#endif

#if defined(WN_WINDOWS)
# define WN_OSSLASH ("\\")
#else
# define WN_OSSLASH ("/")
#endif

/*   say WN_DEPRECATED after a routine extern declaration in the .h file, this
** will generate a compiler warning every time it is called. */
/*   note gcc on sparc doesn't understand this attribute */
#if defined(__GNUC__) && !defined(sparc) && !defined(__ia64__)
# define WN_DEPRECATED __attribute__ ((__deprecated__))
#else
# define WN_DEPRECATED
#endif


#endif /* wnportH */
/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

AUTHOR:

  Will Naylor

****************************************************************************/
#ifndef wnconsH
#define wnconsH



#ifdef local
#undef local
#endif
#define local static


/*	bool is a predefined type in C++, that behaves differently from an
**  int.  For historical and possibly temporary reasons, we are defining
**  it to be an int in C++ as well as C.  This may change.
**	If you are using the STL and possibly other C++ packages, redefining
**  bool to an int will screw things up.  To avoid this, compile with
**  "-DWN_LEAVE_BOOL_ALONE" set on your command line.  Wnlib has been
**  tested using natural C++ bools, and will work.
**	bool must be defined for C, because wnlib uses the type 'bool' all
**  over the place, and C doesn't define a bool type. */

#if !defined(__cplusplus) || !defined(WN_LEAVE_BOOL_ALONE)
# ifdef   bool
#   undef bool
# endif
# define  bool int

# ifdef __STL_CONFIG_H
#   error wnlib: Redefining bool, which can interfere badly with STL
#   error     templates.  See low/wncons.h on how to recover.
# endif
#endif

#ifdef TRUE
#undef TRUE
#endif
#define TRUE    (1)

#ifdef FALSE
#undef FALSE
#endif
#define FALSE   (0)


typedef void *ptr;


#ifndef NULL
#  ifdef __cplusplus
#    define NULL (0)
#  else
#    define NULL ((void *)0)
#  endif
#endif

#endif
/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

AUTHOR:

  Will Naylor

****************************************************************************/
#ifndef wnerrH
#define wnerrH



#define WN_SUBOPTIMAL     1
#define WN_SUCCESS        0
#define WN_SINGULAR       (-1)
#define WN_UNBOUNDED      (-2)
#define WN_BAD_ARGS       (-3)
#define WN_INFEASIBLE     (-4)




#endif


/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

AUTHOR:

  Will Naylor

****************************************************************************/
#ifndef wnconsH
#define wnconsH



#ifdef local
#undef local
#endif
#define local static


/*	bool is a predefined type in C++, that behaves differently from an
**  int.  For historical and possibly temporary reasons, we are defining
**  it to be an int in C++ as well as C.  This may change.
**	If you are using the STL and possibly other C++ packages, redefining
**  bool to an int will screw things up.  To avoid this, compile with
**  "-DWN_LEAVE_BOOL_ALONE" set on your command line.  Wnlib has been
**  tested using natural C++ bools, and will work.
**	bool must be defined for C, because wnlib uses the type 'bool' all
**  over the place, and C doesn't define a bool type. */

#if !defined(__cplusplus) || !defined(WN_LEAVE_BOOL_ALONE)
# ifdef   bool
#   undef bool
# endif
# define  bool int

# ifdef __STL_CONFIG_H
#   error wnlib: Redefining bool, which can interfere badly with STL
#   error     templates.  See low/wncons.h on how to recover.
# endif
#endif

#ifdef TRUE
#undef TRUE
#endif
#define TRUE    (1)

#ifdef FALSE
#undef FALSE
#endif
#define FALSE   (0)


typedef void *ptr;


#ifndef NULL
#  ifdef __cplusplus
#    define NULL (0)
#  else
#    define NULL ((void *)0)
#  endif
#endif

#endif
/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

****************************************************************************/

#ifndef wnparvectH
#define wnparvectH

typedef void *(*wn_parallel_func)(void*);
typedef void **wn_parallel_args;
typedef void (*wn_parallel_callback)(wn_parallel_func, wn_parallel_args,
					char *debug_name);
typedef void (*wn_parallel_vect_operation_func)(int, int, int, void *);
typedef struct wn_parvect_context_struct *wn_parvect_context;

EXTERN
wn_parvect_context wn_parvect_init(wn_parallel_callback par_cb, 
				   int num_threads, int num_vars);
EXTERN void wn_parvect_free(wn_parvect_context context);

EXTERN
void wn_parallelize_vector_operation(wn_parvect_context context,
	  wn_parallel_vect_operation_func operation_func, 
	  void *user_data, int len, char *operation_name);

EXTERN
void wn_add_vect_and_scaled_vect_par(wn_parvect_context context,
				     double *out, double *v1, double *v2,
				     double s2, int len);
EXTERN
void wn_add_scaled_vects_par(wn_parvect_context context,
			     double *out, double *v1, double s1, 
			     double *v2, double s2, int len);
EXTERN
void wn_add_3_scaled_vects_par(wn_parvect_context context,
			       double *out, double *v1, double s1,
			       double *v2, double s2,
			       double *v3, double s3, int len);
EXTERN
double wn_norm2_vect_par(wn_parvect_context context,
			 double *v1, int len);
EXTERN
double wn_dot_vects_par(wn_parvect_context context,
			double *v1, double *v2, int len);
EXTERN
void wn_multiply_vect_by_vect_par(wn_parvect_context context,
				  double *v1, double *v2, int len);
EXTERN
void wn_divide_vect_by_vect_par(wn_parvect_context context,
				double *v1, double *v2, int len);
EXTERN
void wn_zero_vect_par(wn_parvect_context context, double *v1, int len);
EXTERN
void wn_copy_vect_par(wn_parvect_context context,
		      double *v1, double *v2, int len);
EXTERN
double wn_scaled_max_diff_vect_par(wn_parvect_context context, 
				   double *v1, double *v2, double *s, int len);
EXTERN
void wn_make_vect_par(wn_parvect_context context, double **v1, int len);

EXTERN
void wn_add_scaled_vects_and_dot_and_norm_par(wn_parvect_context context,
                             double *out,
                             double *pnorm2_out, double *pdot_v1_out,
			     double *v1, double s1, 
			     double *v2, double s2, int len);

#endif
/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

AUTHOR:

  Will Naylor

****************************************************************************/
#ifndef wnconjdH
#define wnconjdH





#define WN_CONJ_DIR_DBG_NONE        0
#define WN_CONJ_DIR_DBG_PASSES      1
#define WN_CONJ_DIR_DBG_LINESEARCH  2
#define WN_CONJ_DIR_DBG_ALL         3

extern int wn_conj_direction_debug;



EXTERN void wn_conj_direction_method
(
  int *pcode,
  double *pval_min,
  double vect[],
  double initial_coord_x0s[],
  int len,
  double (*pfunction)(double vect[]),
  int max_iterations
);
EXTERN void wn_force_conj_direction_stop(void);



#endif





#ifndef wnmbtrH
#define wnmbtrH


/*
#if !defined(wnmbtrC) && !defined(wnmemgC) && !defined(wnmbtsC)
# error wnmbtr.h included by file other than wnmbtr.c, wnmemg.c, and wnmbts.c
#endif
*/




#define WN_MB_MIN    (-3)
#define WN_MB_MAX    3
#define WN_MB_LT     (-2)
#define WN_MB_GT     2
#define WN_MB_LE     (-1)
#define WN_MB_GE     1
#define WN_MB_EQ     0


typedef struct wn_mbtree_struct *wn_mbtree;
typedef struct wn_mbhandle_struct *wn_mbhandle;
typedef struct wn_free_block_struct *wn_free_block;

struct wn_mbtree_struct
{
  wn_mbhandle handle_tree;

  int (*pcompare_keys_func)(ptr key1,ptr key2);
};

/* A free block consists of a size field, always aligned such that its
** address & 7 == 4, followed by a *wn_handle, which will naturally be
** 8-byte aligned.  We cannot declare one struct like this because on
** 64-bit machines, the compiler will assume the struct is 8 byte aligned
** and pad between the fields to align the pointers accordingly.
**
** while ideally, we would just consolidate the size field with the
** mbhandle field, but as has just been explained, C and C++ won't
** cooperate.  So what we do is declare it in two separate structs and
** provide macros to get between the two. */

struct wn_mbhandle_struct
{
  wn_free_block next_free_block;	/* it is imperative that this
  **	be the 1st field, small free blocks will use this field but
  **	not the rest of the struct */

  wn_free_block free_block_list;	/* This is a header of a list
  **	of all the free blocks of this exact size, not including this one */

  ptr key;	/* size of the segment.  Redundant with the free_block_size
  **	field in the *wn_free_block, easier to leave it in */

  int level, count;

  wn_mbhandle left_child, right_child, parent;
};

struct wn_free_block_struct
{
  int free_block_size;
};

#define WN_FREE_BLOCK_TO_MBHANDLE(free_block) \
	(TRUE ? (wn_mbhandle) ((free_block) + 1) : \
	/**/    (wn_mbhandle)  (free_block)->free_block_size)
	/* the unreachable case here is to provide type checking on
	** free_block, the optimizer will get rid of the code */

#define WN_MBHANDLE_TO_FREE_BLOCK(mbhandle) \
	(TRUE ? ((wn_free_block) mbhandle) - 1 : mbhandle->free_block_list)
	/* again, the unreachable case is just for type checking */

/* the minimum size of a small free_block struct */
#define WN_SIZEOF_FREE_BLOCK_STRUCT	(sizeof(ptr) == 8 ? 12 : \
					(sizeof(ptr) == 4 ?  8 : \
						(wn_assert_notreached(), 0)))

EXTERN void wn_mmkptrbtree(wn_mbtree tree);
EXTERN void wn_mbget(wn_mbhandle *phandle,wn_mbtree tree,ptr key,int compare);
EXTERN void wn_mbins(wn_mbhandle   handle,wn_mbtree tree, ptr key);
EXTERN void wn_mbdel(wn_mbhandle handle,wn_mbtree tree);
EXTERN int wn_mbcount(wn_mbtree tree);

/* don't EXTERN this one - it confuses some (but not all) sparc compilers */
#if !defined(sparc)
  EXTERN void wn_mbact(wn_mbtree tree, void (*paction)(wn_mbhandle handle),
		  ptr low_key,int low_compare,ptr high_key,int high_compare);
#else
	 void wn_mbact(wn_mbtree tree, void (*paction)(wn_mbhandle handle),
		  ptr low_key,int low_compare,ptr high_key,int high_compare);
#endif

EXTERN void wn_mbverify(wn_mbtree tree);

#endif /* wnmbtrH */


#define wnprmH


#define wn_permute_array(_result,_permutation,_array,_size)\
{\
	  int _i,_size_copy;\
	  \
	  _size_copy = (_size);\
	  \
	  wn_assert((_result) != (_array));\
	  wn_assert((int *)(_result) != (_permutation));\
	  \
	  for(_i=0;_i<_size_copy;_i++)\
	  {\
		      (_result)[_i] = (_array)[(_permutation)[_i]];\
		    }\
}


void wn_make_mat(double ***pmat,int len_i,int len_j);
void wn_transpose_mat(double **to_mat,double **from_mat,int len_i,int len_j);
void wn_mult_mats(double **result_mat,double **mat1,double **mat2, int len_i,int len_j,int len_k);
void wn_mult_mat_by_vect(double *result_vect,double **mat,double *vect, int len_i,int len_j);
void wn_copy_mat(double **out_mat,double **in_mat,int len_i,int len_j);
void wn_invert_mat(int *pcode,double **mat,int len_i);
double wn_eval_poly(double x,double coefs[],int len);
