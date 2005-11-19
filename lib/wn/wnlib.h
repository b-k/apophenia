#include <apophenia/headers.h>
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
#ifndef wnlibH
#define wnlibH

#include "wnaux.h"

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
#ifndef wnswapH
#define wnswapH

#define wn_square(_x)      ((_x)*(_x))



#define wn_swap(_a,_b,_type)\
{\
  _type _tmp;\
  \
  _tmp = (_a);\
  (_a) = (_b);\
  (_b) = _tmp;\
}


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

  Will Naylor, Bill Chapman

****************************************************************************/
#ifndef wnasrtH
#define wnasrtH




#if defined(linux) || defined(__linux__)
# define wn_assert_get_string(expr)	__STRING(expr)
#else
# if defined(__STDC__) || defined(WN_WINDOWS)
    /*   the || WN_WINDOWS is because MSVC++ 4.0 doesn't define __STDC__,
    ** yet the preprocessor behaves like __STDC__ */
#   define wn_assert_get_string(expr)	#expr
# else
#   define wn_assert_get_string(expr)	"expr"
# endif
#endif


#if defined(linux) || defined(__linux__)
# if defined(__cplusplus) ? __GNUC_PREREQ (2, 6) : __GNUC_PREREQ (2, 4)
#   define   wn_assert_function_name	__PRETTY_FUNCTION__
# else
#   if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#     define wn_assert_function_name	__func__
#   else
#     define wn_assert_function_name	NULL
#   endif
# endif
#else
# define     wn_assert_function_name	NULL
#endif


#define wn_assert(_cond) \
  ((void) ((_cond) ? (void)0 : wn_assert_routine_func_exp(__FILE__, \
  /**/	__LINE__, wn_assert_function_name, wn_assert_get_string(_cond))))

#define wn_assert_notreached()  wn_assert_notreached_routine_func(__FILE__, \
/**/					__LINE__, wn_assert_function_name)

#define wn_assert_warn(_cond) \
  ((void) ((_cond) ? (void) 0 : wn_assert_warn_routine_func_exp(__FILE__, \
  /**/	__LINE__, wn_assert_function_name, wn_assert_get_string(_cond))))

#define wn_assert_warn_notreached() \
  wn_assert_warn_notreached_r_func(__FILE__, __LINE__, \
  /**/						wn_assert_function_name)

#ifdef WN_FAST
# define wn_fast_assert(_cond)			((void) TRUE)
# define wn_fast_assert_notreached()		((void) TRUE)
# define wn_fast_assert_warn(_cond)		((void) TRUE)
# define wn_fast_assert_warn_notreached()	((void) TRUE)
#else
# define wn_fast_assert(_cond)			wn_assert(_cond)
# define wn_fast_assert_notreached()		wn_assert_notreached()
# define wn_fast_assert_warn(_cond)		wn_assert_warn(_cond)
# define wn_fast_assert_warn_notreached()	wn_assert_warn_notreached()
#endif


WN_EXTERN_BEGIN

extern void wn_abort(void);

extern void wn_set_assert_print
(
  void (*passert_print)(const char string[])
);
extern void wn_default_assert_print(const char string[]);

extern void wn_set_assert_custom_print
(
 void (*passert_custom_print)(const char file_name[], int line_num,
			      const char *func_name, const char *exp_string,
			      bool warn_only, bool notreached)
);
extern void wn_default_assert_custom_print(const char file_name[], 
			      int line_num, const char *func_name, 
                              const char *exp_string,
			      bool warn_only, bool notreached);

extern void wn_set_assert_crash
(
  void (*passert_crash)(void)
);
extern void wn_default_assert_crash(void);

extern void wn_assert_routine(const char file_name[],int line_num);

extern void wn_assert_notreached_routine(const char file_name[],
/**/						int line_num);

extern void wn_assert_warn_routine(const char file_name[], int line_num);

extern void wn_assert_warn_notreached_r(const char file_name[],int line_num);


extern void wn_assert_routine_func_exp(const char file_name[], int line_num,
/**/		const char *func_name, const char *exp_string);

extern void wn_assert_notreached_routine_func(const char *file_name,
/**/			int line_num, const char *func_name);

extern void wn_assert_warn_routine_func_exp(const char file_name[],
/**/	int line_num, const char *func_name, const char *exp_string);

extern void wn_assert_warn_notreached_r_func(const char file_name[],
/**/				int line_num, const char *func_name);

extern void wn_get_assert_override_string(char *output_buffer,
/**/				const char *file_name, int line_num);

WN_EXTERN_END

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
#ifndef wnabsH
#define wnabsH


#define wn_abs(_x)         (  ((_x)<0) ? (-(_x)) : (_x)  ) 

#define wn_sign(_x)        (  ((_x)>0) ? 1 : ( ((_x)==0) ? 0:(-1) )  ) 


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
#ifndef wnmemH
#define wnmemH


#ifndef _STDIO_H
  /* needed for size_t */
#endif
#ifndef wnlibH
#endif

WN_EXTERN_BEGIN

#define WN_NO_FREE           0x1
#define WN_GENERAL_FREE      0x2


typedef struct wn_mem_block_struct *wn_mem_block;
typedef struct wn_memgp_struct *wn_memgp;

struct wn_memgp_struct
{
  int begin_magic;

  wn_memgp next, *plast;
  wn_memgp children, sub_next, *psub_last;

  wn_mem_block current_block;
  ptr block_ptr, block_end_ptr;	/* used for allocating memory.  In no_free
  **				** groups, 8-byte aligned allocation is
  **				** from block_ptr, which goes from start of
  **				** block upward, other memory allocation is
  **				** from block_end_ptr, which goes from top
  **				** of block downward.  In general_free groups,
  **				** all allocation is through block_ptr. */
  int block_mem_left;
  int block_size;
  long int mem_used;

  ptr pad_list, free_list;

  void (*pverify_group)(wn_memgp group),
       (*pfree_piece)(ptr p,wn_memgp group);
  ptr  (*palloc_piece)(int size,wn_memgp group);

  char *label;

  int end_magic;
};

struct wn_mem_block_struct
{
  wn_mem_block next;
  int size, leftover;
};
#define WN_MEM_BLOCK_MEMORY(mem_block) ((ptr) ((mem_block)+1))


#define GROUP_STACK_SIZE    100

typedef struct wn_gpstack_struct *wn_gpstack;
struct wn_gpstack_struct      /* kludge for now */
{
  int stack_ptr;
  wn_memgp *group_stack;
  wn_memgp current_group;
};


EXTERN void wn_print_gp_stack(void);

EXTERN wn_memgp wn_curgp(void);
EXTERN wn_memgp wn_defaultgp(void);
EXTERN void wn_gppush(wn_memgp group);
EXTERN void wn_gppop(void);
EXTERN void wn_gpmake(const char parms[]);
EXTERN void wn_gpmakef(int flags, int block_size);
EXTERN void wn_gpfree(void);

EXTERN ptr wn_zalloc(int size);
EXTERN ptr wn_alloc(int size);
EXTERN void wn_free(ptr p);
EXTERN void wn_realloc(ptr *pp,int old_size,int new_size);
EXTERN void wn_zrealloc(ptr *pp,int old_size,int new_size);

EXTERN void wn_get_current_gpstack(wn_gpstack *pstack);
EXTERN void wn_set_current_gpstack(wn_gpstack stack);

EXTERN ptr wn_system_alloc(int size);
EXTERN void wn_system_free(ptr mem,int size);

EXTERN void wn_gplabel(const char *label);
#define WN_GET_GPLABEL(group) ((group)->label ? (group)->label : "")
EXTERN void wn_gperrfpush(void (*pfunc)(int size));
EXTERN void wn_stack_fill(void);
EXTERN void wn_gp_trap_address(ptr address);
EXTERN void wn_allmem_verify(void);
EXTERN long int wn_mem_used(void);
EXTERN void wn_print_mem_used(void);
EXTERN long int wn_group_mem_used(wn_memgp group);
EXTERN long unsigned wn_total_memory_in_group(wn_memgp group);
EXTERN long unsigned wn_total_wn_memory(void);
#if 0
  /*     dependent on sbrk, which no longer gives accurate heap size and
  ** isn't supported on Windows.  This routine was eliminated. */
  EXTERN long unsigned wn_heapsize(void);
#endif
EXTERN void wn_print_group_mem_used(wn_memgp group);
EXTERN void wn_print_all_groups_mem_used(void);
EXTERN bool wn_mem_in_group(ptr p,wn_memgp group);

/*   note, if you call either of these 2 routines, do so before commencing any
** memory allocation. */
EXTERN void wn_gp_fill(void);
EXTERN void wn_gp_pad(void);

EXTERN void wn_memory_fragmentation_report(void);
EXTERN long int wn_amount_of_free_memory_in_group(wn_memgp group);
EXTERN void wn_print_composition_of_big_blocks_free_list(wn_memgp group);

typedef void *(*wn_system_memory_alloc_func_type)(size_t size);
typedef void  (*wn_system_memory_free_func_type)(ptr p);

/* note the routines passed to these must be of type extern "C" */
EXTERN void wn_set_system_memory_alloc_func(wn_system_memory_alloc_func_type);
EXTERN void wn_set_system_memory_free_func(wn_system_memory_free_func_type);

/* note this is only to be called at the end of the program, after
** all memory allocation is concluded */
EXTERN void wn_abandon_memory_leaks(bool abandon_default_group);

/*     way to do your gpmakes, gpmakefs, & gpfrees such that you
**  check that there haven't been any unbalanced pushes or pops
**  in between. */

#define WN_GPBEGIN(parms_string)			\
  {							\
    wn_memgp _wn_begingp_memgp;				\
							\
    wn_gpmake(parms_string);				\
      _wn_begingp_memgp = wn_curgp();			\
      {

#define WN_GPBEGINF(flags, block_size)			\
  {							\
    wn_memgp _wn_begingp_memgp;				\
							\
    wn_gpmakef(flags, block_size);			\
      _wn_begingp_memgp = wn_curgp();			\
      {

#define WN_GPEND()					\
        ; /* in case we are preceeded by label */	\
      }							\
      wn_fast_assert(wn_curgp() == _wn_begingp_memgp);	\
    wn_gpfree();					\
  }

WN_EXTERN_END

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
#ifndef wnvectH
#define wnvectH




EXTERN void wn_print_vect(double vect[],int len);
EXTERN void wn_enter_vect(double vect[],int len);
EXTERN void wn_random_vect(double vect[],int len);
EXTERN void wn_zero_vect(double vect[],int len);

EXTERN void wn_make_vect(double **pvect,int len);
EXTERN void wn_free_vect(double *vect,int len);

EXTERN void wn_copy_vect(double to_vect[],double from_vect[],int len);

EXTERN double wn_dot_vects(double *vect1,double *vect2,int len);
EXTERN double wn_norm2_vect(double *vect,int len);
EXTERN double wn_norm_vect(double *vect,int len);
EXTERN void wn_unit_vect(double vect[],int len);
EXTERN double wn_dist2_vect(double v1[],double v2[],int len);
EXTERN double wn_dist_vect(double v1[],double v2[],int len);

EXTERN 
void wn_add_scaled_vect(double *to_vect,double *from_vect,double scalar,
			int len);
EXTERN
void wn_add_vect_and_scaled_vect(double *out,
				 double *v1,double *v2,double s2,
			         int len);
EXTERN void wn_add_scaled_vects
(
  double *out,
  double *v1,
  double s1,
  double *v2,
  double s2,
  int len
);
EXTERN void wn_add_3_scaled_vects
(
  double *out,
  double *v1,
  double s1,
  double *v2,
  double s2,
  double *v3,
  double s3,
  int len
);
EXTERN void wn_subtract_vect
(
  double *to_vect,
  double *from_vect,
  int len
);
EXTERN void wn_scale_vect(double *vect,double scalar,int len);

EXTERN void wn_multiply_vect_by_vect(double *v1, double *v2, int len);
EXTERN void wn_add_scaled_vects_and_dot_and_norm
(
  double *out,
  double *pnorm2_out,
  double *pdot_v1_out,
  double *v1,
  double s1,
  double *v2,
  double s2,
  int len
);

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
#ifndef wnconjH
#define wnconjH



typedef enum {WN_CONJ_ILLEGAL = 0,
	      WN_CONJ_DISTANCE, 
	      WN_CONJ_CURVATURE, 
	      WN_CONJ_MUL_FACTOR,
              WN_CONJ_NUM_ENTRIES} wn_jump_len_type;

EXTERN double wn_penalty(double x);
EXTERN double wn_dpenalty(double x);

EXTERN double wn_barrier(double x);
EXTERN double wn_dbarrier(double x);

EXTERN double wn_barrier2(double x);
EXTERN double wn_dbarrier2(double x);

EXTERN double wn_clamp(double x);
EXTERN double wn_dclamp(double x);

EXTERN void wn_conj_gradient_method
(
  int *pcode,
  double *pval_min,
  double vect[],
  int len,
  double (*pfunction)(double vect[]),
  void (*pgradient)(double grad[],double vect[]),
  int max_iterations
);
EXTERN void wn_conj_gradient_diff_method
(
  int *pcode,
  double *pval_min,
  double vect[],
  double delta_vect[],
  int len,
  double (*pfunction)(double vect[]),
  int max_iterations
);
EXTERN void wn_conj_funcgrad_method
(
  int *pcode,
  double *pval_min,
  double vect[],
  int len,
  double (*pfuncgrad)(double grad[],double vect[]),
  bool (*pterminate)(double vect[],double best_f,int iteration),
  double first_jump_len,
  wn_jump_len_type jump_len_flag,
  int no_reset_run_len,
  double reset_expand_factor,
  int max_iterations
);
EXTERN void wn_conj_funcgrad_method_parallel
(
  int *pcode,
  double *pval_min,
  double vect[],
  int num_vars,
  double (*pfuncgrad)(double grad[],double vect[]),
  bool (*pterminate)(double vect[],double best_f,int iteration),
  double first_step_len,
  wn_jump_len_type jump_len_flag,
  int no_reset_run_len,
  double reset_expand_factor,
  int max_iterations,
  wn_parallel_callback par_cb,
  int num_threads
);
EXTERN void wn_lbfgs_funcgrad_method
(
  int *pcode,
  double *pval_min,
  double vect[],
  int num_vars,
  double (*pfuncgrad)(double grad[],double vect[]),
  bool (*pterminate)(double vect[],double best_f),
  double first_jump_len,
  int max_depth, 
  int max_iterations
);

EXTERN void wn_force_conj_gradient_stop(void);

EXTERN void wn_numerical_gradient
(
  double grad[],
  double vect[],
  double delta_vect[],
  int len,
  double (*pfunction)(double vect[])
);


EXTERN void wn_fit_parabola_3p
(
  int *pcode,
  double *pa,double *px0,double *pb,
  double x1,double y1,
  double x2,double y2,
  double x3,double y3
);
EXTERN void wn_fit_parabola_2pd
(
  int *pcode,
  double *pa,double *px0,double *pb,
  double x1,double y1,double dy1,
  double x2,double y2
);
EXTERN void wn_fit_parabola_2d
(
  int *pcode,
  double *pa,double *px0,
  double x1,double dy1,
  double x2,double dy2
);
EXTERN void wn_fit_traditional_parabola_3p
(
  int *pcode,
  double *pa,double *pb, double *pc,
  double x1,double y1,
  double x2,double y2,
  double x3,double y3
);
EXTERN void wn_fit_traditional_parabola_2pa
(
  int *pcode,
  double *pb, double *pc,
  double a,
  double x1,double y1,
  double x2,double y2
);
EXTERN void wn_convert_parabola_traditional_to_centered
(
  int *pcode,
  double *px0,double *pb,
  double a,
  double b,
  double c
);
EXTERN void wn_fit_cubic_2p2d
(
  int *pcode,
  double coef[4],
  double x1,double y1,double dy1,
  double x2,double y2,double dy2
);

EXTERN double wn_clip_f(double);

EXTERN void wn_minimize_1d
(
  int *pcode,
  double *pval_min,
  double *px_min,
  double (*pfunction)(double x),
  int max_iterations
);
EXTERN void wn_minimize_1d_raw
(
  int *pcode,
  double *pf0,
  double *pf1,
  double *pf2,
  double *px0,
  double *px1,
  double *px2,
  double f_goal,
  double (*pfunction)(double x),
  int max_improve_iterations,
  int max_total_iterations
);


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

/****************************************************************************

  "cdn" stands for conjugate-direction noisy

****************************************************************************/
#ifndef wncnjdn2H
#define wncnjdn2H



#define WN_CDN_NOT_STARTED                -1



typedef struct wn_cdn_context_type_struct *wn_cdn_context_type;
typedef struct wn_cdn_srchdir_type_struct *wn_cdn_srchdir_type;

struct wn_cdn_context_type_struct
{
  int num_vars;

  double (*pfunction)(double vect[],int sequence_number);
  void (*pline_eval)(wn_cdn_context_type c,
                     double f[],double x[],int n,
                     double *search_direction);

  wn_cdn_srchdir_type *coord_search_direction_array,
                        *search_direction_array;
  int num_search_directions,max_num_search_directions;

  double *xa,*fa;
  int num_line_samples,min_num_line_samples,
      num_line_samples_mem;

  double *current_vect,*func_call_vect,
         *old_vect,*coord_direction;

  double ob;
  int code;

  int num_func_calls,max_num_func_calls;

  bool force_optimize_stop_flag /* = FALSE */;

  wn_memgp current_group;
};

struct wn_cdn_srchdir_type_struct
{
  bool is_coord_direction;
  int coord;
  double *dir_vect;

  double x_min,x_max;

  double x0;  /* optimimum x from last iteration */
  double curvature;  /* curvature from last iteration */
  double x_width,   /* x width from last iteration */
         max_x_width;
  
  double ratio_df1_noise;
};

EXTERN void wn_compute_func_call_vect
(
  wn_cdn_context_type c,
  double *search_direction,
  double x
);

EXTERN void wn_cdn_make_context
(
  wn_cdn_context_type *pc,
  double (*pfunction)(double vect[],int sequence_number),
  int num_vars
);

EXTERN void wn_cdn_free_context
(
  wn_cdn_context_type c
);

EXTERN void wn_cdn_set_solution
(
  wn_cdn_context_type c,
  double vect[]
);

EXTERN void wn_cdn_get_solution
(
  wn_cdn_context_type c,
  int *pcode,
  double *pob,
  double vect[]
);

EXTERN void wn_cdn_set_line_eval
(
  wn_cdn_context_type c,
  void (*pline_eval)(wn_cdn_context_type c,
                     double f[],double x[],int n,
                     double *search_direction)
);

EXTERN void wn_cdn_set_coord_x_min
(
  wn_cdn_context_type c,
  int coord,
  double x_min
);

EXTERN void wn_cdn_set_coord_x_max
(
  wn_cdn_context_type c,
  int coord,
  double x_max
);

EXTERN void wn_cdn_set_coord_x0
(
  wn_cdn_context_type c,
  int coord,
  double x0
);


EXTERN void wn_cdn_set_coord_x_width
(
  wn_cdn_context_type c,
  int coord,
  double x_width
);


EXTERN void wn_cdn_set_coord_max_x_width
(
  wn_cdn_context_type c,
  int coord,
  double max_x_width
);


EXTERN void wn_cdn_force_stop(wn_cdn_context_type c);


EXTERN void wn_cdn_optimize
(
  wn_cdn_context_type c,
  int num_func_calls
);



#endif

