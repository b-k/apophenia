
/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

****************************************************************************/

/********************************************************************

void conj_gradient_method(pcode,pval_min,vect,
			  len,pfunction,pgradient,max_iterations)

********************************************************************/
#include "wnlib.h"



#define GOLDEN_RATIO          1.618034
#define GOLDEN_SECTION        0.3819660

/* #define INFINITY  WN_FHUGE ** not used, is redefinition sometimes */

double wn_clip_f(double f)
{
  if(!(f < GSL_POSINF))
  {
    return(GSL_POSINF);
  }
  else
  {
    return(f);
  }
}

double wn_penalty(double x)
{
  if(x >= 0.0) 	return(0.0);
  else 		return(x*x);
}


double wn_dpenalty(double x)
{
  if(x >= 0.0) 	return(0.0);
  else 		return(2.0*x);
}


bool wn_force_optimize_stop_flag;

local wn_memgp old_group;

double *buffer_vect;

local double dy1,last_dy1;



void wn_force_conj_gradient_stop(void)
{
  wn_force_optimize_stop_flag = TRUE;
}

double wn_random_bit(gsl_rng *r){
double draw	= gsl_rng_uniform(r);
	if (draw >0.5)	return 1;
	else		return 0;
}


double wn_eval_poly(double x,double coefs[],int len)
{
  double ret;
  int i;

  if(len >= 2)
  {
    ret = coefs[len-2] + x*coefs[len-1];

    for(i=len-3;i>=0;--i)
    {
      ret = coefs[i] + x*ret;
    }
  }
  else
  {
    if(len <= 0)
    {
      ret = 0.0;
    }
    else  /* len == 1 */
    {
      ret = coefs[0];
    }
  }

  return(ret);
}


local bool parabola_fit_improvement_wrong(double y1,double y0,double b,
					  double tolerance)
{
  double expected_improvement;

  expected_improvement = y1-b;
  if(!(expected_improvement > 0.0))  /* weird form to handle Nan problems */
  {
    return(TRUE);
  }

  return(
	  /* not enough improvement */
	  (y0 > y1-(1.0-tolerance)*expected_improvement) 
	    ||
	  /* too much improvement */
	  (y0 < y1-(1.0+tolerance)*expected_improvement) 
        );
}


local void eval_function
(
  double *pval,
  double vect[],
  double direction[],
  double x,
  int len,
  double (*pfunction)(double vect[])
)
{
  wn_add_vect_and_scaled_vect(buffer_vect,vect,direction,x,len);

  *pval = (*pfunction)(buffer_vect);
}


local double *save_vect,*save_direction;
local int save_len;
local double (*save_pfunction)(double vect[]);

local double simple_eval_function(double x)
{
  double ret;

  eval_function(&ret,save_vect,save_direction,x,save_len,save_pfunction);

  return(ret);
}


local bool x0_too_far_out(double x0,double x2,double threshold)
{
  wn_assert(threshold > 1.0);

  if(x2 > 0.0)
  {
    return(x0 > x2*threshold);
  }
  else
  {
    wn_assert(x2 < 0.0);  /* x2 == 0 not allowed */

    return(x0_too_far_out(-x0,-x2,threshold));
  }
}


local bool x0_too_far_in(double x0,double x2,double threshold)
{
  wn_assert(threshold < 1.0);

  if(x2 > 0.0)
  {
    return(x0 < x2*threshold);
  }
  else
  {
    wn_assert(x2 < 0.0);  /* x0 == 0 not allowed */

    return(x0_too_far_in(-x0,-x2,threshold));
  }
}


local void line_minimize_a
(
  int *pcode,
  bool *pfirst_parabolic_fit_succeeded,
  double *pval,
  double *pjump_len,
  double vect[],
  double direction[],
  double last_val,
  double last_g[],
  double last_jump_len,
  int len,
  double (*pfunction)(double vect[])
)
{
  double a,x0,y0,b,x1,y1,x2,y2;
  int code;

  *pcode = WN_SUCCESS;

  x1 = 0.0;
  y1 = last_val;
  dy1 = wn_dot_vects(direction,last_g,len);

  if(last_jump_len == 0.0)
  {
    last_jump_len = 1.0;
  }
  last_jump_len = -wn_sign(dy1)*wn_abs(last_jump_len);
  x2 = last_jump_len;

  eval_function(&y2,vect,direction,x2,len,pfunction);

  wn_fit_parabola_2pd(&code,&a,&x0,&b,x1,y1,dy1,x2,y2);

  /* look for excuses to say that parabolic fit is no good */
  if(code != WN_SUCCESS)
  {
    /*
    printf("parabola fit failed - probably a line.\n");
    */
    x0 = x2+GOLDEN_RATIO*x2;  /* project outward */
    eval_function(&y0,vect,direction,x0,len,pfunction);
    goto parabolic_fit_failed;
  }
  if(!(a > 0))
  {
    /*
    printf("downward facing parabola.\n");
    */
    x0 = x2+GOLDEN_RATIO*x2;  /* project outward */
    eval_function(&y0,vect,direction,x0,len,pfunction);
    goto parabolic_fit_failed;
  }
  if(!(wn_abs(x0) < 10000.0*wn_abs(x2)))
  {
    /*
    printf("x0 too far out.\n");
    */
    x0 = 10000.0*x2;  /* project outward */
    eval_function(&y0,vect,direction,x0,len,pfunction);
    goto parabolic_fit_failed;
  }
  if(!(wn_abs(x0) > (1.0/10000.0)*wn_abs(x2)))
  {
    /*
    printf("x0 too far in.\n");
    */
    x0 = (1.0/10000.0)*x2;  /* project inward */
    eval_function(&y0,vect,direction,x0,len,pfunction);
    goto parabolic_fit_failed;
  }
  if(!(b < y1))  /* no improvement expected,weird form for Nan problems */
  {
    /*
    printf("no improvement expected.\n");
    */
    x0 = GOLDEN_SECTION*x2;
    eval_function(&y0,vect,direction,x0,len,pfunction);
    goto parabolic_fit_failed;
  }

  eval_function(&y0,vect,direction,x0,len,pfunction);

  if(parabola_fit_improvement_wrong(y1,y0,b,0.25))
  {
    /*
    printf("poor parabola fit detected.\n");
    */
    goto parabolic_fit_failed;
  }

  /* parabolic fit succeeded */

  if(y0 > y1)
  {
    x0 = x1;
    y0 = y1;
  }
  *pval = y0;
  *pjump_len = x0;
  wn_copy_vect(vect,buffer_vect,len);

  *pfirst_parabolic_fit_succeeded = TRUE;
  /*
  *pfirst_parabolic_fit_succeeded = 
	 !parabola_fit_improvement_wrong(y1,y0,b,0.25);
  */

  return;

parabolic_fit_failed:
  *pfirst_parabolic_fit_succeeded = FALSE;

  save_vect = vect;
  save_direction = direction;
  save_len = len;
  save_pfunction = pfunction;

  wn_minimize_1d_raw(pcode,&y1,&y0,&y2,&x1,&x0,&x2,y1,(simple_eval_function),
		     3,20);
  if(!((*pcode == WN_SUCCESS)||(*pcode == WN_SUBOPTIMAL)))
  {
    return;
  }
  *pcode = WN_SUCCESS;

  if(y0 <= last_val)
  {
    *pval = y0;
    *pjump_len = x0;
    wn_add_scaled_vect(vect,direction,x0,len);
  }
  else
  {
    *pval = last_val;
    *pjump_len = 0.0;
  }

  return;
}


EXTERN void wn_conj_gradient_method
(
  int *pcode,
  double *pval_min,
  double vect[],
  int len,
  double (*pfunction)(double vect[]),
  void (*pgradient)(double grad[],double vect[]),
  int max_iterations
)
{
  int iteration,no_move_count;
  int stable_satisfy_count;
  double norm2_g,norm2_last_g,g_dot_last_g,val,last_val,beta,
	 jump_len,last_jump_len,alpha;
  double *g,*last_g,*direction;
  bool function_free_method,last_was_function_free_method;
  bool first_parabolic_fit_succeeded;
  double computed_g_dot_dlast;

  old_group = wn_curgp();
  wn_gpmake("no_free");

  wn_force_optimize_stop_flag = FALSE;

  wn_make_vect(&buffer_vect,len);
  wn_make_vect(&g,len);
  wn_make_vect(&last_g,len);
  wn_make_vect(&direction,len);

  function_free_method = FALSE;
  stable_satisfy_count = 0;
  last_was_function_free_method = FALSE;
  jump_len = 1.0;
  no_move_count = 0;
  beta = 0.0;

  wn_gppush(old_group);

  val = (*pfunction)(vect);
  (*pgradient)(g,vect);

  wn_copy_vect(direction,g,len);

  norm2_g = wn_norm2_vect(g,len);

  for(iteration=0;;++iteration)
  {
    last_dy1 = dy1;
    last_jump_len = jump_len;
    wn_swap(last_g,g,double *);  /* move g to last g */
    norm2_last_g = norm2_g;

    if(function_free_method)
    {
      double x2;
      double g0_dot_d0,g1s_dot_d0,dot_diff;

      dy1 = wn_dot_vects(direction,last_g,len);

      last_jump_len = -wn_sign(dy1)*wn_abs(last_jump_len);
      x2 = last_jump_len;

      /*
      printf("last_jump_len = %lg\n",last_jump_len);
      */

      wn_add_vect_and_scaled_vect(buffer_vect,vect,direction,last_jump_len,len);

      (*pgradient)(g,buffer_vect);

      /*
      g0_dot_d0 = wn_dot_vects(last_g,direction,len);
      */
      g0_dot_d0 = dy1;
      g1s_dot_d0 = wn_dot_vects(g,direction,len);
      dot_diff = g0_dot_d0-g1s_dot_d0;
     
      if(!(dot_diff > 0.0))  /* not upward facing parabola */
      {
        stable_satisfy_count = 0;
	goto function_based_method;
      }

      alpha = g0_dot_d0/dot_diff;
      /*
      printf("alpha = %lg\n",alpha);
      */

      if(!(((1.0-10000.0) < alpha)&&(alpha < (1.0+10000.0))))
      {
        stable_satisfy_count = 0;
	goto function_based_method;
      }

      jump_len = alpha*last_jump_len;

      /* g[j] = alpha*g[j] + (1.0-alpha)*last_g[j]; */
      wn_scale_vect(g,alpha,len);
      wn_add_scaled_vect(g,last_g,1.0-alpha,len);

      g_dot_last_g = wn_dot_vects(g,last_g,len);
      if(beta != 0.0)
      {
	computed_g_dot_dlast = -g_dot_last_g/beta;
	/*
	printf("computed_g_dot_dlast=%lg,last_dy1=%lg\n",
	       computed_g_dot_dlast,last_dy1);
	*/
	if(!(wn_abs(computed_g_dot_dlast) < 0.4*wn_abs(last_dy1)))
	{
          stable_satisfy_count = 0;
	  goto function_based_method;
        }
      }

      wn_add_scaled_vect(vect,direction,alpha*last_jump_len,len);

      if(
          wn_force_optimize_stop_flag
	    ||
	  ((max_iterations < WN_IHUGE)&&(iteration >= max_iterations))
        )
      {
        wn_force_optimize_stop_flag = FALSE;
        wn_gppop();
        wn_gpfree();
        val = (*pfunction)(vect);
        *pval_min = val;
        *pcode = WN_SUBOPTIMAL;
        return;
      }
    }
    else  /* function based method */
    {
      function_based_method: 

      if(last_was_function_free_method)
      {
        /* set so that next iteration will succeed */
        val = (*pfunction)(vect);  
      }
      function_free_method = FALSE;

      if(norm2_last_g == 0.0)   /* unlikely */
      {
        wn_gppop();
        wn_gpfree();
        *pval_min = val;
        *pcode = WN_SUCCESS;
        return;
      }

      last_val = val;

      line_minimize_a(pcode,&first_parabolic_fit_succeeded,&val,&jump_len,
		    vect,
		    direction,last_val,last_g,last_jump_len,len,pfunction);
      if(*pcode != WN_SUCCESS)
      {
        wn_gppop();
        wn_gpfree();
        *pval_min = val;
        return;
      }

      if(
	  wn_force_optimize_stop_flag
	    ||
	  ((max_iterations < WN_IHUGE)&&(iteration >= max_iterations))
        )
      {
        wn_force_optimize_stop_flag = FALSE;
        wn_gppop();
        wn_gpfree();
        *pval_min = val;
        *pcode = WN_SUBOPTIMAL;
        return;
      }
      wn_assert(val <= last_val);
      if(val == last_val)
      {
        if(no_move_count >= 2)
        {
          wn_gppop();
          wn_gpfree();
          *pval_min = val;
          *pcode = WN_SUCCESS;
          return;
        }
        else
        {
	  ++no_move_count;
	  jump_len = last_jump_len;
        }
      }
      else
      {
        no_move_count = 0;
      }

      (*pgradient)(g,vect);
      g_dot_last_g = wn_dot_vects(g,last_g,len);

      if((!first_parabolic_fit_succeeded)||(last_jump_len == 0.0))
      {
        stable_satisfy_count = 0;
	goto no_function_method_test_fail;
      }

      alpha = jump_len/last_jump_len;

      if(!(((1.0-3000.0) < alpha)&&(alpha < (1.0+3000.0))))
      {
        stable_satisfy_count = 0;
	goto no_function_method_test_fail;
      }

      if(beta != 0.0)
      {
	computed_g_dot_dlast = -g_dot_last_g/beta;
	/*
	printf("computed_g_dot_dlast=%lg,last_dy1=%lg\n",
	       computed_g_dot_dlast,last_dy1);
	*/
	if(!(wn_abs(computed_g_dot_dlast) < 0.2*wn_abs(last_dy1)))
	{
          stable_satisfy_count = 0;
	  goto no_function_method_test_fail;
        }
      }

      ++stable_satisfy_count;
      if(stable_satisfy_count > 3)
      {
        function_free_method = TRUE;
      }

      no_function_method_test_fail: ;
    }

    norm2_g = wn_norm2_vect(g,len);
    /*
    g_dot_last_g = wn_dot_vects(g,last_g,len);
    */

    beta = (norm2_g - g_dot_last_g)/norm2_last_g;

    wn_add_vect_and_scaled_vect(direction,g,direction,beta,len);

    /*
    printf("norm(g) = %lg\n",wn_norm_vect(g,len));
    printf("ob = %lg,beta = %lg,numerator=%lg,denom=%lg,norm2(direction)=%lg\n",
	   val,beta,numerator,norm2_last_g,wn_norm2(direction));
    printf("iteration = %d,ob = %lg\n",iteration,val);
    */

    last_was_function_free_method = function_free_method;
  }
}





/********************************************************************
Fit parabola to 3 points (x,y).  Parabola is of the form

  y = a*(x-x0)^2 + b

Return a, x0, b.
********************************************************************/
void wn_fit_parabola_3p
(
  int *pcode,
  double *pa,double *px0,double *pb,
  double x1,double y1,
  double x2,double y2,
  double x3,double y3
)
{
  double x12,x23,dx21,dx32,dx2312,dy12,dy23,ddy,diff;

  /* sort the x's */
  if(!(x1 < x3))
  {
    wn_swap(x1,x3,double);
    wn_swap(y1,y3,double);
  }
  if(x2 < x1)
  {
    wn_swap(x1,x2,double);
    wn_swap(y1,y2,double);
  }
  else if(x3 < x2) 
  {
    wn_swap(x2,x3,double);
    wn_swap(y2,y3,double);
  }

  dx21 = x2-x1;
  dx32 = x3-x2;

  if(!((x1 != x3)&&(dx21 != 0.0)&&(dx32 != 0.0)))
  {
    *pcode = WN_SINGULAR;
    return;
  }

  x12 = 0.5*(x1+x2);
  x23 = 0.5*(x2+x3);

  dx2312 = x23-x12;
  if(!(dx2312 != 0.0))
  {
    *pcode = WN_SINGULAR;
    return;
  }

  dy12 = (y2-y1)/dx21;
  dy23 = (y3-y2)/dx32;

  ddy = dy23-dy12;
  *pa = 0.5*ddy/dx2312;

  if(ddy != 0.0)
  {
    *px0 = (dy23*x12-dy12*x23)/ddy;

    diff = x2-(*px0);
    *pb = y2-(*pa)*diff*diff;
  }
  else
  {
    *px0 = 0.0;
    *pb = y2;
  }

  *pcode = WN_SUCCESS;
}


/********************************************************************
Fit parabola to 2 points (x,y) and a derivative at one point.  
Parabola is of the form

  y = a(x-x0)^2 + b

Return a, x0, b.
********************************************************************/
void wn_fit_parabola_2pd
(
  int *pcode,
  double *pa,double *px0,double *pb,
  double x1,double y1,double dy1,
  double x2,double y2
)
{
  double x12,dx21,dx121,dy12,ddy,diff;

  dx21 = x2-x1;

  if(!(dx21 != 0.0))
  {
    *pcode = WN_SINGULAR;
    return;
  }

  x12 = 0.5*(x1+x2);

  dy12 = (y2-y1)/dx21;

  ddy = dy12-dy1;

  dx121 = x12-x1;
  if(!(dx121 != 0.0))
  {
    *pcode = WN_SINGULAR;
    return;
  }
  *pa = 0.5*ddy/dx121;

  if(ddy != 0.0)
  {
    *px0 = (dy12*x1-dy1*x12)/ddy;

    diff = x1-(*px0);
    *pb = y1-(*pa)*diff*diff;
  }
  else
  {
    *px0 = 0.0;
    *pb = y1;
  }

  *pcode = WN_SUCCESS;
}


/********************************************************************
Fit parabola to derivatives at 2 points.
Parabola is of the form

  y = a(x-x0)^2 + b

Return a, x0.  b is impossible to compute so it is not returned.
********************************************************************/
void wn_fit_parabola_2d
(
  int *pcode,
  double *pa,double *px0,
  double x1,double dy1,
  double x2,double dy2
)
{
  double dx,ddy;

  dx = x2-x1;

  if(!(dx != 0.0))
  {
    *pcode = WN_SINGULAR;
    return;
  }

  ddy = dy2-dy1;

  if(!(ddy != 0.0))
  {
    *pcode = WN_SINGULAR;
    return;
  }

  *pa = 0.5*ddy/dx;
  *px0 = (dy2*x1-dy1*x2)/ddy;

  *pcode = WN_SUCCESS;
}


/********************************************************************
Fit parabola to 3 points (x,y).  Parabola is of the form

  y = a*x^2 + b*x + c

Return a, b, c.
********************************************************************/
void wn_fit_traditional_parabola_3p
(
  int *pcode,
  double *pa,double *pb, double *pc,
  double x1,double y1,
  double x2,double y2,
  double x3,double y3
)
{
  double x12,x23,x13,dx21,dx32,dx31,dx2312,dy21,dy32,dy31,ddy;

  /* sort the x's */
  if(!(x1 < x3))
  {
    wn_swap(x1,x3,double);
    wn_swap(y1,y3,double);
  }
  if(x2 < x1)
  {
    wn_swap(x1,x2,double);
    wn_swap(y1,y2,double);
  }
  else if(x3 < x2) 
  {
    wn_swap(x2,x3,double);
    wn_swap(y2,y3,double);
  }

  dx21 = x2-x1;
  dx32 = x3-x2;
  dx31 = x3-x1;

  if(!((x1 != x3)&&(dx21 != 0.0)&&(dx32 != 0.0)&&(dx31 != 0)))
  {
    *pcode = WN_SINGULAR;
    return;
  }

  x12 = 0.5*(x1+x2);
  x23 = 0.5*(x2+x3);
  x13 = 0.5*(x1+x3);

  dx2312 = x23-x12;
  if(!(dx2312 != 0.0))
  {
    *pcode = WN_SINGULAR;
    return;
  }

  dy21 = (y2-y1)/dx21;
  dy32 = (y3-y2)/dx32;
  dy31 = (y3-y1)/dx31;

  ddy = dy32-dy21;
  *pa = 0.5*ddy/dx2312;

  *pb = dy31 - 2.0*(*pa)*x13;

  *pc = y2 - (*pa)*x2*x2 - (*pb)*x2;

  *pcode = WN_SUCCESS;
}


/********************************************************************
Fit parabola to 2 points (x,y), given that curvature a is already
known.  Parabola is of the form

  y = a*x^2 + b*x + c

Return b, c.
********************************************************************/
void wn_fit_traditional_parabola_2pa
(
  int *pcode,
  double *pb, double *pc,
  double a,
  double x1,double y1,
  double x2,double y2
)
{
  double dx,sx,mult_dxsx;

  /* sort the x's */
  if(x2 < x1)
  {
    wn_swap(x1,x2,double);
    wn_swap(y1,y2,double);
  }

  dx = x2-x1;

  if(dx == 0.0)
  {
    *pcode = WN_SINGULAR;
    return;
  }

  sx = x2+x1;
  mult_dxsx = dx*sx;

  *pb = ((y2-y1) - a*mult_dxsx)/dx;
  *pc = 0.5*((y2+y1) - a*(x2*x2+x1*x1) - (*pb)*sx);

  *pcode = WN_SUCCESS;
}


/********************************************************************
Convert traditional parabola of the form

  y = a*x^2 + b*x + c

to centered parabola of form 

  y = a*(x-x0)^2 + b

Return x0, b.
********************************************************************/
void wn_convert_parabola_traditional_to_centered
(
  int *pcode,
  double *px0,double *pb,
  double a,
  double b,
  double c
)
{
  if(a == 0.0)
  {
    *pcode = WN_SINGULAR;
    return;
  }

  *px0 = -0.5*b/a;
  *pb = c - 0.25*b*b/a;

  *pcode = WN_SUCCESS;
}


/********************************************************************
Fit traditional cubic to  2 points (x,y) and derivatives
Cubic is of the form

  y = c3*x^3 + c2*x^2 + c1*x + c0

Return c3, c2, c1, c0
********************************************************************/
void wn_fit_cubic_2p2d
(
  int *pcode,
  double coef[4],
  double x1,double y1,double dy1,
  double x2,double y2,double dy2
)
{
  double c[4];
  double m,b;
  double xform[2];
  double accum[5],next_accum[5];
  int i,j;

  *pcode = WN_SUCCESS;

  if(x1 == x2)
  {
    *pcode = WN_SINGULAR;
    return;
  }

  /* compute mapping from x space to -1 +1 space */
  m = 2.0/(x2-x1);
  b = 1.0 - m*x2;

  dy1 /= m;
  dy2 /= m;

  /* compute coefs assuming x1==-1,x2==1 */
  c[2] = (1.0/4.0)*(dy2-dy1);
  c[0] = (1.0/2.0)*(y2+y1) - c[2];
  c[3] = (1.0/2.0)*(y1+dy1 + c[2] - c[0]);
  c[1] = dy2 - 3.0*c[3] - 2.0*c[2];

  /* use mapping to transform cubic */
  wn_zero_vect(coef,4);
  wn_zero_vect(accum,4);
  accum[0] = 1.0;

  xform[0] = b;
  xform[1] = m;

  for(i=0;i<4;++i)
  {
    for(j=0;j<=i;++j)
    {
      coef[j] += c[i]*accum[j];
    }

    wn_mult_polys(next_accum,accum,i+1,xform,2);
    wn_copy_vect(accum,next_accum,i+2);
  }

  *pcode = WN_SUCCESS;
}
/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

****************************************************************************/

/********************************************************************
  void wn_conj_direction_method
  (pcode, pval_min, vect, num_vars, pfunction, max_func_calls)
********************************************************************/


#define MAX_EXPAND   (100.0)
#define MIN_CONTRACT (1.0/10.0)

int wn_conj_direction_debug = WN_CONJ_DIR_DBG_NONE;

/* local bool show_linesearch=FALSE; ** unused - bchapman 041111 */

local int num_vars;

local double **search_directions;
local int num_search_directions,max_num_search_directions;

local double *coord_x0s,*search_direction_x0s;
local double *coord_as,*search_direction_as;

local int num_func_calls;

local bool last_line_function_x_valid;
local double last_line_function_x,last_line_function_ret;

local double *save_vect,*save_direction;
local double (*save_pfunction)(double vect[]);

local bool force_optimize_stop_flag=FALSE;

local double sqrt_tolerance;


local double fract_diff(double n1,double n2)
{
  n1 = wn_abs(n1);
  n2 = wn_abs(n2);

  if(n1 > n2)
  {
    return(1.0-n2/n1);
  }
  else if(n2 > n1)
  {
    return(1.0-n1/n2);
  }
  else if(n1 == n2)
  {
    return(0.0);
  }
  else
  {
    return(1.0);
  }
}


local bool is_valid_number(double x)
{
  return((-WN_FHUGE < x)&&(x < WN_FHUGE));
}


local bool too_close(double x1, double x2)
{
  if((x2 == 0.0) && (wn_abs(x1) < sqrt_tolerance))
  {
    return(TRUE);
  }
  else if((x1 == 0.0) && (wn_abs(x2) < sqrt_tolerance))
  {
    return(TRUE);
  }
  else if(
	  (wn_sign(x1) == wn_sign(x2))
	   &&
	  (fract_diff(x1, x2) < sqrt_tolerance)
	 )
  {
   return(TRUE);
  }
  else if(wn_abs(x1 - x2) < sqrt_tolerance)
  {
   return(TRUE);
  }

  return(FALSE);
}


void wn_force_conj_direction_stop(void)
{
  force_optimize_stop_flag = TRUE;
}


local double line_function(double x)
{
  double ret;

  if(last_line_function_x_valid)
  {
    if(x == last_line_function_x)
    {
      return(last_line_function_ret);
    }
  }

  wn_copy_vect(save_vect,buffer_vect,num_vars);
  wn_add_scaled_vect(save_vect,save_direction,x,num_vars);

  if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_ALL)
  {
    printf("function call %d at ",num_func_calls);
    wn_print_vect(save_vect,num_vars);
  }

  ++num_func_calls;

  ret = wn_clip_f((*save_pfunction)(save_vect));

  if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_ALL)
  {
    printf("function value is %lg\n",ret);
    fflush(stdout);
  }

  last_line_function_x_valid = TRUE;
  last_line_function_x = x;
  last_line_function_ret = ret;

  return(ret);
}


local void fit_parabola_2pa(int *pcode,double *px0,double *pb,
			    double a,
			    double x1,double y1,
			    double x2,double y2)
{
  double b,c;

  if(a <= 0.0)
  {
    *pcode = WN_SINGULAR;
    return;
  }

  wn_fit_traditional_parabola_2pa(pcode,&b,&c,a,x1,y1,x2,y2);

  if(*pcode != WN_SUCCESS)
  {
    return;
  }

  wn_convert_parabola_traditional_to_centered(pcode,px0,pb,a,b,c);
}


local void line_minimize_b
(
  double vect[],
  double direction[],
  double *pval_min,
  double *psave_x0,
  double *psave_a,
  double (*pfunction)(double vect[])
)
{
  double ax,bx,cx,x0,fa,fb,fc,fx0;
  double a,b;
  double old_x0,old_a;
  int code;
  gsl_rng *r	= gsl_rng_alloc(gsl_rng_taus);

  if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_LINESEARCH)
  {
    printf("start line minimize.\n");
  }

  last_line_function_x_valid = FALSE;

  wn_copy_vect(buffer_vect,vect,num_vars);
  save_vect = vect;
  save_direction = direction;
  save_pfunction = pfunction;

  old_x0 = *psave_x0;
  old_a = *psave_a;

  bx = 0.0;
  fb = *pval_min;
  if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_ALL)
  {
    printf("First point at %lg, function value = %lg\n", bx, fb);
  }

  if(old_x0 == 0.0)
  {
    old_x0 = 1.0;
  }

  ax = old_x0*apop_random_double(0.9,1.1,r);
  if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_ALL)
  {
    printf("Second point at %lg (old_x0 = %lg)\n", ax, old_x0);
  }
  fa = line_function(ax);

  if(!(old_a > 0.0))
  {
    goto simple_parabola_fit;
  }

  /* the curvature along a search direction is constant for a 
     quadratic function, therefore, try to use the curvature
     from the last search */
  fit_parabola_2pa(&code,&x0,&b,old_a,ax,fa,bx,fb);
  if(
     (code != WN_SUCCESS)
       ||
     (!(wn_abs(x0)<MAX_EXPAND*wn_abs(old_x0)) && (*psave_x0 != 0.0))
       ||
     too_close(x0, ax) || too_close(x0, bx)
       ||
     !is_valid_number(x0) || !is_valid_number(ax) || !is_valid_number(bx)
    )
  {
    goto simple_parabola_fit;
  }   

  cx = x0;
  if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_ALL)
  {
    printf("Third point at %lg\n", cx);
  }
  fc = line_function(cx);

  wn_fit_parabola_3p(&code,&a,&x0,&b,ax,fa,bx,fb,cx,fc);

  if((code != WN_SUCCESS)||(!(a > 0.0))||
     (!(wn_abs(x0)<MAX_EXPAND*wn_abs(old_x0))&&(*psave_x0 != 0.0)))
  {
    goto full_linesearch;
  }   

  if(!(b < fb))
  {
    if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_ALL)
    {
      printf("Doing slow line search (parabola fit returned suspect min value).\n");
    }
    goto full_linesearch;
  }
  if((!(fc < fb))||(!(fc < fa)))
  {
    /* evaluate one more point */
    goto evaluate_x0;
  }

  /* is it economical to evaluate one more point? */
  if((fb-b) <= 1.5*(fb-fc))
  { 
    /* do not evaluate one more point */
    wn_swap(fb,fc,double);
    wn_swap(bx,cx,double);
    goto finish;
  }
  else
  {
    /* evaluate one more point */
    goto evaluate_x0;
  }

simple_parabola_fit:
  if(fa < fb)
  {
    cx = 2.0*ax*apop_random_double(0.8,1.2,r);
  }
  else
  {
    cx = -1.0*ax*apop_random_double(0.8,1.2,r);
  }

  fc = line_function(cx);

  wn_fit_parabola_3p(&code,&a,&x0,&b,ax,fa,bx,fb,cx,fc);

  if((code != WN_SUCCESS)||(!(a > 0.0))||
     (!(wn_abs(x0)<MAX_EXPAND*wn_abs(old_x0))&&(*psave_x0 != 0.0)))
  {

    if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_ALL)
    {
      printf("Parabola fit failed. Switching to slow line search mode.\n");
    }
    goto full_linesearch;
  }   

evaluate_x0:
  fx0 = line_function(x0);

  if(!(fx0 <= fb))
  {
    if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_ALL)
    {
      printf("Doing a slow line search because f(x0) is too large (x0 = %lg).\n", x0);
    }
    goto full_linesearch;
  }

  fb = fx0;
  bx = x0;

  if(!(fa <= fc))
  {
    wn_swap(fa,fc,double);
    wn_swap(ax,cx,double);
  }
  if(!(fb <= fa))
  {
    wn_swap(fb,fa,double);
    wn_swap(bx,ax,double);
  }

  goto finish;

full_linesearch: ;

  /*
  printf("now.\n");
  */
  do
  {
    if(ax == bx)
    {
      if(wn_random_bit(r))
      {
        ax += apop_random_double(-1.0,1.0,r);
        fa = line_function(ax);
      }
      else
      {
        bx += apop_random_double(-1.0,1.0,r);
        fb = line_function(bx);
      }
    }
    if(ax == cx)
    {
      if(wn_random_bit(r))
      {
        ax += apop_random_double(-1.0,1.0,r);
        fa = line_function(ax);
      }
      else
      {
        cx += apop_random_double(-1.0,1.0,r);
        fc = line_function(cx);
      }
    }
    if(bx == cx)
    {
      if(wn_random_bit(r))
      {
        bx += apop_random_double(-1.0,1.0,r);
        fb = line_function(bx);
      }
      else
      {
        cx +=  apop_random_double(-1.0,1.0,r);
        fc = line_function(cx);
      }
    }
  } while((ax == bx)||(ax == cx)||(bx == cx));
  wn_minimize_1d_raw(&code,&fa,&fb,&fc,&ax,&bx,&cx,fb,&line_function,1,20);
  /*
  printf("l = %lf\n",bx);
  */

finish: ;

  /*
  if(show_linesearch)
  {
    printf("ax=%lg,bx=%lg,cx=%lg,old_x0=%lg\n",ax,bx,cx,old_x0);
  }
  */

  wn_copy_vect(vect,buffer_vect,num_vars);

  /* compute *psave_x0 */
  if(wn_abs(bx) < MIN_CONTRACT*wn_abs(old_x0))
  {
    if(bx < 0.0)
    {
      *psave_x0 = -MIN_CONTRACT*wn_abs(old_x0);
    }
    else
    {
      *psave_x0 = MIN_CONTRACT*wn_abs(old_x0);
    }
  }
  else
  {
    *psave_x0 = bx;
  }

  /* compute *psave_a */
  wn_fit_parabola_3p(&code,&a,&x0,&b,ax,fa,bx,fb,cx,fc);

  if((code != WN_SUCCESS)||(!(a > 0.0)))
  {
    *psave_a = 0.0;
  }   
  else
  {
    *psave_a = a;
  }

  if(*pval_min == fb)
  {
    if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_LINESEARCH)
    {
      printf("finish line minimize.\n");
      fflush(stdout);
    }
    return;  /* do not move if no improvement */
  }

  wn_add_scaled_vect(vect,direction,bx,num_vars);

  *pval_min = fb;

  if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_LINESEARCH)
  {
    printf("finish line minimize.\n");
    fflush(stdout);
  }
}


EXTERN void wn_conj_direction_method
(
  int *pcode,
  double *pval_min,
  double vect[],
  double initial_coord_x0s[],
  int passed_num_vars,
  double (*pfunction)(double vect[]),
  int max_func_calls
)
{
  int i,j,iteration;
  double *old_vect,*coord_direction;
  double *new_search_direction;
  double old_val_min;

  wn_memgp conj_dir_memgp;

  wn_gpmake("no_free");

  force_optimize_stop_flag = FALSE;

  num_vars = passed_num_vars;
  max_num_search_directions = num_vars;

  wn_make_vect(&buffer_vect,num_vars);
  search_directions = (double **)wn_zalloc(
		  max_num_search_directions*sizeof(double *));
  wn_make_vect(&old_vect,num_vars);
  wn_make_vect(&coord_direction,num_vars);
  wn_make_vect(&coord_x0s,num_vars);
  wn_make_vect(&search_direction_x0s,max_num_search_directions);
  wn_make_vect(&coord_as,num_vars);
  wn_make_vect(&search_direction_as,max_num_search_directions);
  if(initial_coord_x0s == NULL)
  {
    wn_zero_vect(coord_x0s,num_vars);
  }
  else
  {
    wn_copy_vect(coord_x0s,initial_coord_x0s,num_vars);
  }
  wn_zero_vect(search_direction_x0s,max_num_search_directions);
  wn_zero_vect(coord_as,num_vars);
  wn_zero_vect(search_direction_as,max_num_search_directions);

  /* name and pop the memory group we created */

  conj_dir_memgp = wn_curgp();
  wn_gppop();

  sqrt_tolerance = sqrt(GSL_DBL_EPSILON);

  num_search_directions = 0;

  num_func_calls = 0;

  last_line_function_x_valid = FALSE;
  *pval_min = wn_clip_f((*pfunction)(vect));
  ++num_func_calls;

  for(iteration=0;;++iteration)
  {
    if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_PASSES)
    {
      printf("iteration = %d ********************************\n",iteration);
      printf("ob = %lg\n",*pval_min);
      fflush(stdout);
    }
    /*
    (void)getchar();
    */

    old_val_min = *pval_min;
    wn_copy_vect(old_vect,vect,num_vars);

    /* minimize along acceleration search directions */
    if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_LINESEARCH)
    {
      printf("start acceleration line minimizations ------------------\n");
    }
    for(i=0;i<num_search_directions;++i)
    {
      if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_LINESEARCH)
      {
        printf("acceleration line search %d\n",i);
      }

      line_minimize_b(vect,search_directions[i],pval_min,
		    &(search_direction_x0s[i]),&(search_direction_as[i]),
		    pfunction);

      if(
	  ((max_func_calls < WN_IHUGE)&&(num_func_calls > max_func_calls))
	    ||
          force_optimize_stop_flag
	)
      {
        *pcode = WN_SUBOPTIMAL;
	goto finish;
      }
    }

    /* minimize along coordinate directions */
    if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_LINESEARCH) 
    {
      printf("start coordinate line minimizations ------------------\n");
    }
    for(i=0;i<num_vars;++i)
    {
      if(wn_conj_direction_debug >= WN_CONJ_DIR_DBG_LINESEARCH)
      {
        printf("coord line search %d\n",i);
      }

      coord_direction[i] = 1.0;

      line_minimize_b(vect,coord_direction,pval_min,
		    &(coord_x0s[i]),&(coord_as[i]),
		    pfunction);

      coord_direction[i] = 0.0;

      if(
	  ((max_func_calls < WN_IHUGE)&&(num_func_calls > max_func_calls))
	    ||
          force_optimize_stop_flag
	)
      {
        *pcode = WN_SUBOPTIMAL;
	goto finish;
      }
    }

    if(*pval_min >= old_val_min)
    {
      wn_assert(*pval_min == old_val_min);

      *pcode = WN_SUCCESS;
      break;
    }

    /* compute new acceleration search direction */
    if(num_search_directions < max_num_search_directions)
    {
      wn_gppush(conj_dir_memgp);
      wn_make_vect(&new_search_direction,num_vars);
      wn_gppop();
      for(i=num_search_directions;i>0;--i)
      {
        search_directions[i] = search_directions[i-1];
	search_direction_x0s[i] = search_direction_x0s[i-1];
	search_direction_as[i] = search_direction_as[i-1];
      }
      search_directions[0] = new_search_direction;
      search_direction_x0s[0] = 1.0;
      search_direction_as[0] = 0.0;

      ++num_search_directions;
    }
    else
    {
      new_search_direction = search_directions[max_num_search_directions-1];
      for(i=max_num_search_directions-1;i>0;--i)
      {
        search_directions[i] = search_directions[i-1];
	search_direction_x0s[i] = search_direction_x0s[i-1];
	search_direction_as[i] = search_direction_as[i-1];
      }
      search_directions[0] = new_search_direction;
      search_direction_x0s[0] = 1.0;
      search_direction_as[0] = 0.0;
    }

    for(j=0;j<num_vars;++j)
    {
      new_search_direction[j] = vect[j] - old_vect[j];
    }
  }

finish: ;

  force_optimize_stop_flag = FALSE;
  last_line_function_x_valid = FALSE;

  wn_gppush(conj_dir_memgp);
  wn_gpfree();
}
