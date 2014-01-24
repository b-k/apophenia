/** \file 
  adaptive rejection metropolis sampling */

/** (C) Wally Gilks; see documentation below for details.
  Adaptations for Apophenia (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"

#define XEPS  0.00001            /* critical relative x-value difference */
#define YEPS  0.1                /* critical y-value difference */
#define EYEPS 0.001              /* critical relative exp(y) difference */
#define YCEIL 50.                /* maximum y avoiding overflow in exp(y) */

/* declarations for functions defined in this file  (minus those in arms.h). */
void invert(double prob, arms_state *env, POINT *p);
int test(arms_state *state, POINT *p,  apop_arms_settings *params, gsl_rng *r);
int update(arms_state *state, POINT *p,  apop_arms_settings *params);
static void cumulate(arms_state *env);
int meet (POINT *q, arms_state *state, apop_arms_settings *params);
double area(POINT *q);
double expshift(double y, double y0);
double logshift(double y, double y0);
double perfunc(apop_arms_settings*, double x);
void display(FILE *f, arms_state *env, apop_arms_settings *);
int initial (apop_arms_settings* params, arms_state *state);

Apop_settings_copy(apop_arms,
    out->state = malloc(sizeof(arms_state));
    *out->state = *in->state;
)

Apop_settings_free(apop_arms,
    if (in->state){
        free(in->state->p);
        free(in->state);
	}
)

Apop_settings_init(apop_arms,
    if ((in.xl || in.xr) && !in.xinit)
        out->xinit = (double []) {in.xl+GSL_DBL_EPSILON, (in.xl+in.xr)/2., in.xr-GSL_DBL_EPSILON};
    else{
        //Apop_varad_set(xinit, ((double []) {0, 0.5, 1}));
        Apop_varad_set(xinit, ((double []) {-1, 0, 1}));
    }
    Apop_varad_set(ninit, 3);
    Apop_varad_set(xl, GSL_MIN(out->xinit[0]/10., out->xinit[0]*10)-.1);
    Apop_varad_set(xr, GSL_MAX(out->xinit[out->ninit-1]/10., out->xinit[out->ninit-1]*10)+.1);
    Apop_varad_set(convex, 0);
    Apop_varad_set(npoint, 100);
    Apop_varad_set(do_metro, 'y');
    Apop_varad_set(xprev, (out->xinit[0]+out->xinit[out->ninit-1])/2.);
    Apop_varad_set(neval, 1000);
    Apop_assert(out->model, "the model input (e.g.: .model = parent_model) is mandatory.");

  // allocate the state 
    out->state = malloc(sizeof(arms_state));
    Apop_assert(out->state, "Malloc failed. Out of memory?");
    *out->state = (arms_state) { };
    int err = initial(out, out->state);
    Apop_assert_c(!err, NULL, 0, "init failed, error %i. Returning NULL", err);

  /* finish setting up metropolis struct (can only do this after setting up env) */
    if(out->do_metro=='y'){
        /* I don't understand why this is needed.
          if((params->xprev < params->xl) || (params->xprev > params->xr))
            apop_assert(0, 1007, 0, 's', "previous Markov chain iterate out of range")*/
        out->state->metro_xprev = out->xprev;
        out->state->metro_yprev = perfunc(out,out->xprev);
        assert(isfinite(out->state->metro_xprev));
        assert(isfinite(out->state->metro_yprev));
    }
)

void distract_doxygen_arms(){/*Doxygen gets thrown by the settings macros. This decoy function is a workaround. */}

/** \brief Adaptive rejection metropolis sampling.

This is a function to make random draws from any univariate distribution (more or less).

The author, Wally Gilks, explains on 
http://www.amsta.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html, that
``ARS works by constructing an envelope function of the log of the target density, which is then used in rejection sampling (see, for example,  Ripley, 1987). Whenever a point is rejected by ARS, the envelope is updated to correspond more closely to the true log density, thereby reducing the chance of rejecting subsequent points. Fewer ARS rejection steps implies fewer point-evaluations of the log density.''

\li It accepts only functions with univariate inputs. I.e., it will put a single value in the vector part of a \ref apop_data set, and then evaluate the log likelihood at that point.

\li It is currently the default for the \ref apop_draw function, so you can just call that if you prefer.

\li There are a great number of parameters, in the \c apop_arms_settings structure.  The structure also holds a history of the points tested to date. That means that the system will be more accurate as more draws are made. It also means that if the parameters change, or you use \ref apop_model_copy, you should call <tt>Apop_settings_rm_group(your_model, apop_arms)</tt> to clear the model of points that are not valid for a different situation.

\li See \ref apop_arms_settings for the list of parameters that you may want to set, via a form like <tt>apop_model_add_group(your_model, apop_arms, .model=your_model, .xl=8, .xr =14);</tt>.  The \c model element is mandatory; you'll get a run-time complaint if you forget it.
  */
int apop_arms_draw (double *out, gsl_rng *r, apop_model *m){
    apop_arms_settings *params = Apop_settings_get_group(m, apop_arms);
    if (!params) params = Apop_model_add_group(m, apop_arms, .model=m);
  POINT pwork;        /* a working point, not yet incorporated in envelope */
  int msamp=0;        /* the number of x-values currently sampled */
  arms_state *state = params->state; 
  /* now do adaptive rejection */
  do {
    // Sample a new point from piecewise exponential envelope 
    double prob = gsl_rng_uniform(r);
    /* get x-value correponding to a cumulative probability prob */
    assert(isfinite(state->p->x));
    assert(isfinite(state->p->y));
    invert(prob,state,&pwork);

    /* perform rejection (and perhaps metropolis) tests */
    int i = test(state,&pwork, params, r);
    if (i == 1){ // point accepted 
        Apop_notify(3, " point accepted.");
        *out = pwork.x;
        assert(isfinite(pwork.x));
        return 0;
    } else 
      Apop_stopif(i!=0, return 1,-5, "envelope error - violation without metropolis");
    msamp ++;
    Apop_notify(3, " point rejected.");
  } while (msamp < 1e3);
  Apop_notify(1, "I just rejected 1,000 samples. Something is wrong.");
  return 0;
}

int initial (apop_arms_settings* params,  arms_state *env){
// to set up initial envelope

  POINT *q;
  int mpoint = 2*params->ninit + 1;

  Apop_assert_c(params->ninit>=3, 1001, 0, "too few initial points");
  Apop_assert_c(params->npoint >= mpoint, 1002, 0, "too many initial points");
  Apop_assert_c((params->xinit[0] >= params->xl) && (params->xinit[params->ninit-1] <= params->xr),
                     1003, 0, "initial points do not satisfy bounds");
  for(int i=1; i<params->ninit; i++)
      Apop_assert_c(params->xinit[i] > params->xinit[i-1], 1004, 0, "data not ordered");
  Apop_assert_c(params->convex >= 0.0, 1008, 0, "negative convexity parameter");

  env->convex = &params->convex; // copy convexity address to env
  params->neval = 0; // initialise current number of function evaluations

  /* set up space for envelope POINTs */
  env->npoint = params->npoint;
  env->p = malloc(params->npoint*sizeof(POINT));
  Apop_assert(env->p, "malloc of env->p failed. Out of space?");

  /* set up envelope POINTs */
  q = env->p;
  q->x = params->xl; // left bound
  q->f = 0;
  q->pl = NULL;
  q->pr = q+1;
  for(int j=1, k=0; j<mpoint-1; j++){
    q++;
    if(j%2){
        /* point on log density */
        q->x = params->xinit[k++];
        q->y = perfunc(params,q->x);
        Apop_assert(isfinite(q->x), "the initial param is %g", q->x);
        Apop_assert(isfinite(q->y), "f(an initial parameter)= %g", q->y);
        q->f = 1;
    } else // intersection point
        q->f = 0;
    q->pl = q-1;
    q->pr = q+1;
  }
  /* right bound */
  q++;
  q->x = params->xr;
  q->f = 0;
  q->pl = q-1;
  q->pr = NULL;

  assert(isfinite(q->x));
  /* calculate intersection points */
  q = env->p;
  for (int j=0; j<mpoint; j=j+2, q=q+2)
    Apop_assert_c(!meet(q,env, params), 2000, 0, "envelope violation without metropolis");

  cumulate(env); // exponentiate and integrate envelope
  env->cpoint = mpoint; // note number of POINTs currently in envelope
  return 0;
}

void invert(double prob, arms_state *env, POINT *p){
/* to obtain a point corresponding to a given cumulative probability   
   prob    : cumulative probability under envelope   
   *env    : envelope attributes   
   *p      : a working POINT to hold the sampled value */

  double u,xl=0,xr=0,yl,yr,eyl,eyr,prop;

  /* find rightmost point in envelope */
  POINT *q = env->p;
  while(q->pr != NULL)q = q->pr;

  /* find exponential piece containing point implied by prob */
  u = prob * q->cum;
  while(q->pl->cum > u)q = q->pl;

  /* piece found: set left and right POINTs of p, etc. */
  p->pl = q->pl;
  p->pr = q;
  p->f = 0;
  p->cum = u;

  /* calculate proportion of way through integral within this piece */
  prop = (u - q->pl->cum) / (q->cum - q->pl->cum);

  /* get the required x-value */
  if (q->pl->x == q->x){
    /* interval is of zero length */
    p->x = q->x;
    p->y = q->y;
    p->ey = q->ey;
  } else {
    xl = q->pl->x;
    xr = q->x;
    yl = q->pl->y;
    yr = q->y;
    eyl = q->pl->ey;
    eyr = q->ey;
    if(fabs(yr - yl) < YEPS){
      /* linear approximation was used in integration in function cumulate */
      if(fabs(eyr - eyl) > EYEPS*fabs(eyr + eyl))
        p->x = xl + ((xr - xl)/(eyr - eyl)) * (-eyl + sqrt((1. - prop)*eyl*eyl + prop*eyr*eyr));
      else 
        p->x = xl + (xr - xl)*prop;
      p->ey = ((p->x - xl)/(xr - xl)) * (eyr - eyl) + eyl;
      p->y = logshift(p->ey, env->ymax);
    } else {
      /* piece was integrated exactly in function cumulate */
      p->x = xl + ((xr - xl)/(yr - yl))
	      * (-yl + logshift(((1.-prop)*eyl + prop*eyr), env->ymax));
      p->y = ((p->x - xl)/(xr - xl)) * (yr - yl) + yl;
      p->ey = expshift(p->y, env->ymax);
    }
  }
  assert(isfinite(p->x));
  assert(isfinite(p->y));
  assert(isfinite(q->x));
  assert(isfinite(q->y));

  /* guard against imprecision yielding point outside interval */
  if ((p->x < xl) || (p->x > xr)) exit(1);
}

int test(arms_state *env, POINT *p, apop_arms_settings *params, gsl_rng *r){
/* to perform rejection, squeezing, and metropolis tests   
   *env        : state data
   *p            : point to be tested   */
assert(p->pl && p->pr);

  double u,y,ysqueez,ynew,yold,znew,zold,w;
  POINT *ql,*qr;
  
  /* for rejection test */
  u = gsl_rng_uniform(r) * p->ey;
  y = logshift(u,env->ymax);

  if(params->do_metro !='y' && (p->pl->pl != NULL) && (p->pr->pr != NULL)){
    /* perform squeezing test */
    ql = p->pl->f ? p->pl : p->pl->pl;
    qr = p->pr->f ? p->pr : p->pr->pr;
    ysqueez = (qr->y * (p->x - ql->x) + ql->y * (qr->x - p->x))
               /(qr->x - ql->x);
    if(y <= ysqueez) // accept point at squeezing step
        return 1;
  }

  /* evaluate log density at point to be tested */
  ynew = perfunc(params,p->x);
assert(isfinite(p->x));
assert(p->pl && p->pr);
Apop_notify(3, "tested (%g, %g); ", p->x, ynew);
  
  /* perform rejection test */
  if(params->do_metro != 'y' || (params->do_metro == 'y' && (y >= ynew))){
    /* update envelope */
    p->y = ynew;
    p->ey = expshift(p->y,env->ymax);
    p->f = 1;
    if(update(env,p, params)) 
        Apop_assert_c(0, -1, 0, "envelope violation without metropolis");
    /* perform rejection test: accept iff y < ynew */
    return (y < ynew);
  }

  /* continue with metropolis step */
  yold = env->metro_yprev;
  /* find envelope piece containing metrop->xprev */
  ql = env->p;
  while(ql->pl != NULL) ql = ql->pl;
  while(ql->pr->x < env->metro_xprev) ql = ql->pr;
  qr = ql->pr;
  /* calculate height of envelope at metrop->xprev */
  w = (env->metro_xprev - ql->x)/(qr->x - ql->x);
  zold = ql->y + w*(qr->y - ql->y);
  znew = p->y;
  if(yold < zold)zold = yold;
  if(ynew < znew)znew = ynew;
  w = ynew-znew-yold+zold;
  w = GSL_MIN(w, 0.0);
  w = (w > -YCEIL) ?  exp(w) : 0.0;
  u = gsl_rng_uniform(r);
  if(u > w){
      /* metropolis says don't move, so replace current point with previous */
      /* markov chain iterate */
      p->x = env->metro_xprev;
      p->y = env->metro_yprev;
      Apop_notify(3, "metro step (%g) rejected with w=%g, "
                "ynew=%g, yold=%g, znew = %g, zold=%g; ", p->x, w, ynew, yold, znew, zold);
      p->ey = expshift(p->y,env->ymax);
assert(isfinite(p->x));
assert(isfinite(p->y));
assert(isfinite(p->ey));
      p->f = 1;
      p->pl = ql;
      p->pr = qr;
  } else {
      /* trial point accepted by metropolis, so update previous markov */
      /* chain iterate */
      env->metro_xprev = p->x;
      env->metro_yprev = ynew;
  }
  return 1;
}

int update(arms_state *env, POINT *p, apop_arms_settings *params){
/* to update envelope to incorporate new point on log density
   *env          : state information
   *p            : point to be incorporated 
*/

  POINT *m,*ql,*qr,*q;

  if(!(p->f) || (env->cpoint > env->npoint - 2))
    /* y-value has not been evaluated or no room for further points */
    return 0; // ignore this point

  /* copy working POINT p to a new POINT q */
  q = env->p + env->cpoint++;
  q->x = p->x;
  q->y = p->y;
  q->f = 1;

  /* allocate an unused POINT for a new intersection */
  m = env->p + env->cpoint++;
  m->f = 0;
  if((p->pl->f) && !(p->pr->f)){
    /* left end of piece is on log density; right end is not */
    /* set up new intersection in interval between p->pl and p */
    m->pl = p->pl;
    m->pr = q;
    q->pl = m;
    q->pr = p->pr;
    m->pl->pr = m;
    q->pr->pl = q;
  } else if (!(p->pl->f) && (p->pr->f)){
    /* left end of interval is not on log density; right end is */
    /* set up new intersection in interval between p and p->pr */
    m->pr = p->pr;
    m->pl = q;
    q->pr = m;
    q->pl = p->pl;
    m->pr->pl = m;
    q->pl->pr = q;
  } else exit(10);// this should be impossible 

  /* now adjust position of q within interval if too close to an endpoint */
  ql = q->pl->pl ? q->pl->pl : q->pl;
  qr = q->pr->pr ? q->pr->pr : q->pr;
  if (q->x < (1. - XEPS) * ql->x + XEPS * qr->x){
    /* q too close to left end of interval */
    q->x = (1. - XEPS) * ql->x + XEPS * qr->x;
    q->y = perfunc(params,q->x);
  } else if (q->x > XEPS * ql->x + (1. - XEPS) * qr->x){
    /* q too close to right end of interval */
    q->x = XEPS * ql->x + (1. - XEPS) * qr->x;
    q->y = perfunc(params,q->x);
  }

  /* revise intersection points */
  if(meet(q->pl,env, params) /* envelope violations without metropolis */
        || meet(q->pr,env, params) 
        || (q->pl->pl != NULL && meet(q->pl->pl->pl,env, params))
        || (q->pr->pr != NULL && meet(q->pr->pr->pr,env, params)))
     return 1;

  /* exponentiate and integrate new envelope */
  cumulate(env);
  return 0;
}

static void cumulate(arms_state *env){
/* to exponentiate and integrate envelope */
/* *env     : envelope attributes */

  POINT *q,*qlmost;

  qlmost = env->p;
  /* find left end of envelope */
  while(qlmost->pl) qlmost = qlmost->pl;

  /* find maximum y-value: search envelope */
  env->ymax = qlmost->y;
  for(q = qlmost->pr; q != NULL; q = q->pr)
    if(q->y > env->ymax)
        env->ymax = q->y;

  /* exponentiate envelope */
  for(q = qlmost; q != NULL; q = q->pr)
      q->ey = expshift(q->y,env->ymax);

  /* integrate exponentiated envelope */
  qlmost->cum = 0.;
  for(q = qlmost->pr; q != NULL; q = q->pr)
      q->cum = q->pl->cum + area(q);
}

int meet (POINT *q, arms_state *env, apop_arms_settings *params){ 
/* To find where two chords intersect 
   q         : to store point of intersection 
   *env      : state attributes 
*/
  double gl=0,gr=0,grl=0,dl=0,dr=0;
  int il=0,ir=0,irl=0;

  Apop_assert(!(q->f), "error 30: this is not an intersection point.");

  /* calculate coordinates of point of intersection */
  if ((q->pl != NULL) && (q->pl->pl->pl != NULL)){
      /* chord gradient can be calculated at left end of interval */
      gl = (q->pl->y - q->pl->pl->pl->y)/(q->pl->x - q->pl->pl->pl->x);
      il = 1;
  } else // no chord gradient on left 
      il = 0;

  if ((q->pr != NULL) && (q->pr->pr->pr != NULL)){
      /* chord gradient can be calculated at right end of interval */
      gr = (q->pr->y - q->pr->pr->pr->y)/(q->pr->x - q->pr->pr->pr->x);
      ir = 1;
  } else // no chord gradient on right
      ir = 0;

  if ((q->pl != NULL) && (q->pr != NULL)){
      /* chord gradient can be calculated across interval */
      grl = (q->pr->y - q->pl->y)/(q->pr->x - q->pl->x);
      irl = 1;
  } else 
      irl = 0;

  if(irl && il && (gl<grl)){
    /* convexity on left exceeds current threshold */
    if(params->do_metro !='y') // envelope violation without metropolis
        return 1;
    gl = gl + (1.0 + *(env->convex)) * (grl - gl); // adjust left gradient 
  }

  if(irl && ir && (gr>grl)){
    /* convexity on right exceeds current threshold */
    if(params->do_metro !='y') // envelope violation without metropolis 
        return 1;
    gr = gr + (1.0 + *(env->convex)) * (grl - gr); // adjust right gradient 
  }

  if(il && irl){
    dr = (gl - grl) * (q->pr->x - q->pl->x);
    if(dr < YEPS) // adjust dr to avoid numerical problems
      dr = YEPS;
  }

  if(ir && irl){
    dl = (grl - gr) * (q->pr->x - q->pl->x);
    if(dl < YEPS) // adjust dl to avoid numerical problems 
      dl = YEPS;
  }

  if(il && ir && irl){
    /* gradients on both sides */
    q->x = (dl * q->pr->x + dr * q->pl->x)/(dl + dr);
    q->y = (dl * q->pr->y + dr * q->pl->y + dl * dr)/(dl + dr);
  } else if (il && irl){
    /* gradient only on left side, but not right hand bound */
    q->x = q->pr->x;
    q->y = q->pr->y + dr;
  } else if (ir && irl){
    /* gradient only on right side, but not left hand bound */
    q->x = q->pl->x;
    q->y = q->pl->y + dl;
  } else if (il)
    q->y = q->pl->y + gl * (q->x - q->pl->x); // right hand bound 
  else if (ir)
    q->y = q->pr->y - gr * (q->pr->x - q->x); // left hand bound 
  else 
     Apop_assert(0, "error 31: gradient on neither side - should be impossible.");
  if(((q->pl != NULL) && (q->x < q->pl->x)) ||
     ((q->pr != NULL) && (q->x > q->pr->x))){
     Apop_assert(0, "error 32: intersection point outside interval (through imprecision)");
  }
  return 0; // successful exit : intersection has been calculated
}

double area(POINT *q){
/* To integrate piece of exponentiated envelope to left of POINT q */ 

  if(q->pl == NULL) // this is leftmost point in envelope 
      exit(1);
  if(q->pl->x == q->x) // interval is zero length
      return 0.;
  if (fabs(q->y - q->pl->y) < YEPS) // integrate straight line piece
      return 0.5*(q->ey + q->pl->ey)*(q->x - q->pl->x);
  // integrate exponential piece 
  return ((q->ey - q->pl->ey)/(q->y - q->pl->y))*(q->x - q->pl->x);
}

double expshift(double y, double y0) {
/* to exponentiate shifted y without underflow */
  if (y - y0 > -2.0 * YCEIL)
      return exp(y - y0 + YCEIL);
  else
      return 0.0;
}

double logshift(double y, double y0){
/* inverse of function expshift */
  return (log(y) + y0 - YCEIL);
}

double perfunc(apop_arms_settings *params, double x){
// to evaluate log density and increment count of evaluations 
    Staticdef( apop_data *, d , apop_data_alloc(1));
    d->vector->data[0] = x;
  double y = apop_log_likelihood(d, params->model);
  Apop_assert(isfinite(y), "Evaluating the log likelihood of %g returned %g.", x, y);
  (params->neval)++; // increment count of function evaluations
  return y;
}
