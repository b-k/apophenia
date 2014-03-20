/* Functions to calculate loess regressions. 

Provenance:
   * Originally written in FORTRAN `77 by Cleveland, Devlin, Grosse, and Shyu, 1988
   * Some C wrapper written by Cleveland, Grosse, and Shyu, 1992
   * Punched through f2c, munged into one file, and heavily edited by Klemens, 2009, 2011
   * Legal & ref.s: Documentation from CG&S state that their code is
            public domain. See http://netlib.org/a/cloess.pdf (esp. if
            you hope to modify the below). You can still get the `92
            version from http://netlib.org/a/dloess (a shell script
            that unpacks into everything you need). Most BK edits (c)
            2009, licensed under the modified GPLv2; see COPYING and
            COPYING2. Those BK edits made during time working as a gov't
            employee are public domain.

    By the way, search the code for execnt: many functions will let you
    query how many times they have been hit, which you might find to be useful.

\amodel apop_loess Regression via loess smoothing

This uses a somewhat black-box routine, first written by Chamberlain, Devlin, Grosse,
and Shyu in 1988, to fit a smoothed series of quadratic curves to the input data,
thus producing a curve more closely fitting than a simple regression would.

The curve is basically impossible to describe using a short list of parameters, so the
representation is in the form of the \c predicted vector of the \c expected data set;
see below.

From the 1992 manual for the package:
``The method we will use to fit local regression models is called {\em loess}, which
is short for local regression, and was chosen as the name since a loess is a deposit
of fine clay or silt along a river valley, and thus is a surface of sorts. The word
comes from the German löss, and is pronounced löíss.''

\adoc    Input_format  
The data is basically OLS-like:                     
the first column of the data is the dependent variable to be explained; subsequent
variables are the independent explanatory variables.  Thus, your input data can either
have a dependent vector plus explanatory matrix, or a matrix where the first column
is the dependent variable.

Unlike with OLS, I won't move your original data, and I won't add a <b>1</b>, because
that's not really the loess custom. You can of course set up your data that way if
you like.

If your data set has a weights vector, I'll use it.

In any case, all data is copied into the model's \ref apop_loess_settings. The code
is primarily FORTRAN code from 1988 converted to C; the data thus has to be converted
into a relatively obsolete internal format.


\adoc    Parameter_format  The parameter vector is unused. 
\adoc    estimated_parameters None.  
\adoc    estimated_settings
The \ref apop_loess_settings is filled with results (and internal processing cruft). The
\c out_model->info data set has a table giving the actual, \c predicted, and \c residual
columns, which is probably what you were looking for.  Try:
        \code
        apop_data_show(apop_data_get_page(output_model->info, "<Predicted>"));
        \endcode
        
\adoc    Predict 
Fills in the zeroth column (ignoring and overwriting any data there), and at the data's <tt>->more</tt> pointer, adds an \ref
apop_data set named "Confidence" (i.e., 
\code
!strcmp(outdata->more->names->title, "Confidence") == 1.
\endcode 

This routine is in beta testing.

\adoc    settings \ref apop_loess_settings */

#include "apop_internal.h"

////////////a few lines from f2c.h
#define TRUE_ (1)
#define FALSE_ (0)
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

typedef long int logical;
typedef long int integer;

#define Calloc(n,t)	(t *)calloc((unsigned)(n),sizeof(t))   // From #include "S.h"
#define Warning(msg) Apop_assert_c(0, , 0 , "%s", msg);

static integer c__0 = 0;
static integer c__1 = 1;
static integer c__15 = 15;
static integer c__2 = 2;
static integer c__21 = 21;
double doublepluszero = 0;

//I'm using the GSL's blas system. These are pass-through functions that
//save me the trouble of slogging through the code and making substitutions.
void dswap(const integer N, double *x, double *y){ cblas_dswap(N, x, 1, y, 1);}

double dnrm2(const integer *N, const double *x){ return cblas_dnrm2(*N, x, 1); }

double ddot_(const integer *N, const double *x, const integer *incx, const double *y, const integer *incy){
		return cblas_ddot(*N, x, *incx, y, *incy);}

void daxpy_(integer *N, const double *alpha, const double *x, integer *incx, double *y, integer *incy){
		cblas_daxpy(*N, *alpha, x, *incx, y, *incy);}

void dcopy_(const integer *N, const double *x, const integer incx, double *y, const integer incy){
		cblas_dcopy(*N, x, incx, y, incy); }

void dscal(integer *N, const double alpha, double *x){ cblas_dscal(*N, alpha, x, 1);}

static void drotg_(double *a,double *b, double *c, double *s){ cblas_drotg(a,b,c,s);}

void drot_(const integer *N, double *x, const integer *incx, double *y, const integer *incy, const double *c, const double *s){
		cblas_drot(*N, x, *incx, y, *incy, *c, *s);}

//Procedures from libf2c.
static integer pow_ii(integer x, integer n) { return gsl_pow_int(x, n); }

static double d_sign(double *a, double *b) {
    double x = (*a >= 0 ? *a : - *a);
    return( *b >= 0 ? x : -x);
}

/*//// dqrsl.f -- translated by f2c (version 20061008).

     dqrsl applies the output of dqrdc to compute coordinate transformations,
     projections, and least squares solutions.  for k .le. min(n,p), let xk be the matrix

            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k))) 

     formed from columns jpvt(1), ... ,jpvt(k) of the original n x p matrix x that was
     input to dqrdc (if no pivoting was done, xk consists of the first k columns of x
     in their original order).  dqrdc produces a factored orthogonal matrix q and an
     upper triangular matrix r such that

              xk = q * (r)
                         (0)             

     this information is contained in coded form in the arrays x and qraux. 

     on entry 

        x      double precision(ldx,p). 
               x contains the output of dqrdc. 

        ldx    integer. 
               ldx is the leading dimension of the array x. 

        n      integer. 
               n is the number of rows of the matrix xk.  it must have the same value
               as n in dqrdc.

        k      integer. 
               k is the number of columns of the matrix xk.  k must nnot be greater
               than min(n,p), where p is the same as in the calling sequence to dqrdc.

        qraux  double precision(p). 
               qraux contains the auxiliary output from dqrdc. 

        y      double precision(n) 
               y contains an n-vector that is to be manipulated by dqrsl. 

        job    integer. 
               job specifies what is to be computed.  job has the decimal expansion
               abcde, with the following meaning.

                    if a.ne.0, compute qy.
                      if b,c,d, or e .ne. 0, compute qty.
                      if c.ne.0, compute b.
                      if d.ne.0, compute rsd.
                      if e.ne.0, compute xb. 

               note that a request to compute b, rsd, or xb automatically triggers
               the computation of qty, for which an array must be provided in the
               calling sequence.

     on return 

        qy     double precision(n). 
               qy contains q*y, if its computation has been requested. 

        qty    double precision(n). 
               qty contains trans(q)*y, if its computation has been requested.
               here trans(q) is the transpose of the matrix q.

        b      double precision(k) 
               b contains the solution of the least squares problem 

                    minimize norm2(y - xk*b), 

               if its computation has been requested.  (note that if pivoting was
               requested in dqrdc, the j-th component of b will be associated with
               column jpvt(j) of the original matrix x that was input into dqrdc.)

        rsd    double precision(n). 
               rsd contains the least squares residual y - xk*b, if its computation
               has been requested.  rsd is also the orthogonal projection of y onto
               the orthogonal complement of the column space of xk.

        xb     double precision(n). 
               xb contains the least squares approximation xk*b, if its computation
               has been requested.  xb is also the orthogonal projection of y onto
               the column space of x.

        info   integer. 
               info is zero unless the computation of b has been requested and r is
               exactly singular.  in this case, info is the index of the first zero
               diagonal element of r and b is left unaltered.

     The parameters qy, qty, b, rsd, and xb are not referenced if their computation
     is not requested and in this case can be replaced by dummy variables in the
     calling program.  to save storage, the user may in some cases use the same array
     for different parameters in the calling sequence.  a frequently occurring example
     is when one wishes to compute any of b, rsd, or xb and does not need y or qty.
     in this case one may identify y, qty, and one of b, rsd, or xb, while providing
     separate arrays for anything else that is to be computed.  thus the calling sequence

          call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info) 

     will result in the computation of b and rsd, with rsd overwriting y.  More
     generally, each item in the following list contains groups of permissible
     identifications for a single calling sequence.

          1. (y,qty,b) (rsd) (xb) (qy)
            2. (y,qty,rsd) (b) (xb) (qy)
            3. (y,qty,xb) (b) (rsd) (qy)
            4. (y,qy) (qty,b) (rsd) (xb)
            5. (y,qy) (qty,rsd) (b) (xb)
            6. (y,qy) (qty,xb) (b) (rsd) 

     in any group the value returned in the array allocated to the group corresponds
     to the last member of the group.

     linpack. this version dated 08/14/78 . 
     g.w. stewart, university of maryland, argonne national lab. */
static void dqrsl_(double *x, integer *ldx, integer *n, integer * k, double *qraux, double *y, 
        double *qy, double *qty, double *b, double *rsd, double *xb, integer job, integer *info) {

    integer x_dim1, i__1, i__2;
    static integer i__, j;
    static double t, temp;
    static integer jj, ju, kp1;
    static logical cb, cr, cxb, cqy, cqty;

    x_dim1 = *ldx;
    x -= 1 + x_dim1;
    --qraux;
    --y;
    --qy;
    --qty;
    --b;
    --rsd;
    --xb;

    *info = 0; /*     set info flag. */

    /*     determine what is to be computed. */
    cqy = job / 10000 != 0;
    cqty = job % 10000 != 0;
    cb = job % 1000 / 100 != 0;
    cr = job % 100 / 10 != 0;
    cxb = job % 10 != 0;
    ju = min(*k,*n - 1);

    /*     special action when n=1. */
    if (ju != 0) goto L40;
    if (cqy)     qy[1] = y[1];
    if (cqty)    qty[1] = y[1];
    if (cxb)     xb[1] = y[1];
    if (!cb)     goto L30;
    if (x[x_dim1 + 1] != 0.) goto L10;
    *info = 1;
    goto L20;
L10:
    b[1] = y[1] / x[x_dim1 + 1];
L20:
L30:
    if (cr) rsd[1] = 0.;
    return;
L40:
    /*        set up to compute qy or qty. */
    if (cqy)  dcopy_(n, &y[1], 1, &qy[1], 1);
    if (cqty) dcopy_(n, &y[1], 1, &qty[1], 1);
    if (!cqy) goto L70;

    /*           compute qy. */
    for (jj = 1; jj <= ju; ++jj) {
        j = ju - jj + 1;
        if (qraux[j] == 0.)
            continue;
        temp = x[j + j * x_dim1];
        x[j + j * x_dim1] = qraux[j];
        i__2 = *n - j + 1;
        t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &qy[j], &c__1) / x[j + j * x_dim1];
        i__2 = *n - j + 1;
        daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qy[j], &c__1); x[j + j * x_dim1] = temp;
    }
L70:
    if (cqty) /*           compute trans(q)*y. */
        for (j = 1; j <= ju; ++j) {
            if (qraux[j] == 0.)
                continue;
            temp = x[j + j * x_dim1];
            x[j + j * x_dim1] = qraux[j];
            i__2 = *n - j + 1;
            t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &qty[j], &c__1) / x[j + j * x_dim1];
            i__2 = *n - j + 1;
            daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qty[j], &c__1);
            x[j + j * x_dim1] = temp;
        }

    /*        set up to compute b, rsd, or xb. */
    if (cb) dcopy_(k, &qty[1], 1, &b[1], 1);
    kp1 = *k + 1;
    if (cxb) dcopy_(k, &qty[1], 1, &xb[1], 1);
    if (cr && *k < *n) {
        i__1 = *n - *k;
        dcopy_(&i__1, &qty[kp1], 1, &rsd[kp1], 1);
    }
    if (! cxb || kp1 > *n)
        goto L120;
    for (i__ = kp1; i__ <= *n; ++i__)
        xb[i__] = 0.;
L120:
    if (! cr)
        goto L140;
    for (i__ = 1; i__ <= *k; ++i__)
        rsd[i__] = 0.;
L140:
    if (! cb)
        goto L190;

    /*           compute b. */
    for (jj = 1; jj <= *k; ++jj) {
        j = *k - jj + 1;
        if (x[j + j * x_dim1] != 0.)
            goto L150;
        *info = j;
        break;
    L150:
        b[j] /= x[j + j * x_dim1];
        if (j == 1)
            continue;
        t = -b[j];
        i__2 = j - 1;
        daxpy_(&i__2, &t, &x[j * x_dim1 + 1], &c__1, &b[1], &c__1);
    }
L190:
    if (! cr && ! cxb)
        return;
    /*           compute rsd or xb as required. */
    for (jj = 1; jj <= ju; ++jj) {
        j = ju - jj + 1;
        if (qraux[j] == 0.)
            continue;
        temp = x[j + j * x_dim1];
        x[j + j * x_dim1] = qraux[j];
        if (cr) {
            i__2 = *n - j + 1;
            t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1) / x[j + j * x_dim1];
            i__2 = *n - j + 1;
            daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1);
        }
        if (cxb) {
            i__2 = *n - j + 1;
            t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &xb[j], &c__1) / x[j + j * x_dim1];
            i__2 = *n - j + 1;
            daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &xb[j], &c__1);
        }
        x[j + j * x_dim1] = temp;
    }
} /* dqrsl_ */

/*//// dsvdc.f -- translated by f2c (version 20061008).  

     dsvdc is a subroutine to reduce a double precision nxp matrix x by orthogonal
     transformations u and v to diagonal form.  the diagonal elements s(i) are the
     singular values of x.  the columns of u are the corresponding left singular vectors,
     and the columns of v the right singular vectors.

     on entry 

         x         double precision(ldx,p), where ldx.ge.n. 
                   x contains the matrix whose singular value decomposition is to
                   be computed.  x is destroyed by dsvdc.

         ldx       integer. 
                   ldx is the leading dimension of the array x. 

         n         integer. 
                   n is the number of rows of the matrix x. 

         p         integer. 
                   p is the number of columns of the matrix x. 

         ldu       integer. 
                   ldu is the leading dimension of the array u.  (see below).

         ldv       integer. 
                   ldv is the leading dimension of the array v.  (see below).

         work      double precision(n). 
                   work is a scratch array. 

         job       integer. 
                   job controls the computation of the singular vectors.  It has the
                   decimal expansion ab with the following meaning

                    a.eq.0    do not compute the left singular vectors.
                      a.eq.1    return the n left singular vectors in u.
                      a.ge.2    return the first min(n,p) singular vectors in u. 
                      b.eq.0    do not compute the right singular vectors. 
                      b.eq.1    return the right singular vectors in v. 

     on return 

         s         double precision(mm), where mm=min(n+1,p).
                   the first min(n,p) entries of s contain the singular values of x
                   arranged in descending order of magnitude.

         e         double precision(p), 
                   e ordinarily contains zeros.  however see the discussion of info
                   for exceptions.

         u         double precision(ldu,k), where ldu.ge.n.  if 
                                   joba.eq.1 then k.eq.n, if joba.ge.2 
                                   then k.eq.min(n,p). 
                   u contains the matrix of left singular vectors.  u is not referenced
                   if joba.eq.0.  if n.le.p or if joba.eq.2, then u may be identified
                   with x in the subroutine call.

         v         double precision(ldv,p), where ldv.ge.p. 
                   v contains the matrix of right singular vectors.  v is not referenced
                   if job.eq.0.  if p.le.n, then v may be identified with x in the
                   subroutine call.

         info      integer. 
                   the singular values (and their corresponding singular vectors)
                   s(info+1),s(info+2),...,s(m) are correct (here m=min(n,p)).  thus if
                   info.eq.0, all the singular values and their vectors are correct.
                   in any event, the matrix b = trans(u)*x*v is the bidiagonal matrix
                   with the elements of s on its diagonal and the elements of e on its
                   super-diagonal (trans(u) is the transpose of u).  thus the singular
                   values of x and b are the same.

     linpack. this version dated 08/14/78 . 
              correction made to shift 2/84. 
     g.w. stewart, university of maryland, argonne national lab. 

     Modified 2000-12-28 to use a relative convergence test, 
     as this was infinite-looping on ix86. 

     dsvdc uses the following functions and subprograms. */

static void dsvdc_(double *x, integer *ldx, integer *n, integer * p, double *s,
        double *e, double *u, integer *ldu,
        double *v, integer *ldv, double *work, integer *job, integer * info) {
    integer x_dim1, u_dim1, v_dim1, i__2, i__3;
    double d__1;
    static double b, c__, f, g, t, t1, el, cs, sl, sm, sn, acc, emm1, smm1;
    static double test, scale, shift, ztest;
    static integer i__, j, k, l, m, kk, ll, mm, ls, lu, lm1, mm1, lp1, mp1, nct, ncu, lls, nrt;
    static integer kase, jobu, iter, nctp1, nrtp1, maxit;
    static logical wantu, wantv;

    x_dim1 = *ldx;
    x -= 1 + x_dim1;
    --s;
    --e;
    u_dim1 = *ldu;
    u -= 1 + u_dim1;
    v_dim1 = *ldv;
    v -= 1 + v_dim1;
    --work;
    l = 0;
    ls = 0;
    maxit = 30; // set the maximum number of iterations.

    /*     determine what is to be computed. */
    wantu = FALSE_;
    wantv = FALSE_;
    jobu = *job % 100 / 10;
    ncu = *n;
    if (jobu > 1)
        ncu = min(*n,*p);
    if (jobu != 0)
        wantu = TRUE_;
    if (*job % 10 != 0)
        wantv = TRUE_;

    /*     reduce x to bidiagonal form, storing the diagonal elements */
    /*     in s and the super-diagonal elements in e. */
    *info = 0;
    nct = min(*n - 1,*p);
    nrt = max(0,min(*p - 2,*n));
    lu = max(nct,nrt);
    if (lu < 1)
        goto L170;
    for (l = 1; l <= lu; ++l) {
        lp1 = l + 1;
        if (l <= nct) {
            /*           compute the transformation for the l-th column and */
            /*           place the l-th diagonal in s(l). */
            i__2 = *n - l + 1;
            s[l] = dnrm2(&i__2, &x[l + l * x_dim1]);
            if (s[l] != 0.) {
                if (x[l + l * x_dim1] != 0.) {
                    s[l] = d_sign(&s[l], &x[l + l * x_dim1]);
                }
                i__2 = *n - l + 1;
                dscal(&i__2, 1./s[l], &x[l + l * x_dim1]);
                x[l + l * x_dim1] += 1.;
            }
            s[l] = -s[l];
        }
        if (*p >= lp1)
            for (j = lp1; j <= *p; ++j) {
                if ((l <= nct) && (s[l] != 0.)) {
                    /*              apply the transformation. */
                    i__3 = *n - l + 1;
                    t = -ddot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &c__1) / x[l + l * x_dim1];
                    i__3 = *n - l + 1;
                    daxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &c__1);
                }
                /*           place the l-th row of x into  e for the */
                /*           subsequent calculation of the row transformation. */
                e[j] = x[l + j * x_dim1];
            }
        if (!(! wantu || l > nct)) /* place the transformation in u for subsequent back multiplication. */
            for (i__ = l; i__ <= *n; ++i__)
                u[i__ + l * u_dim1] = x[i__ + l * x_dim1];
        if (l > nrt)
            continue;

        /*           compute the l-th row transformation and place the */
        /*           l-th super-diagonal in e(l). */
        i__2 = *p - l;
        e[l] = dnrm2(&i__2, &e[lp1]);
        if (e[l] != 0.) {
            if (e[lp1] != 0.)
                e[l] = d_sign(&e[l], &e[lp1]);
            i__2 = *p - l;
            dscal(&i__2, 1./e[l], &e[lp1]);
            e[lp1] += 1.;
        }
        e[l] = -e[l];
        if (!(lp1 > *n || e[l] == 0.)) {
            /*              apply the transformation. */
            for (i__ = lp1; i__ <= *n; ++i__)
                work[i__] = 0.;
            for (j = lp1; j <= *p; ++j) {
                i__3 = *n - l;
                daxpy_(&i__3, &e[j], &x[lp1 + j * x_dim1], &c__1, &work[lp1], &c__1);
            }
            for (j = lp1; j <= *p; ++j) {
                i__3 = *n - l;
                d__1 = -e[j] / e[lp1];
                daxpy_(&i__3, &d__1, &work[lp1], &c__1, &x[lp1 + j * x_dim1], &c__1);
            }
        }
        if (! wantv)
            continue;

        /*  place the transformation in v for subsequent back multiplication. */
        for (i__ = lp1; i__ <= *p; ++i__)
            v[i__ + l * v_dim1] = e[i__];
    }
L170:

    /*     set up the final bidiagonal matrix or order m. */
    m = GSL_MIN(*p, *n + 1);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;
    if (nct < *p)
        s[nctp1] = x[nctp1 + nctp1 * x_dim1];
    if (*n < m)
        s[m] = 0.;
    if (nrtp1 < m)
        e[nrtp1] = x[nrtp1 + m * x_dim1];
    e[m] = 0.;

    /*     if required, generate u. */
    if (wantu) {
        if (ncu >= nctp1)
            for (j = nctp1; j <= ncu; ++j) {
                for (i__ = 1; i__ <= *n; ++i__)
                    u[i__ + j * u_dim1] = 0.;
                u[j + j * u_dim1] = 1.;
            }
        if (nct >= 1)
            for (ll = 1; ll <= nct; ++ll) {
                l = nct - ll + 1;
                if (s[l] != 0.) {
                    lp1 = l + 1;
                    if (ncu >= lp1)
                        for (j = lp1; j <= ncu; ++j) {
                            i__3 = *n - l + 1;
                            t = -ddot_(&i__3, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], & c__1) / u[l + l * u_dim1];
                            i__3 = *n - l + 1;
                            daxpy_(&i__3, &t, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], & c__1);
                        }
                    i__2 = *n - l + 1;
                    dscal(&i__2, -1., &u[l + l * u_dim1]);
                    u[l + l * u_dim1] += 1.;
                    lm1 = l - 1;
                    if (lm1 >= 1)
                        for (i__ = 1; i__ <= lm1; ++i__)
                            u[i__ + l * u_dim1] = 0.;
                    continue;
                }
                for (i__ = 1; i__ <= *n; ++i__)
                    u[i__ + l * u_dim1] = 0.;
                u[l + l * u_dim1] = 1.;
            }
    }

    /*     if it is required, generate v. */
    if (wantv)
        for (ll = 1; ll <= *p; ++ll) {
            l = *p - ll + 1;
            lp1 = l + 1;
            if ((l <= nrt) && (e[l] != 0.))
                for (j = lp1; j <= *p; ++j) {
                    i__3 = *p - l;
                    t = -ddot_(&i__3, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * v_dim1], &c__1) / v[lp1 + l * v_dim1];
                    i__3 = *p - l;
                    daxpy_(&i__3, &t, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * v_dim1], &c__1);
                }
            for (i__ = 1; i__ <= *p; ++i__)
                v[i__ + l * v_dim1] = 0.;
            v[l + l * v_dim1] = 1.;
        }

    /*     main iteration loop for the singular values. */
    mm = m;
    iter = 0;
L360:
    // quit if all the singular values have been found.
    if (m == 0)
        return;
    if (iter < maxit) // if too many iterations have been performed, set  flag and return.
        goto L370;
    *info = m;
    return;
L370:

/*      This section of the program inspects for negligible elements in the s and
        e arrays.  On completion the variables kase and l are set as follows.

           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m 
           kase = 2     if s(l) is negligible and l.lt.m 
           kase = 3     if e(l-1) is negligible, l.lt.m, and 
                        s(l), ..., s(m) are not negligible (qr step). 
           kase = 4     if e(m-1) is negligible (convergence). */

    for (ll = 1; ll <= m; ++ll) {
        l = m - ll;
        if (l == 0)
            break;
        test = abs(s[l]) + abs(s[l + 1]);
        ztest = test + abs(e[l]);
        acc = abs(test - ztest) / (test + 1e-100);
        if (acc > 1e-15)
            continue;
    /*            if (ztest .ne. test) go to 380 */
        e[l] = 0.;
        break;
    }
    if (l != m - 1)
        goto L410;
    kase = 4;
    goto L480;
L410:
    lp1 = l + 1;
    mp1 = m + 1;
    for (lls = lp1; lls <= mp1; ++lls) {
        ls = m - lls + lp1;
        if (ls == l)
            break;
        test = 0.;
        if (ls != m)
            test += abs(e[ls]);
        if (ls != l + 1)
            test += abs(e[ls - 1]);
        ztest = test + abs(s[ls]);
        /* 1.0d-100 is to guard against a zero matrix, hence zero test */
        acc = abs(test - ztest) / (test + 1e-100);
        if (acc > 1e-15)
            continue;
        /*               if (ztest .ne. test) go to 420 */
        s[ls] = 0.;
        break;
    }
    if (ls != l)
        goto L450;
    kase = 3;
    goto L470;
L450:
    if (ls != m)
        goto L460;
    kase = 1;
    goto L470;
L460:
    kase = 2;
    l = ls;
L470:
L480:
    ++l;

    switch (kase) { /*        perform the task indicated by kase. */
        case 1:  goto L490;
        case 2:  goto L520;
        case 3:  goto L540;
        case 4:  goto L570;
    }

/*        deflate negligible s(m). */
L490:
    mm1 = m - 1;
    f = e[m - 1];
    e[m - 1] = 0.;
    for (kk = l; kk <= mm1; ++kk) {
        k = mm1 - kk + l;
        t1 = s[k];
        drotg_(&t1, &f, &cs, &sn);
        s[k] = t1;
        if (k != l) {
            f = -sn * e[k - 1];
            e[k - 1] = cs * e[k - 1];
        }
        if (wantv)
            drot_(p, &v[k * v_dim1 + 1], &c__1, &v[m * v_dim1 + 1], &c__1, & cs, &sn);
    }
    goto L610;

/*        split at negligible s(l). */
L520:
    f = e[l - 1];
    e[l - 1] = 0.;
    for (k = l; k <= m; ++k) {
        t1 = s[k];
        drotg_(&t1, &f, &cs, &sn);
        s[k] = t1;
        f = -sn * e[k];
        e[k] = cs * e[k];
        if (wantu)
            drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(l - 1) * u_dim1 + 1], & c__1, &cs, &sn);
    }
    goto L610;

/*        perform one qr step. */
L540:

    /*           calculate the shift. */
    scale = GSL_MAX(abs(s[m]),  //chain binary maxes to find largest in the group.
                GSL_MAX(abs(s[m - 1]),
                GSL_MAX(abs(e[m - 1]), 
                GSL_MAX(abs(s[l]), 
                GSL_MAX(abs(e[l]), abs(s[m - 1]))))));
    sm = s[m] / scale;
    smm1 = s[m - 1] / scale;
    emm1 = e[m - 1] / scale;
    sl = s[l] / scale;
    el = e[l] / scale;
    b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.;
    c__ = gsl_pow_2( sm * emm1);
    shift = 0.;
    if (b == 0. && c__ == 0.)
        goto L550;
    shift = sqrt(b * b + c__);
    if (b < 0.)
        shift = -shift;
    shift = c__ / (b + shift);
L550:
    f = (sl + sm) * (sl - sm) + shift;
    g = sl * el;

    /*           chase zeros. */
    mm1 = m - 1;
    for (k = l; k <= mm1; ++k) {
        drotg_(&f, &g, &cs, &sn);
        if (k != l)
            e[k - 1] = f;
        f = cs * s[k] + sn * e[k];
        e[k] = cs * e[k] - sn * s[k];
        g = sn * s[k + 1];
        s[k + 1] = cs * s[k + 1];
        if (wantv)
            drot_(p, &v[k * v_dim1 + 1], &c__1, &v[(k + 1) * v_dim1 + 1], & c__1, &cs, &sn);
        drotg_(&f, &g, &cs, &sn);
        s[k] = f;
        f = cs * e[k] + sn * s[k + 1];
        s[k + 1] = -sn * e[k] + cs * s[k + 1];
        g = sn * e[k + 1];
        e[k + 1] = cs * e[k + 1];
        if (wantu && k < *n)
            drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(k + 1) * u_dim1 + 1], & c__1, &cs, &sn);
    }
    e[m - 1] = f;
    ++iter;
    goto L610;

/*        convergence. */
L570:

/*           make the singular value  positive. */

    if (s[l] < 0.) {
        s[l] = -s[l];
        if (wantv)
            dscal(p, -1., &v[l * v_dim1 + 1]);
    }
/*           order the singular value. */

    while ((l != mm) && (s[l] < s[l + 1])) {
        t = s[l];
        s[l] = s[l + 1];
        s[l + 1] = t;
        if (wantv && l < *p)
            dswap(*p, &v[l * v_dim1 + 1], &v[(l + 1) * v_dim1 + 1]);
        if (wantu && l < *n)
            dswap(*n, &u[l * u_dim1 + 1], &u[(l + 1) * u_dim1 + 1]);
        ++l;
    }
    iter = 0;
    --m;
L610:
    goto L360;
} /* dsvdc_ */

////// loessc.c
#define	GAUSSIAN	1
#define SYMMETRIC	0

static long	*iv, liv, lv, tau;
static double *v;

/* begin ehg's FORTRAN-callable C-codes */

static void loess_error(int i){ //used to be ehg182.
  char *mess, mess2[50];
    switch(i){
case 101: mess="d>dMAX in ehg131.  Need to recompile with increased dimensions."; break;
case 102: mess="liv too small.   (Discovered by lowesd)"; break;
case 103: mess="lv too small.    (Discovered by lowesd)"; break;
case 104: mess="span too small.  fewer data values than degrees of freedom."; break;
case 105: mess="k>d2MAX in ehg136.  Need to recompile with increased dimensions."; break;
case 106: mess="lwork too small"; break;
case 110: mess="not enough extra workspace for robustness calculation"; break;
case 120: mess="zero-width neighborhood. make span bigger"; break;
case 121: mess="all data on boundary of neighborhood. make span bigger"; break;
case 123: mess="ihat=1 (diag L) in l2fit only makes sense if z=x (eval=data)."; break;
case 171: mess="lowesd must be called first."; break;
case 172: mess="lowesf must not come between lowesb and lowese, lowesr, or lowesl."; break;
case 173: mess="lowesb must come before lowese, lowesr, or lowesl."; break;
case 174: mess="lowesb need not be called twice."; break;
case 175: mess="need setLf=.true. for lowesl."; break;
case 180: mess="nv>nvmax in cpvert."; break;
case 182: mess="svddc failed in l2fit."; break;
case 185: mess="trouble descending to leaf in vleaf."; break;
case 186: mess="insufficient workspace for lowesf."; break;
case 187: mess="insufficient stack space"; break;
case 193: mess="workspace in loread appears to be corrupted"; break;
case 194: mess="trouble in l2fit/l2tr"; break;
case 195: mess="only constant, linear, or quadratic local models allowed"; break;
case 196: mess="degree must be at least 1 for vertex influence matrix"; break;
default: sprintf(mess=mess2,"Assert failed; error code %d\n",i); break;
    }
    Apop_assert_n(0, "%s", mess);
}

static void ehg183_(char *s, integer *i, integer n, integer inc) {
  char mess[4000], num[20];
  strcpy(mess,s);
  for (int j=0; j<n; j++) {
    sprintf(num," %ld",i[j * inc]);
    strcat(mess,num);
  }
  Warning(mess);
}

static void ehg184_(char *s, double *x, integer n, integer inc) {
  char mess[4000], num[30];
  strcpy(mess,s);
  for (int j=0; j< n; j++) {
    sprintf(num," %.5g",x[j * inc]);
    strcat(mess,num);
  }
  Warning(mess);
}

////// loessf.f -- translated by f2c (version 20061008). 

static void dqrdc_(double *x, integer *ldx, integer *n, integer *p, double *qraux, 
        integer *jpvt, double *work, integer job) ;
static double ehg128_(double *, integer *, integer *, integer *, integer *, 
        double *, integer *, integer *, integer *, double *, integer *, double *);
static void ehg129_(integer *, integer *, integer *, double *, integer *, integer, double *),
ehg139_(double *, integer *, integer *, integer *, integer *, integer *, double *, double *,
        integer *, integer *, double *, double *, double *, integer *, integer *,
	    double *, double *, double *, double *, integer *, double *, double *,
        double *, integer *, integer *, integer *, double *, integer *, integer *, integer *,
        integer *, double *, integer *, integer *, integer *, integer *, integer *, double *,
        logical *, double *),
ehg197(integer deg, integer d__, double f, integer *dk, double *trl);
static double ehg176_(double *);

static integer ifloor(double x) {
    integer ret_val = x;
    if ((double) ret_val > x)
        --ret_val;
    return ret_val;
}

static void ehg126_(integer *d__, integer *n, integer *vc, double *x, double *v, integer *nvmax) {
    integer v_dim1, x_dim1;
    static integer execnt = 0, i__, j, k;
    static double t, mu, beta, alpha, machin;

    x_dim1 = *n;
    x -= 1 + x_dim1;
    v_dim1 = *nvmax;
    v -= 1 + v_dim1;

    ++execnt;
    if (execnt == 1)
        machin = DBL_MAX;
/*     fill in vertices for bounding box of $x$ */
/*     lower left, upper right */
    for (k = 1; k <= *d__; ++k) {
        alpha = machin;
        beta = -machin;
        for (i__ = 1; i__ <= *n; ++i__) {
            t = x[i__ + k * x_dim1];
            alpha = min(alpha,t);
            beta = max(beta,t);
        }
    /*        expand the box a little */
        mu = .005 * GSL_MAX(beta - alpha,
                            GSL_MAX(abs(alpha), abs(beta)) * 1e-10 + 1e-30);
        alpha -= mu;
        beta += mu;
        v[k * v_dim1 + 1] = alpha;
        v[*vc + k * v_dim1] = beta;
    }
/*     remaining vertices */
    for (i__ = 2; i__ <= *vc - 1; ++i__) {
        j = i__ - 1;
        for (k = 1; k <= *d__; ++k) {
            v[i__ + k * v_dim1] = v[j % 2 * (*vc - 1) + 1 + k * v_dim1];
            j = (integer) ((double) j / 2.);
        }
    }
} /* ehg126_ */

static void ehg125_(integer *p, integer *nv, double *v, integer *vhit, integer nvmax, integer d__,
       integer k, double *t, integer *r__, integer *s, integer *f, integer *l, integer *u) {

    integer f_dim1, l_dim1, u_dim1, v_dim1;
    static integer h__, i__, j, m, i3, mm, execnt = 0;
    static logical match;

    --vhit;
    v_dim1 = nvmax;
    v -= 1 + v_dim1;
    u_dim1 = *r__;
    u -= 1 + (u_dim1 << 1);
    l_dim1 = *r__;
    l -= 1 + (l_dim1 << 1);
    f_dim1 = *r__;
    f -= 1 + (f_dim1 << 1);

    ++execnt;
    h__ = *nv;
    for (i__ = 1; i__ <= *r__; ++i__)
        for (j = 1; j <= *s; ++j) {
            ++h__;
            for (i3 = 1; i3 <= d__; ++i3)
                v[h__ + i3 * v_dim1] = v[f[i__ + (j << 1) * f_dim1] + i3 * v_dim1];
            v[h__ + k * v_dim1] = *t;
    /*           check for redundant vertex */
            match = FALSE_;
            m = 1;
            while (!match && (m <=*nv)){
                match = v[m + v_dim1] == v[h__ + v_dim1];
                for (mm = 2 ;match && ( mm <= d__); mm++)
                    match = v[m + mm * v_dim1] == v[h__ + mm * v_dim1];
                ++m;
                }
            --m;
            if (match)
                --h__;
            else {
                m = h__;
                if (vhit[1] >= 0)
                    vhit[m] = *p;
            }
            l[i__ + (j << 1) * l_dim1] = f[i__ + (j << 1) * f_dim1];
            l[i__ + ((j << 1) + 1) * l_dim1] = m;
            u[i__ + (j << 1) * u_dim1] = m;
            u[i__ + ((j << 1) + 1) * u_dim1] = f[i__ + ((j << 1) + 1) * f_dim1];
        }
    *nv = h__;
    if (! (*nv <= nvmax))
        loess_error(180);
} /* ehg125_ */

static void find_kth_smallest(integer il, integer ir, integer k, integer nk, double *p, integer *pi) {
    //Formerly ehg106
    integer p_dim1;
    static integer execnt = 0, i__, j, l, r, ii;
    static double t;

    --pi;
    p_dim1 = nk;
    p -= 1 + p_dim1;

    ++execnt;
    /*     find the $k$-th smallest of $n$ elements */
    /*     Floyd+Rivest, CACM Mar '75, Algorithm 489 */
    l = il;
    r = ir;
    while (l < r) {
        /*  to avoid recursion, sophisticated partition deleted */
        /*  partition $x sub {l..r}$ about $t$ */
        t = p[pi[k] * p_dim1 + 1];
        i__ = l;
        j = r;
        ii = pi[l];
        pi[l] = pi[k];
        pi[k] = ii;
        if (t < p[pi[r] * p_dim1 + 1]) {
            ii = pi[l];
            pi[l] = pi[r];
            pi[r] = ii;
        }
        while (i__ < j) {
            ii = pi[i__];
            pi[i__] = pi[j];
            pi[j] = ii;
            ++i__;
            --j;
            while (p[pi[i__] * p_dim1 + 1] < t)
                ++i__;
            while (t < p[pi[j] * p_dim1 + 1])
                --j;
        }
        if (p[pi[l] * p_dim1 + 1] == t) {
            ii = pi[l];
            pi[l] = pi[j];
            pi[j] = ii;
        } else {
            ++j;
            ii = pi[r];
            pi[r] = pi[j];
            pi[j] = ii;
        }
        if (j <= k)
            l = j + 1;
        if (k <= j)
            r = j - 1;
    }
} /* ehg106_ */

static integer idamax(integer n, double *dx, integer incx) {
    /*     Finds the index of element having max. absolute value. */
    /*     jack dongarra, linpack, 3/11/78. */
    int ret_val = 1;
    static integer i__, ix;
    static double dmax__;
    --dx;

    if (n < 1)
        return 0;
    if (n == 1)
        return 1;
    if (incx != 1) { // code for increment not equal to 1
        ix = 1;
        dmax__ = abs(dx[1]);
        ix += incx;
        for (i__ = 2; i__ <= n; ++i__, ix += incx)
            if ( abs(dx[ix]) > dmax__){
                ret_val = i__;
                dmax__ = abs(dx[ix]);
            }
    } else { // code for increment equal to 1
        dmax__ = abs(dx[1]);
        for (i__ = 2; i__ <= n; ++i__)
            if ( abs(dx[i__]) > dmax__){
                ret_val = i__;
                dmax__ = abs(dx[i__]);
            }
    }
    return ret_val;
} /* idamax */

static void ehg124_(integer *ll, integer *uu, integer d__, integer n, integer *nv, integer *nc,
        integer *ncmax, integer *vc, double * x, integer *pi, integer *a, double *xi,
        integer *lo, integer *hi, integer *c__, double *v, integer *vhit, integer nvmax, integer *
        fc, double *fd, integer *dd) {
    integer c_dim1, v_dim1, v_offset, x_dim1, x_offset, i__1, i__3;
    static integer execnt = 0, k, l, m, p, u, i4, check, lower, upper, inorm2, offset;
    static logical i1, i2, leaf;
    static double diag[8], diam, sigma[8];

    --pi; --hi; --lo; --xi; --a; --vhit;
    x_dim1 = n;
    x -= x_offset = 1 + x_dim1;
    c_dim1 = *vc;
    c__ -= 1 + c_dim1;
    v_dim1 = nvmax;
    v -= v_offset = 1 + v_dim1;

    ++execnt;
    p = 1;
    l = *ll;
    u = *uu;
    lo[p] = l;
    hi[p] = u;
    while (p <= *nc){
        for (i4 = 1; i4 <= *dd; ++i4)
            diag[i4 - 1] = v[c__[*vc + p * c_dim1] + i4 * v_dim1] - v[c__[p * c_dim1 + 1] + i4 * v_dim1];
        diam = 0.;
        for (inorm2 = 1; inorm2 <= *dd; ++inorm2)
            diam += gsl_pow_2(diag[inorm2 - 1]);
        diam = sqrt(diam);
        i1 = (u - l + 1 <= *fc)
             ? TRUE_
             : (diam <= *fd);
        if (i1)
            leaf = TRUE_;
        else {
            if (*ncmax < *nc + 2)
                i2 = TRUE_;
            else
                i2 = (double) (nvmax) < *nv + (double) (*vc) / 2.;
            leaf = i2;
        }
        if (! leaf) {
            ehg129_(&l, &u, dd, &x[x_offset], &pi[1], n, sigma);
            k = idamax(*dd, sigma, 1);
            m = (integer) ((double) (l + u) / 2.);
            find_kth_smallest(l, u, m, 1, &x[k * x_dim1 + 1], &pi[1]);
            /* bug fix from btyner@gmail.com 2006-07-20 */
            offset = 0;
            while (! (m + offset >= u || m + offset < l)) {
                if (offset < 0) {
                    lower = l;
                    check = m + offset;
                    upper = check;
                } else {
                    lower = m + offset + 1;
                    check = lower;
                    upper = u;
                }
                find_kth_smallest(lower, upper, check, 1, &x[k * x_dim1 + 1], &pi[1]);
                if (x[pi[m + offset] + k * x_dim1] == x[pi[m + offset + 1] + k * x_dim1]) {
                    offset = -offset;
                    if (offset >= 0)
                        ++offset;
                } else {
                    m += offset;
                    break;
                }
            }
            if (v[c__[p * c_dim1 + 1] + k * v_dim1] == x[pi[m] + k * x_dim1])
                leaf = TRUE_;
            else
                leaf = v[c__[*vc + p * c_dim1] + k * v_dim1] == x[pi[m] + k * x_dim1];
        }
        if (leaf)
            a[p] = 0;
        else {
            a[p] = k;
            xi[p] = x[pi[m] + k * x_dim1];
            /*           left son */
            ++(*nc);
            lo[p] = *nc;
            lo[*nc] = l;
            hi[*nc] = m;
            /*           right son */
            ++(*nc);
            hi[p] = *nc;
            lo[*nc] = m + 1;
            hi[*nc] = u;
            i__1 = pow_ii(2, k - 1);
            i__3 = pow_ii(2, d__ - k);
            ehg125_(&p, nv, &v[v_offset], &vhit[1], nvmax, d__, k, &xi[p], &i__1,
                 &i__3, &c__[p * c_dim1 + 1], &c__[lo[p] * c_dim1 + 1], &c__[hi[p] * c_dim1 + 1]);
        }
        ++p;
        l = lo[p];
        u = hi[p];
    }
} /* ehg124_ */

static void ehg127_(double *q, integer *n, integer *d__, integer *nf, double *f, double *x,
        integer *psi, double *y, double *rw, integer *kernel, integer *k, double *dist,
        double *eta, double *b, integer *od, double *w, double *rcond, integer *sing,
        double *sigma, double *u, double *e, double *dgamma, double *qraux,
        double * work, double *tol, integer *dd, integer *tdeg, integer *cdeg, double *s) {

    integer b_dim1, x_dim1, b_offset;
    double d__1;
    static integer execnt = 0, i__, j, i3, i9, jj, info, jpvt, inorm2, column;
    static double g[15], i2, rho, scal, machep, colnor[15];

    --rw; --y; --psi;
    x_dim1 = *n;
    x -= 1 + x_dim1;
    --q;
    --w;
    --eta;
    b_dim1 = *nf;
    b -= b_offset = 1+b_dim1;
    --sigma;
    u -= 16;
    e -= 16;
    --dgamma; --qraux; --work; --cdeg;

    ++execnt;
    if (execnt == 1)
        machep = DBL_EPSILON;
    /*     sort by distance */
    for (i3 = 1; i3 <= *n; ++i3)
        dist[i3] = 0.;
    for (j = 1; j <= *dd; ++j)
        for (i3 = 1; i3 <= *n; ++i3)
            dist[i3] += gsl_pow_2(x[i3 + j * x_dim1] - q[j]);
    find_kth_smallest(1, *n, *nf, 1, &dist[1], &psi[1]);
    rho = dist[psi[*nf]] * max(1.,*f);
    if (! (0. < rho))
        loess_error(120);
    /*     compute neighborhood weights */
    if (*kernel == 2) {
        for (i__ = 1; i__ <= *nf; ++i__)
            w[i__] = (dist[psi[i__]] < rho)
                    ? sqrt(rw[psi[i__]])
                    : 0.;
    } else {
        for (i3 = 1; i3 <= *nf; ++i3)
            w[i3] = sqrt(dist[psi[i3]] / rho);
        for (i3 = 1; i3 <= *nf; ++i3)
            w[i3] = sqrt(rw[psi[i3]] * gsl_pow_3(1 - gsl_pow_3(w[i3])) );
    }
    if (abs(w[idamax(*nf, &w[1], 1)]) == 0.) { //why |x|==0, and not just x == 0 ?
        ehg184_("at ", &q[1], *dd, 1);
        ehg184_("radius ", &rho, 1,1);
        loess_error(121);
    }
    /*     fill design matrix */
    column = 1;
    for (i3 = 1; i3 <= *nf; ++i3)
        b[i3 + column * b_dim1] = w[i3];
    if (*tdeg >= 1)
        for (j = 1; j <= *d__; ++j)
            if (cdeg[j] >= 1) {
                ++column;
                for (i3 = 1; i3 <= *nf; ++i3)
                    b[i3 + column * b_dim1] = w[i3] * (x[psi[i3] + j * x_dim1] - q[j]);
            }
    if (*tdeg >= 2) {
        for (j = 1; j <= *d__; ++j)
            if (cdeg[j] >= 1) {
                if (cdeg[j] >= 2) {
                    ++column;
                    for (i3 = 1; i3 <= *nf; ++i3)
                        b[i3 + column * b_dim1] = w[i3] * gsl_pow_2(x[psi[i3] + j * x_dim1] - q[j]);
                }
                for (jj = j + 1; jj <= *d__; ++jj)
                    if (cdeg[jj] >= 1) {
                        ++column;
                        for (i3 = 1; i3 <= *nf; ++i3)
                            b[i3 + column * b_dim1] = w[i3] * (x[psi[i3] + j * x_dim1] - q[j])
                                                        * (x[psi[i3] + jj * x_dim1] - q[jj]);
                    }
            }
        *k = column;
    }
    for (i3 = 1; i3 <= *nf; ++i3)
        eta[i3] = w[i3] * y[psi[i3]];
    /*     equilibrate columns */
    for (j = 1; j <= *k; ++j) {
        scal = 0.;
        for (inorm2 = 1; inorm2 <= *nf; ++inorm2)
            scal += gsl_pow_2(b[inorm2 + j * b_dim1]);
        scal = sqrt(scal);
        if (0. < scal) {
            for (i3 = 1; i3 <= *nf; ++i3)
                b[i3 + j * b_dim1] /= scal;
            colnor[j - 1] = scal;
        } else
            colnor[j - 1] = 1.;
    }
/*     singular value decomposition */
    dqrdc_(&b[b_offset], nf, nf, k, &qraux[1], &jpvt, &work[1], 0);
    dqrsl_(&b[b_offset], nf, nf, k, &qraux[1], &eta[1], &work[1], &eta[1], & eta[1], &work[1], &work[1], 1000, &info);
    for (i9 = 1; i9 <= *k; ++i9)
        for (i3 = 1; i3 <= *k; ++i3)
            u[i3 + i9 * 15] = 0.;
    for (i__ = 1; i__ <= *k; ++i__)
        for (j = i__; j <= *k; ++j)
            u[i__ + j * 15] = b[i__ + j * b_dim1];
    dsvdc_(&u[16], &c__15, k, k, &sigma[1], g, &u[16], &c__15, &e[16], &c__15, &work[1], &c__21, &info);
    if (info != 0)
        loess_error(182);
    *tol = sigma[1] * (machep * 100);
    *rcond = GSL_MIN(*rcond, sigma[*k]/sigma[1]);
    if (sigma[*k] <= *tol) {
        ++(*sing);
        if (*sing == 1) {
            ehg184_("Warning. pseudoinverse used at", &q[1], *d__, 1);
            d__1 = sqrt(rho);
            ehg184_("neighborhood radius", &d__1, 1, 1);
            ehg184_("reciprocal condition number ", rcond, 1, 1);
        } else if (*sing == 2)
            ehg184_("There are other near singularities as well.", &rho, 1, 1);
    }
/*     compensate for equilibration */
    for (j = 1; j <= *k; ++j) 
        for (i3 = 1; i3 <= *k; ++i3)
            e[j + i3 * 15] /= colnor[j - 1];
/*     solve least squares problem */
    for (j = 1; j <= *k; ++j) {
        if (*tol < sigma[j])
            i2 = ddot_(k, &u[j * 15 + 1], &c__1, &eta[1], &c__1) / sigma[j];
        else
            i2 = 0.;
        dgamma[j] = i2;
    }
    for (j = 0; j <= *od; ++j)  //bug fix 2006-07-04 for k=1, od>1.   (thanks btyner@gmail.com) */
        if (j < *k)
            s[j] = ddot_(k, &e[j + 16], &c__15, &dgamma[1], &c__1);
        else
            s[j] = 0.;
} /* ehg127_ */

static void ehg129_(integer *l, integer *u, integer *d__, double *x, integer *pi, integer n, double *sigma) {
    integer x_dim1;
    static integer execnt = 0;
    static double t, beta, alpha, machin;
    --sigma;
    --pi;
    x_dim1 = n;
    x -= 1 + x_dim1;
    ++execnt;
    if (execnt == 1)
        machin = DBL_MAX;
    for (integer k = 1; k <= *d__; ++k) {
        alpha = machin;
        beta = -machin;
        for (integer i__ = *l; i__ <= *u; ++i__) {
            t = x[pi[i__] + k * x_dim1];
            alpha = GSL_MIN(alpha, x[pi[i__] + k * x_dim1]);
            beta = max(beta,t);
        }
        sigma[k] = beta - alpha;
    }
} /* ehg129_ */

static void ehg131_(double *x, double *y, double *rw, double *trl, 
        double *diagl, integer *kernel, integer *k, integer *n, integer *d__, 
        integer *nc, integer *ncmax, integer *vc, integer *nv, integer *nvmax, 
        integer *nf, double *f, integer *a, integer *c__, integer *hi, integer *lo, 
        integer *pi, integer *psi, double *v, integer *vhit, double *vval, double *xi,
        double *dist, double *eta, double *b, integer *ntol, double *fd, 
        double *w, double *vval2, double *rcond, integer *sing, integer *dd, 
        integer *tdeg, integer *cdeg, integer *lq, double *lf, logical *setlf) {

    integer lq_dim1, lq_offset, c_dim1, c_offset, lf_dim1, lf_dim2, lf_offset,
	     v_dim1, v_offset, vval_dim1, vval_offset, vval2_dim1, vval2_offset, x_dim1, x_offset;

    static integer execnt = 0, j, i1, i2;
    static double delta[8];
    static integer identi;

    --psi; --pi; 
    x_dim1 = *n;
    x -= x_offset = 1 + x_dim1;
    --xi;
    --lo;
    --hi;
    --a;
    c_dim1 = *vc;
    c__ -= c_offset = 1 + c_dim1;
    vval2_dim1 = *d__ - 0 + 1;
    vval2 -= vval2_offset = 0 + vval2_dim1;
    vval_dim1 = *d__ - 0 + 1;
    vval -= vval_offset = 0 + vval_dim1;
    --vhit;
    v_dim1 = *nvmax;
    v -= v_offset = 1 + v_dim1;
    lf_dim1 = *d__ - 0 + 1;
    lf_dim2 = *nvmax;
    lf -= lf_offset = 0 + lf_dim1 * (1 + lf_dim2);
    lq_dim1 = *nvmax;
    lq -= lq_offset = 1 + lq_dim1;
    --w; --eta; --b; --cdeg;

    ++execnt;
    if (! (*d__ <= 8))
        loess_error(101);
/*     build $k$-d tree */
    ehg126_(d__, n, vc, &x[x_offset], &v[v_offset], nvmax);
    *nv = *vc;
    *nc = 1;
    for (j = 1; j <= *vc; ++j) {
        c__[j + *nc * c_dim1] = j;
        vhit[j] = 0;
    }
    for (i1 = 1; i1 <= *d__; ++i1)
        delta[i1 - 1] = v[*vc + i1 * v_dim1] - v[i1 * v_dim1 + 1];
    *fd *= dnrm2(d__, delta);
    for (identi = 1; identi <= *n; ++identi)
        pi[identi] = identi;
    ehg124_(&c__1, n, *d__, *n, nv, nc, ncmax, vc, &x[x_offset], &pi[1], &a[1],
	    &xi[1], &lo[1], &hi[1], &c__[c_offset], &v[v_offset], &vhit[1], *nvmax, ntol, fd, dd);
    // Smooth
    if (*trl != 0.)
        for (i2 = 1; i2 <= *nv; ++i2)
            for (i1 = 0; i1 <= *d__; ++i1)
                vval2[i1 + i2 * vval2_dim1] = 0.;
    ehg139_(&v[v_offset], nvmax, nv, n, d__, nf, f, &x[x_offset], &pi[1], &psi[1], y,
            rw, trl, kernel, k, dist, dist, &eta[1], &b[1], d__, &w[1], diagl,
            &vval2[vval2_offset], nc, vc, &a[1], &xi[1], &lo[1], &hi[1], &c__[c_offset], &vhit[1], rcond, sing,
            dd, tdeg, &cdeg[1], &lq[lq_offset], &lf[lf_offset], setlf, &vval[vval_offset]);
} /* ehg131_ */

static void ehg133_(integer *n, integer *d__, integer *vc, integer *nvmax, integer *nc, 
        integer *ncmax, integer *a, integer *c__, integer *hi, integer *lo, double *v, 
        double *vval, double *xi, integer m, double *z__, double *s) {
    integer c_dim1, c_offset, v_dim1, v_offset, vval_dim1, vval_offset, z_dim1, z_offset;

    static integer execnt = 0, i__, i1;
    static double delta[8];

    vval_dim1 = *d__ - 0 + 1;
    vval -= vval_offset = 0 + vval_dim1;
    --s;
    v_dim1 = *nvmax;
    v -= v_offset = 1 + v_dim1;
    c_dim1 = *vc;
    c__ -= c_offset = 1 + c_dim1;
    z_dim1 = m;
    z__ -= z_offset = 1 + z_dim1;

    ++execnt;
    for (i__ = 1; i__ <= m; ++i__) {
        for (i1 = 1; i1 <= *d__; ++i1)
            delta[i1 - 1] = z__[i__ + i1 * z_dim1];
        s[i__] = ehg128_(delta, d__, ncmax, vc, a, xi, lo, hi,
             &c__[c_offset], &v[v_offset], nvmax, &vval[vval_offset]);
    }
}

static void set_cs(integer d, integer i, double *c1, double *c2, double *c3){
    static double c[48] = { .297162,.380266,.5886043,.4263766,.3346498,
	    .6271053,.5241198,.3484836,.6687687,.6338795,.4076457,.7207693,
	    .1611761,.3091323,.4401023,.2939609,.3580278,.5555741,.397239,
	    .4171278,.6293196,.4675173,.469907,.6674802,.2848308,.2254512,
	    .2914126,.5393624,.251723,.389897,.7603231,.2969113,.474013,
	    .9664956,.3629838,.5348889,.207567,.2822574,.2369957,.3911566,
	    .2981154,.3623232,.5508869,.3501989,.4371032,.7002667,.4291632,
	    .493037 };
    if (d <= 4) {
        *c1 = c[i - 1];
        *c2 = c[i];
        *c3 = c[i + 1];
    } else {
        *c1 = c[i - 1] + (d - 4) * (c[i - 1] - c[i - 4]);
        *c2 = c[i]     + (d - 4) * (c[i]     - c[i - 3]);
        *c3 = c[i + 1] + (d - 4) * (c[i + 1] - c[i - 2]);
    }
}

static void ehg141_(double *trl, integer *n, integer *deg, integer *k, integer *d,
        integer *nsing, integer *dk, double * delta1, double *delta2) {

    static integer i;
    static double z, c1, c2, c3, c4, corx;

/*     coef, d, deg, del */
    if (*deg == 0)
        *dk = 1;
    if (*deg == 1)
        *dk = *d + 1;
    if (*deg == 2)
        *dk = ((double) ((*d + 2) * (*d + 1)) / 2.);
    corx = sqrt(*k / (double) (*n));
    z = (sqrt(*k / *trl) - corx) / (1 - corx);
    if ((*nsing == 0 && 1. < z) ||(z < 0.) )
        ehg184_("Chernobyl! trL<k", trl, 1, 1);
    z = GSL_MIN(1.,  GSL_MAX(0.,z));
    c4 = exp(ehg176_(&z));
    i = (min(*d,4) - 1 + ((*deg - 1) << 2)) * 3 + 1;
    set_cs(*d, i, &c1, &c2, &c3);
    *delta1 = *n - *trl * exp(c1 * pow(z, c2) * pow(1-z, c3) * c4);
    i += 24;
    set_cs(*d, i, &c1, &c2, &c3);
    *delta2 = *n - *trl * exp(c1 * pow(z, c2) * pow(1-z, c3) * c4);
} /* ehg141_ */

static void lowesc_(integer *n, double *l, double *ll, double *trl, double *delta1, double *delta2) {
    static integer execnt = 0, i__, j;
    integer l_dim1, ll_dim1;

    ll_dim1 = *n;
    ll -= 1 + ll_dim1;
    l_dim1 = *n;
    l -= 1 + l_dim1;

    ++execnt;
/*     compute $LL~=~(I-L)(I-L)'$ */
    for (i__ = 1; i__ <= *n; ++i__)
        --l[i__ + i__ * l_dim1];
    for (i__ = 1; i__ <= *n; ++i__)
        for (j = 1; j <= i__; ++j)
            ll[i__ + j * ll_dim1] = ddot_(n, &l[i__ + l_dim1], n, &l[j + l_dim1], n);
    for (i__ = 1; i__ <= *n; ++i__)
        for (j = i__ + 1; j <= *n; ++j)
            ll[i__ + j * ll_dim1] = ll[j + i__ * ll_dim1];
    for (i__ = 1; i__ <= *n; ++i__)
        ++l[i__ + i__ * l_dim1];
/*     accumulate first two traces */
    *trl = 0.;
    *delta1 = 0.;
    for (i__ = 1; i__ <= *n; ++i__) {
        *trl += l[i__ + i__ * l_dim1];
        *delta1 += ll[i__ + i__ * ll_dim1];
    }
/*     $delta sub 2 = "tr" LL sup 2$ */
    *delta2 = 0.;
    for (i__ = 1; i__ <= *n; ++i__)
        *delta2 += ddot_(n, &ll[i__ + ll_dim1], n, &ll[i__ * ll_dim1 + 1], & c__1);
} /* lowesc_ */

static void ehg169_(integer d__, integer *vc, integer *nc, integer *ncmax, integer *nv, 
        integer nvmax, double *v, integer *a, double *xi, integer *c__, integer *hi, integer *lo) {
    integer c_dim1, v_dim1, v_offset, i__1, i__3;
    static integer execnt = 0, i__, j, k, p, mc, mv, novhit[1];

    --lo;
    --hi;
    c_dim1 = *vc;
    c__ -= 1 + c_dim1;
    --xi;
    --a;
    v_dim1 = nvmax;
    v -= v_offset = 1 + v_dim1;

    ++execnt;
    /*     as in bbox */
    /*     remaining vertices */
    for (i__ = 2; i__ <= *vc - 1; ++i__) {
        j = i__ - 1;
        for (k = 1; k <= d__; ++k) {
            v[i__ + k * v_dim1] = v[j % 2 * (*vc - 1) + 1 + k * v_dim1];
            j = ifloor((double)j / 2.);
        }
    }
    /*     as in ehg131 */
    mc = 1;
    mv = *vc;
    novhit[0] = -1;
    for (j = 1; j <= *vc; ++j)
        c__[j + mc * c_dim1] = j;
    /*     as in rbuild */
    for (p=1; p <= *nc; p++)
        if (a[p] != 0) {
            k = a[p];
            // left son
            ++mc;
            lo[p] = mc;
            // right son
            ++mc;
            hi[p] = mc;
            i__1 = pow_ii(2, k-1);
            i__3 = pow_ii(2, d__ - k);
            ehg125_(&p, &mv, &v[v_offset], novhit, nvmax, d__, k, &xi[p], &i__1,
                &i__3, &c__[p * c_dim1 + 1], &c__[lo[p] * c_dim1 + 1], &c__[hi[p] * c_dim1 + 1]);
        }
    if (! (mc == *nc) || ! (mv == *nv))
        loess_error(193);
}

static double ehg176_(double *z) {
    static integer d__ = 1;
    static integer vc = 2;
    static integer nv = 10;
    static integer nc = 17;
    static integer a[17] = { 1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0 };
    static struct {
        integer e_1[7];
        integer fill_2[7];
        integer e_3;
        integer fill_4[2];
	} equiv_94 = { {3, 5, 7, 9, 11, 13, 15}, {0}, 17 };//seven numbers, seven zeros, a number.

#define hi ((integer *)&equiv_94)

    static struct {
        integer e_1[7];
        integer fill_2[7];
        integer e_3;
        integer fill_4[2];
	} equiv_95 = { {2, 4, 6, 8, 10, 12, 14}, {0}, 16 };

#define lo ((integer *)&equiv_95)

    static struct {
        double e_1[7];
        double fill_2[7];
        double e_3;
        double fill_4[2];
	} equiv_96 = { {.3705, .2017, .5591, .1204, .2815, .4536, .7132}, {0}, .8751 };

#define xi ((double *)&equiv_96)

    static integer c__[34]	/* was [2][17] */ = { 1,2,1,3,3,2,1,4,4,3,3,5,
	    5,2,1,6,6,4,4,7,7,3,3,8,8,5,5,9,9,2,9,10,10,2 };
    static double vval[20]	/* was [2][10] */ = { -.090572,4.4844,
	    -.010856,-.7736,-.053718,-.3495,.026152,-.7286,-.058387,.1611,
	    .095807,-.7978,-.031926,-.4457,-.06417,.032813,-.020636,.335,
	    .040172,-.041032 };
    static double v[10]	/* was [10][1] */ = { -.005,1.005,.3705,.2017,
	    .5591,.1204,.2815,.4536,.7132,.8751 };

    return ehg128_(z, &d__, &nc, &vc, a, xi, lo, hi, c__, v, &nv, vval);
}

#undef xi
#undef lo
#undef hi

static void lowesa_(double *trl, integer *n, integer *d__,
            integer *tau, integer *nsing, double *delta1, double *delta2) {
    static integer execnt = 0, dka, dkb;
    static double d1a, d1b, d2a, d2b, alpha;

    ++execnt;
    ehg141_(trl, n, &c__1, tau, d__, nsing, &dka, &d1a, &d2a);
    ehg141_(trl, n, &c__2, tau, d__, nsing, &dkb, &d1b, &d2b);
    alpha = (double) (*tau - dka) / (double) (dkb - dka);
    *delta1 = (1 - alpha) * d1a + alpha * d1b;
    *delta2 = (1 - alpha) * d2a + alpha * d2b;
} /* lowesa_ */

static void ehg191_(integer *m, double *z__, double *l, integer *d__, integer *n, integer *nf,
        integer *nv, integer *ncmax, integer *vc, integer *a, double *xi, integer *lo, integer *hi,
        integer *c__, double *v, integer *nvmax, double *vval2, double *lf, integer *lq) {

    integer lq_dim1, c_offset, l_dim1, lf_dim1, lf_dim2, v_offset, vval2_dim1, vval2_offset, z_dim1;

    static integer execnt = 0, i__, j, p, i1, i2, lq1;
    static double zi[8];
    z_dim1 = *m;
    z__ -= 1 + z_dim1;
    l_dim1 = *m;
    l -=  1 + l_dim1;
    --hi;
    --lo;
    --xi;
    --a;
    c__ -= c_offset = 1 + *vc;
    lq_dim1 = *nvmax;
    lq -= 1 + lq_dim1;
    lf_dim1 = *d__ - 0 + 1;
    lf_dim2 = *nvmax;
    lf -= 0 + lf_dim1 * (1 + lf_dim2);
    vval2_dim1 = *d__ - 0 + 1;
    vval2 -= vval2_offset = 0 + vval2_dim1;
    v -= v_offset = 1 + *nvmax;

    ++execnt;
    for (j = 1; j <= *n; ++j) {
        for (i2 = 1; i2 <= *nv; ++i2)
            for (i1 = 0; i1 <= *d__; ++i1)
                vval2[i1 + i2 * vval2_dim1] = 0.;
        for (i__ = 1; i__ <= *nv; ++i__) { // linear search for i in Lq
            lq1 = lq[i__ + lq_dim1];
            lq[i__ + lq_dim1] = j;
            p = *nf;
            while (lq[i__ + p * lq_dim1] != j)
                --p;
            lq[i__ + lq_dim1] = lq1;
            if (lq[i__ + p * lq_dim1] == j) //BK: doesn't this always hit?  
                for (i1 = 0; i1 <= *d__; ++i1)
                    vval2[i1 + i__ * vval2_dim1] = lf[i1 + (i__ + p * lf_dim2) * lf_dim1];
        }
        for (i__ = 1; i__ <= *m; ++i__) {
            for (i1 = 1; i1 <= *d__; ++i1)
                zi[i1 - 1] = z__[i__ + i1 * z_dim1];
            l[i__ + j * l_dim1] = ehg128_(zi, d__, ncmax, vc, &a[1], &xi[1], &
                lo[1], &hi[1], &c__[c_offset], &v[v_offset], nvmax, & vval2[vval2_offset]);
        }
    }
} /* ehg191_ */

static void ehg196_(integer tau, integer d__, double f, double *trl) {
    static integer execnt = 0, dka, dkb;
    static double trla, trlb, alpha;

    ++execnt;
    ehg197(1, d__, f, &dka, &trla);
    ehg197(2, d__, f, &dkb, &trlb);
    alpha = (double) (tau - dka) / (double) (dkb - dka);
    *trl = (1 - alpha) * trla + alpha * trlb;
}

static void ehg197(integer deg, integer d__, double f, integer *dk, double *trl) {
    *dk = 0;
    if (deg == 1)
        *dk = d__ + 1;
    if (deg == 2)
        *dk = (integer) ((double) ((d__ + 2) * (d__ + 1)) / 2.);
    float g1 = (d__ * -.08125 + .13) * d__ + 1.05;
    *trl = *dk * (GSL_MAX( 0., (g1 - f) / f) + 1);
}

void hermite_prep(double h, double *phi0, double *phi1, double *psi0, double *psi1){
    *phi0 = (1-h)*(1-h)* (h * 2 + 1);
    *phi1 = h *h * (3 - h * 2);
    *psi0 = h * (1-h)*(1-h);
    *psi1 = h * h * (h - 1);
}

int xibar_search(const integer *a, const integer t[], const double* xi, const double xibar, integer nt){
    int m = nt -1;
    int done;
    while (1){
        if (m == 0)
            done = TRUE_;
        else {
            if (a[t[m - 1]] == 2)
                done = xi[t[m - 1]] == xibar;
            else
                done = FALSE_;
        }
        if (done)
            return m;
        --m;
    }
}

static double ehg128_(double *z__, integer *d__, integer *ncmax, integer *vc,
        integer *a, double *xi, integer *lo, integer *hi, integer *c__,
        double *v, integer *nvmax, double *vval) {
    integer c_dim1, v_dim1, vval_dim1;
    static double g[2304]	/* was [9][256] */, h__;
    static logical i2;
    static integer execnt = 0, t[20], i__, j, m, i1, i11, i12, ig, ii, lg, ll, nt, ur;
    static double g0[9], g1[9], s, v0, v1, ge, gn, gs, gw;
    static double gpe, gpn, gps, gpw, sew, sns, phi0, phi1, psi0, psi1, xibar;

    --z__; --hi; --lo; --xi; --a;
    c_dim1 = *vc;
    c__ -= 1 + c_dim1;
    vval_dim1 = *d__ - 0 + 1;
    vval -= 0 + vval_dim1;
    v_dim1 = *nvmax;
    v -= 1 + v_dim1;

    ++execnt;
    /*     locate enclosing cell */
    nt = 1;
    t[nt - 1] = 1;
    j = 1;
    while (a[j] != 0) {
        ++nt;
        /*     bug fix 2006-07-18 (thanks, btyner@gmail.com) */
        if (z__[a[j]] <= xi[j])
            i1 = lo[j];
        else
            i1 = hi[j];
        t[nt - 1] = i1;
        Apop_assert(nt < 20, "nt>=20 in eval.");
        j = t[nt - 1];
    }
    /*     tensor */
    for (i12 = 1; i12 <= *vc; ++i12)
        for (i11 = 0; i11 <= *d__; ++i11)
                g[i11 + i12 * 9 - 9] = vval[i11 + c__[i12 + j * c_dim1] * vval_dim1];
    lg = *vc;
    ll = c__[j * c_dim1 + 1];
    ur = c__[*vc + j * c_dim1];
    for (i__ = *d__; i__ >= 1; --i__) {
        h__ = (z__[i__] - v[ll + i__ * v_dim1]) / (v[ur + i__ * v_dim1] - v[ ll + i__ * v_dim1]);
        if (h__ < -.001) {
            ehg184_("eval ", &z__[1], *d__, 1);
            ehg184_("lowerlimit ", &v[ll + v_dim1], *d__, *nvmax);
        } else if (1.001 < h__) {
                ehg184_("eval ", &z__[1], *d__, 1);
                ehg184_("upperlimit ", &v[ur + v_dim1], *d__, *nvmax);
            }
        if (-.001 <= h__)
            i2 = h__ <= 1.001;
        else
            i2 = FALSE_;
        Apop_assert(i2, "extrapolation not allowed with blending.");
        lg = (integer) ((double) lg / 2.);
        for (ig = 1; ig <= lg; ++ig) {
            // Hermite basis
            hermite_prep(h__, &phi0, &phi1, &psi0, &psi1);
            g[ig * 9 - 9] = phi0 * g[ig * 9 - 9] + phi1 * g[(ig + lg) * 9 - 9]
                 + (psi0 * g[i__ + ig * 9 - 9] + psi1 * g[i__ + (ig + lg)
                * 9 - 9]) * (v[ur + i__ * v_dim1] - v[ll + i__ * v_dim1]);
            for (ii = 1; ii <= i__ - 1; ++ii)
                g[ii + ig * 9 - 9] = phi0 * g[ii + ig * 9 - 9] + phi1 * g[ii + (ig + lg) * 9 - 9];
        }
    }
    s = g[0];
/*     blending */
    if (*d__ == 2) {
    /*        ----- North ----- */
        v0 = v[ll + v_dim1];
        v1 = v[ur + v_dim1];
        for (i11 = 0; i11 <= *d__; ++i11) 
            g0[i11] = vval[i11 + c__[j * c_dim1 + 3] * vval_dim1];
        for (i11 = 0; i11 <= *d__; ++i11) 
            g1[i11] = vval[i11 + c__[j * c_dim1 + 4] * vval_dim1];
        xibar = v[ur + (v_dim1 << 1)];
        m= xibar_search(a, t, xi, xibar, nt);
        if (m >= 1) {
            m = hi[t[m - 1]];
            while (a[m] != 0)
                if (z__[a[m]] <= xi[m])
                    m = lo[m];
                else
                    m = hi[m];
            if (v0 < v[c__[m * c_dim1 + 1] + v_dim1]) {
                v0 = v[c__[m * c_dim1 + 1] + v_dim1];
                for (i11 = 0; i11 <= *d__; ++i11)
                    g0[i11] = vval[i11 + c__[m * c_dim1 + 1] * vval_dim1];
            }
            if (v[c__[m * c_dim1 + 2] + v_dim1] < v1) {
                v1 = v[c__[m * c_dim1 + 2] + v_dim1];
                for (i11 = 0; i11 <= *d__; ++i11)
                    g1[i11] = vval[i11 + c__[m * c_dim1 + 2] * vval_dim1];
            }
        }
        h__ = (z__[1] - v0) / (v1 - v0);
    /*        Hermite basis */
        hermite_prep(h__, &phi0, &phi1, &psi0, &psi1);
        gn = phi0 * g0[0] + phi1 * g1[0] + (psi0 * g0[1] + psi1 * g1[1]) * (v1 - v0);
        gpn = phi0 * g0[2] + phi1 * g1[2];
    /*        ----- South ----- */
        v0 = v[ll + v_dim1];
        v1 = v[ur + v_dim1];
        for (i11 = 0; i11 <= *d__; ++i11)
            g0[i11] = vval[i11 + c__[j * c_dim1 + 1] * vval_dim1];
        for (i11 = 0; i11 <= *d__; ++i11)
            g1[i11] = vval[i11 + c__[j * c_dim1 + 2] * vval_dim1];
        xibar = v[ll + (v_dim1 << 1)];
        m= xibar_search(a, t,  xi, xibar, nt);
        if (m >= 1) {
            m = lo[t[m - 1]];
            while  (a[m] != 0)
                if (z__[a[m]] <= xi[m])
                    m = lo[m];
                else
                    m = hi[m];
            if (v0 < v[c__[m * c_dim1 + 3] + v_dim1]) {
                v0 = v[c__[m * c_dim1 + 3] + v_dim1];
                for (i11 = 0; i11 <= *d__; ++i11)
                    g0[i11] = vval[i11 + c__[m * c_dim1 + 3] * vval_dim1];
            }
            if (v[c__[m * c_dim1 + 4] + v_dim1] < v1) {
                v1 = v[c__[m * c_dim1 + 4] + v_dim1];
                for (i11 = 0; i11 <= *d__; ++i11)
                    g1[i11] = vval[i11 + c__[m * c_dim1 + 4] * vval_dim1];
            }
        }
        h__ = (z__[1] - v0) / (v1 - v0);
    /*        Hermite basis */
        hermite_prep(h__, &phi0, &phi1, &psi0, &psi1);
        gs = phi0 * g0[0] + phi1 * g1[0] + (psi0 * g0[1] + psi1 * g1[1]) * (v1 - v0);
        gps = phi0 * g0[2] + phi1 * g1[2];
    /*        ----- East ----- */
        v0 = v[ll + (v_dim1 << 1)];
        v1 = v[ur + (v_dim1 << 1)];
        for (i11 = 0; i11 <= *d__; ++i11)
            g0[i11] = vval[i11 + c__[j * c_dim1 + 2] * vval_dim1];
        for (i11 = 0; i11 <= *d__; ++i11)
            g1[i11] = vval[i11 + c__[j * c_dim1 + 4] * vval_dim1];
        xibar = v[ur + v_dim1];
        m= xibar_search(a, t,  xi, xibar, nt);
        if (m >= 1) {
            m = hi[t[m - 1]];
            while (a[m] != 0)
                if (z__[a[m]] <= xi[m])
                    m = lo[m];
                else
                    m = hi[m];
            if (v0 < v[c__[m * c_dim1 + 1] + (v_dim1 << 1)]) {
                v0 = v[c__[m * c_dim1 + 1] + (v_dim1 << 1)];
                for (i11 = 0; i11 <= *d__; ++i11)
                    g0[i11] = vval[i11 + c__[m * c_dim1 + 1] * vval_dim1];
            }
            if (v[c__[m * c_dim1 + 3] + (v_dim1 << 1)] < v1) {
                v1 = v[c__[m * c_dim1 + 3] + (v_dim1 << 1)];
                for (i11 = 0; i11 <= *d__; ++i11)
                    g1[i11] = vval[i11 + c__[m * c_dim1 + 3] * vval_dim1];
            }
        }
        h__ = (z__[2] - v0) / (v1 - v0);
    /*        Hermite basis */
        hermite_prep(h__, &phi0, &phi1, &psi0, &psi1);
        ge = phi0 * g0[0] + phi1 * g1[0] + (psi0 * g0[2] + psi1 * g1[2]) * ( v1 - v0);
        gpe = phi0 * g0[1] + phi1 * g1[1];
    /*        ----- West ----- */
        v0 = v[ll + (v_dim1 << 1)];
        v1 = v[ur + (v_dim1 << 1)];
        for (i11 = 0; i11 <= *d__; ++i11)
            g0[i11] = vval[i11 + c__[j * c_dim1 + 1] * vval_dim1];
        for (i11 = 0; i11 <= *d__; ++i11)
            g1[i11] = vval[i11 + c__[j * c_dim1 + 3] * vval_dim1];
        xibar = v[ll + v_dim1];
        m = xibar_search(a, t,  xi, xibar, nt);
        if (m >= 1) {
            m = lo[t[m - 1]];
            while (a[m] != 0)
                if (z__[a[m]] <= xi[m])
                    m = lo[m];
                else
                    m = hi[m];
            if (v0 < v[c__[m * c_dim1 + 2] + (v_dim1 << 1)]) {
                v0 = v[c__[m * c_dim1 + 2] + (v_dim1 << 1)];
                for (i11 = 0; i11 <= *d__; ++i11)
                    g0[i11] = vval[i11 + c__[m * c_dim1 + 2] * vval_dim1];
            }
            if (v[c__[m * c_dim1 + 4] + (v_dim1 << 1)] < v1) {
                v1 = v[c__[m * c_dim1 + 4] + (v_dim1 << 1)];
                for (i11 = 0; i11 <= *d__; ++i11)
                    g1[i11] = vval[i11 + c__[m * c_dim1 + 4] * vval_dim1];
            }
        }
        h__ = (z__[2] - v0) / (v1 - v0);
    /*        Hermite basis */
        hermite_prep(h__, &phi0, &phi1, &psi0, &psi1);
        gw = phi0 * g0[0] + phi1 * g1[0] + (psi0 * g0[2] + psi1 * g1[2]) * (
            v1 - v0);
        gpw = phi0 * g0[1] + phi1 * g1[1];
    /*        NS */
        h__ = (z__[2] - v[ll + (v_dim1 << 1)]) / (v[ur + (v_dim1 << 1)] - v[
            ll + (v_dim1 << 1)]);
    /*        Hermite basis */
        hermite_prep(h__, &phi0, &phi1, &psi0, &psi1);
        sns = phi0 * gs + phi1 * gn + (psi0 * gps + psi1 * gpn) * (v[ur + (v_dim1 << 1)] - v[ll + (v_dim1 << 1)]);
    /*        EW */
        h__ = (z__[1] - v[ll + v_dim1]) / (v[ur + v_dim1] - v[ll + v_dim1]);
    /*        Hermite basis */
        hermite_prep(h__, &phi0, &phi1, &psi0, &psi1);
        sew = phi0 * gw + phi1 * ge + (psi0 * gpw + psi1 * gpe) * (v[ur + v_dim1] - v[ll + v_dim1]);
        s = sns + sew - s;
    }
    return s;
}

static void ehg136_(double *u, integer *lm, integer *m, integer *n, integer *d__, integer *nf, 
        double *f, double *x, integer *psi, double *y, double *rw, integer *kernel, integer *k, 
        double *dist, double *eta, double *b, integer *od, double *o, integer *ihat, double *w, 
        double *rcond, integer *sing, integer *dd, integer *tdeg, integer *cdeg, double * s) {

    static integer execnt = 0;
    integer o_dim1, b_dim1, b_offset, s_dim1, u_dim1, x_dim1, x_offset;
    static integer i__, j, l, i1, info, identi;
    static double q[8], tol, work[15], scale, sigma[15], qraux[15], dgamma[15];
    static double e[225]	/* was [15][15] */, g[225]	/* was [15][15] */;

    o_dim1 = *m;
    o -= 1 + o_dim1;
    --dist;
    --rw;
    --psi;
    x_dim1 = *n;
    x -= x_offset = 1 + x_dim1;
    u_dim1 = *lm;
    u -= 1 + u_dim1;
    --w;
    --eta;
    b_dim1 = *nf;
    b -= b_offset = 1 + b_dim1;
    s_dim1 = *od - 0 + 1;
    s -= 0 + s_dim1;
    --cdeg;

    ++execnt;
    if (! (*k <= *nf - 1))
        loess_error(104);
    if (! (*k <= 15))
        loess_error(105);
    for (identi = 1; identi <= *n; ++identi)
        psi[identi] = identi;
    for (l = 1; l <= *m; ++l) {
        for (i1 = 1; i1 <= *d__; ++i1)
            q[i1 - 1] = u[l + i1 * u_dim1];
        ehg127_(q, n, d__, nf, f, &x[x_offset], &psi[1], y, &rw[1],
            kernel, k, &dist[1], &eta[1], &b[b_offset], od, &w[1], rcond,
            sing, sigma, e, g, dgamma, qraux, work, &tol, dd, tdeg, &cdeg[1], &s[l * s_dim1]);
        if (*ihat == 1) {
    /*           $L sub {l,l} = */
    /*           V sub {1,:} SIGMA sup {+} U sup T */
    /*           (Q sup T W e sub i )$ */
            if (! (*m == *n))
                loess_error(123);
    /*           find $i$ such that $l = psi sub i$ */
            i__ = 1;
            while  (l != psi[i__]) {
                ++i__;
                if (! (i__ < *nf))
                    loess_error(123);
            }
            for (i1 = 1; i1 <= *nf; ++i1)
                eta[i1] = 0.;
            eta[i__] = w[i__];
    /*           $eta = Q sup T W e sub i$ */
            dqrsl_(&b[b_offset], nf, nf, k, qraux, &eta[1], &eta[1], &eta[1],
                &eta[1], &eta[1], &eta[1], 1000, &info);
    /*           $gamma = U sup T eta sub {1:k}$ */
            for (i1 = 1; i1 <= *k; ++i1)
                dgamma[i1 - 1] = 0.;
            for (j = 1; j <= *k; ++j)
                for (i1 = 1; i1 <= *k; ++i1)
                        dgamma[i1 - 1] += eta[j] * e[j + i1 * 15 - 16];
    /*           $gamma = SIGMA sup {+} gamma$ */
            for (j = 1; j <= *k; ++j)
                if (tol < sigma[j - 1])
                         dgamma[j - 1] /= sigma[j - 1];
                else
                        dgamma[j - 1] = 0.;
            o[l + o_dim1] = ddot_(k, g, &c__15, dgamma, &c__1);
        } else if (*ihat == 2) {
            /*     $L sub {l,:} = */
            /*     V sub {1,:} SIGMA sup {+} */
            /*     ( U sup T Q sup T ) W $ */
            for (i1 = 1; i1 <= *n; ++i1)
                    o[l + i1 * o_dim1] = 0.;
            for (j = 1; j <= *k; ++j) {
                for (i1 = 1; i1 <=  *nf; ++i1)
                    eta[i1] = 0.;
                for (i1 = 1; i1 <=  *k; ++i1)
                    eta[i1] = e[i1 + j * 15 - 16];
                dqrsl_(&b[b_offset], nf, nf, k, qraux, &eta[1], &eta[1],
                    work, work, work, work, 10000, &info);
                if (tol < sigma[j - 1])
                    scale = 1. / sigma[j - 1];
                else
                    scale = 0.;
                for (i1 = 1; i1 <=  *nf; ++i1)
                    eta[i1] *= scale * w[i1];
                for (i__ = 1; i__ <= *nf; ++i__)
                    o[l + psi[i__] * o_dim1] += g[j * 15 - 15] * eta[i__];
            }
        }
    }
} /* ehg136_ */

static void ehg137_(double *z__, integer *kappa, integer *leaf, integer *nleaf, integer *d__, 
        integer *nv, integer *nvmax, integer * ncmax, integer *a, double *xi, integer *lo, integer *hi) {
    static integer execnt = 0, p, pstack[20], stackt;

    --leaf;
    --z__;
    --hi;
    --lo;
    --xi;
    --a;
    /*     stacktop -> stackt */
    ++execnt;
    /*     find leaf cells affected by $z$ */
    stackt = 0;
    p = 1;
    *nleaf = 0;
    while (0 < p) {
        if (a[p] == 0) {
            // leaf
            ++(*nleaf);
            leaf[*nleaf] = p;
            // Pop
            if (stackt >= 1)
                p = pstack[stackt - 1];
            else
                p = 0;
            stackt = GSL_MAX(0, stackt-1);
        } else {
            if (z__[a[p]] == xi[p]) {
                // Push
                ++stackt;
                if (! (stackt <= 20))
                    loess_error(187);
                pstack[stackt - 1] = hi[p];
                p = lo[p];
            } else {
                p = (z__[a[p]] <= xi[p])
                        ? lo[p]
                        : hi[p];
            }
        }
    }
    if (! (*nleaf <= 256))
        loess_error(185);
} /* ehg137_ */


//BK: Changed phi (the 17th input) from integer to double, because this fn is only called
//once, and there, dist==phi.
static void ehg139_(double *v, integer *nvmax, integer *nv, integer *n, integer *d__, 
        integer *nf, double *f, double *x, integer *pi, integer *psi, double *y, 
        double *rw, double *trl, integer *kernel, integer *k, double *dist, double *phi,
        double *eta, double *b, integer *od, double *w, double *diagl, double *vval2, 
        integer *ncmax, integer *vc, integer *a, double *xi, integer *lo, integer *hi, integer *c__,
        integer *vhit, double *rcond, integer *sing, integer *dd, integer
        *tdeg, integer *cdeg, integer *lq, double *lf, logical *setlf, double *s) {

    integer lq_dim1, c_dim1, c_offset, lf_dim1, lf_dim2, b_dim1, b_offset, 
            s_dim1, v_dim1, v_offset, vval2_dim1, vval2_offset, x_dim1, x_offset, i__1, i__3;

    static integer execnt = 0;
    static double e[225]	/* was [15][15] */;
    static double q[8], u[225]	/* was [15][15] */, z__[8], i4, i7, tol;
    static integer i__, j, l, i5, i6, ii, leaf[256], info, ileaf, nleaf, identi;
    static double term, work[15], scale, sigma[15], qraux[15], dgamma[15];

    --vhit; --diagl; --phi; --dist; --rw; --y;
    --psi; --pi; --w; --eta; --hi; --lo; --xi; --cdeg;
    vval2_dim1 = *d__ - 0 + 1;
    vval2 -= vval2_offset = 0 + vval2_dim1;
    x_dim1 = *n;
    x -= x_offset = 1 + x_dim1;
    v_dim1 = *nvmax;
    v -= v_offset = 1 + v_dim1;
    lf_dim1 = *d__ - 0 + 1;
    lf_dim2 = *nvmax;
    lf -=  0 + lf_dim1 * (1 + lf_dim2);
    lq_dim1 = *nvmax;
    lq -= 1 + lq_dim1;
    b_dim1 = *nf;
    b -= b_offset = 1 + b_dim1;
    s_dim1 = *od - 0 + 1;
    s -= 0 + s_dim1;
    c_dim1 = *vc;
    c__ -= c_offset = 1 + c_dim1;

    ++execnt;
    /*     l2fit with trace(L) */
    if (! (*k <= *nf - 1))
        loess_error(104);
    if (! (*k <= 15))
        loess_error(105);
    if (*trl != 0.) {
        for (i5 = 1; i5 <= *n; ++i5)
            diagl[i5] = 0.;
        for (i6 = 1; i6 <= *nv; ++i6)
            for (i5 = 0; i5 <= *d__; ++i5)
                vval2[i5 + i6 * vval2_dim1] = 0.;
    }
    for (identi = 1; identi <= *n; ++identi)
        psi[identi] = identi;
    i__1 = *nv;
    for (l = 1; l <= i__1; ++l) {
        for (i5 = 1; i5 <= *d__; ++i5)
            q[i5 - 1] = v[l + i5 * v_dim1];
        ehg127_(q, n, d__, nf, f, &x[x_offset], &psi[1], &y[1], &rw[1],
            kernel, k, dist, &eta[1], &b[b_offset], od, &w[1], rcond,
            sing, sigma, u, e, dgamma, qraux, work, &tol, dd, tdeg, &cdeg[1], &s[l * s_dim1]);
        if (*trl != 0.) { // invert $psi$
            for (i5 = 1; i5 <= *n; ++i5)
                phi[i5] = 0;
            for (i__ = 1; i__ <= *nf; ++i__)
                phi[psi[i__]] = i__;
            for (i5 = 1; i5 <= *d__; ++i5)
                z__[i5 - 1] = v[l + i5 * v_dim1];
            ehg137_(z__, &vhit[l], leaf, &nleaf, d__, nv, nvmax, ncmax, a, &xi[1], &lo[1], &hi[1]);
            for (ileaf = 1; ileaf <= nleaf; ++ileaf) {
                i__3 = hi[leaf[ileaf - 1]];
                for (ii = lo[leaf[ileaf - 1]]; ii <= i__3; ++ii) {
                    i__ = phi[pi[ii]];
                    if (i__ != 0) {
                        if (! (psi[i__] == pi[ii]))
                            loess_error(194);
                        for (i5 = 1; i5 <= *nf; ++i5)
                            eta[i5] = 0.;
                        eta[i__] = w[i__];
                        /*                    $eta = Q sup T W e sub i$ */
                        dqrsl_(&b[b_offset], nf, nf, k, qraux, &eta[1], work,
                            &eta[1], &eta[1], work, work, 1000, &info);
                        for (j = 1; j <= *k; ++j) {
                            i4 = (tol < sigma[j - 1])
                                ? ddot_(k, &u[j * 15 - 15], &c__1, &eta[1], &c__1) / sigma[j - 1]
                                : 0.;
                            dgamma[j - 1] = i4;
                        }
                        for (j = 1; j <= *d__ + 1; ++j) // bug fix 2006-07-15 for k=1, od>1.   (thanks btyner@gmail.com) */
                            vval2[j - 1 + l * vval2_dim1] = (j <= *k)
                                        ? ddot_(k, &e[j - 1], &c__15, dgamma, &c__1)
                                        : 0.;
                        for (i5 = 1; i5 <= *d__; ++i5)
                            z__[i5 - 1] = x[pi[ii] + i5 * x_dim1];
                        term = ehg128_(z__, d__, ncmax, vc, a, &xi[1], &
                            lo[1], &hi[1], &c__[c_offset], &v[v_offset],
                            nvmax, &vval2[vval2_offset]);
                        diagl[pi[ii]] += term;
                        for (i5 = 0; i5 <= *d__; ++i5)
                            vval2[i5 + l * vval2_dim1] = 0.;
                    }
                }
            }
        }
        if (*setlf) {
            /*           $Lf sub {:,l,:} = V SIGMA sup {+} U sup T Q sup T W$ */
            if (! (*k >= *d__ + 1))
                loess_error(196);
            for (i5 = 1; i5 <= *nf; ++i5)
                lq[l + i5 * lq_dim1] = psi[i5];
            for (i6 = 1; i6 <= *nf; ++i6)
                for (i5 = 0; i5 <= *d__; ++i5)
                    lf[i5 + (l + i6 * lf_dim2) * lf_dim1] = 0.;
            for (j = 1; j <= *k; ++j) {
                for (i5 = 1; i5 <= *nf; ++i5)
                    eta[i5] = 0.;
                for (i5 = 1; i5 <= *k; ++i5)
                    eta[i5] = u[i5 + j * 15 - 16];
                dqrsl_(&b[b_offset], nf, nf, k, qraux, &eta[1], &eta[1], work, work, work, work, 10000, &info);
                scale = (tol < sigma[j - 1])
                        ? 1. / sigma[j - 1]
                        : 0.;
                for (i5 = 1; i5 <= *nf; ++i5)
                    eta[i5] *= scale * w[i5];
                for (i__ = 1; i__ <= *nf; ++i__) {
                    i7 = eta[i__];
                    for (i5 = 0; i5 <= *d__; ++i5)
                        if (i5 < *k)
                            lf[i5 + (l + i__ * lf_dim2) * lf_dim1] += e[i5 + 1 + j * 15 - 16] * i7;
                        else
                            lf[i5 + (l + i__ * lf_dim2) * lf_dim1] = 0.;
                }
            }
        }
    }
    if (*trl != 0.) {
        if (*n <= 0)
            *trl = 0.;
        else {
            *trl = diagl[*n];
            for (int i2 = *n - 1; i2 >= 1; --i2)
                *trl += diagl[i2];
        }
    }
}

static void dqrdc_(double *x, integer *ldx, integer *n, integer *p, double *qraux, 
        integer *jpvt, double *work, integer job) {
    integer x_dim1, i__2, i__3;
    double d__2;

    static integer j, l, jj, jp, pl, pu, lp1, lup, maxj;
    static logical negj, swapj;
    static double t, tt, nrmxl, maxnrm;

/*     dqrdc uses householder transformations to compute the qr 
     factorization of an n by p matrix x.  column pivoting 
     based on the 2-norms of the reduced columns may be 
     performed at the users option. 

     on entry 

        x       double precision(ldx,p), where ldx .ge. n. 
                x contains the matrix whose decomposition is to be 
                computed. 

        ldx     integer. 
                ldx is the leading dimension of the array x. 

        n       integer. 
                n is the number of rows of the matrix x. 

        p       integer. 
                p is the number of columns of the matrix x. 

        jpvt    integer(p). 
                jpvt contains integers that control the selection 
                of the pivot columns.  the k-th column x(k) of x 
                is placed in one of three classes according to the 
                value of jpvt(k). 

                   if jpvt(k) .gt. 0, then x(k) is an initial 
                                      column. 

                   if jpvt(k) .eq. 0, then x(k) is a free column. 

                   if jpvt(k) .lt. 0, then x(k) is a final column. 

                before the decomposition is computed, initial columns 
                are moved to the beginning of the array x and final 
                columns to the end.  both initial and final columns 
                are frozen in place during the computation and only 
                free columns are moved.  at the k-th stage of the 
                reduction, if x(k) is occupied by a free column 
                it is interchanged with the free column of largest 
                reduced norm.  jpvt is not referenced if 
                job .eq. 0. 

        work    double precision(p). 
                work is a work array.  work is not referenced if 
                job .eq. 0. 

        job     integer. 
                job is an integer that initiates column pivoting. 
                if job .eq. 0, no pivoting is done. 
                if job .ne. 0, pivoting is done. 

     on return 

        x       x contains in its upper triangle the upper 
                triangular matrix r of the qr factorization. 
                below its diagonal x contains information from 
                which the orthogonal part of the decomposition 
                can be recovered.  note that if pivoting has 
                been requested, the decomposition is not that 
                of the original matrix x but that of x 
                with its columns permuted as described by jpvt. 

        qraux   double precision(p). 
                qraux contains further information required to recover 
                the orthogonal part of the decomposition. 

        jpvt    jpvt(k) contains the index of the column of the 
                original matrix that has been interchanged into 
                the k-th column, if pivoting was requested. 

     linpack. this version dated 08/14/78 . 
     g.w. stewart, university of maryland, argonne national lab. 

     dqrdc uses the following functions and subprograms. 

     blas daxpy,ddot,dscal,dswap,dnrm2 
     fortran dabs,dmax1,min0,dsqrt */

/*     internal variables */

    x_dim1 = *ldx;
    x -= 1 + x_dim1;
    --qraux;
    --jpvt;
    --work;

    pl = 1;
    pu = 0;
    if (job == 0)
        goto L60;

    /*        pivoting has been requested.  rearrange the columns */
    /*        according to jpvt. */
    for (j = 1; j <= *p; ++j) {
        swapj = jpvt[j] > 0;
        negj = jpvt[j] < 0;
        jpvt[j] = j;
        if (negj)
            jpvt[j] = -j;
        if (! swapj)
            goto L10;
        if (j != pl)
            dswap(*n, &x[pl * x_dim1 + 1], &x[j * x_dim1 + 1]);
        jpvt[j] = jpvt[pl];
        jpvt[pl] = j;
        ++pl;
    L10:
        ;
    }
    pu = *p;
    for (jj = 1; jj <= *p; ++jj) {
        j = *p - jj + 1;
        if (jpvt[j] >= 0)
            continue;
        jpvt[j] = -jpvt[j];
        if (j == pu)
            goto L30;
        dswap(*n, &x[pu * x_dim1 + 1], &x[j * x_dim1 + 1]);
        jp = jpvt[pu];
        jpvt[pu] = jpvt[j];
        jpvt[j] = jp;
    L30:
        --pu;
    }
L60:

    /*     compute the norms of the free columns. */
    if (!(pu < pl) )
        for (j = pl; j <= pu; ++j) {
            qraux[j] = dnrm2(n, &x[j * x_dim1 + 1]);
            work[j] = qraux[j];
        }

    /*     perform the householder reduction of x. */
    lup = min(*n,*p);
    for (l = 1; l <= lup; ++l) {
        if (l < pl || l >= pu)
            goto L120;

        /*           locate the column of largest norm and bring it */
        /*           into the pivot position. */
        maxnrm = 0.;
        maxj = l;
        for (j = l; j <= pu; ++j) {
            if (qraux[j] <= maxnrm)
                continue;
            maxnrm = qraux[j];
            maxj = j;
        }
        if (maxj != l) {
            dswap(*n, &x[l * x_dim1 + 1], &x[maxj * x_dim1 + 1]);
            qraux[maxj] = qraux[l];
            work[maxj] = work[l];
            jp = jpvt[maxj];
            jpvt[maxj] = jpvt[l];
            jpvt[l] = jp;
        }
    L120:
        qraux[l] = 0.;
        if (l == *n)
            continue;
        /*           compute the householder transformation for column l. */

        i__2 = *n - l + 1;
        nrmxl = dnrm2(&i__2, &x[l + l * x_dim1]);
        if (nrmxl == 0.)
            continue;
        if (x[l + l * x_dim1] != 0.)
            nrmxl = d_sign(&nrmxl, &x[l + l * x_dim1]);
        i__2 = *n - l + 1;
        dscal(&i__2, 1./nrmxl, &x[l + l * x_dim1]);
        x[l + l * x_dim1] += 1.;

        /*  apply the transformation to the remaining columns,  updating the norms. */
        lp1 = l + 1;
        if (!(*p < lp1))
            for (j = lp1; j <= *p; ++j) {
                i__3 = *n - l + 1;
                t = -ddot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
                    c__1) / x[l + l * x_dim1];
                i__3 = *n - l + 1;
                daxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
                    c__1);
                if ((j < pl || j > pu) || (qraux[j] == 0.) )
                    continue;
                /* Computing 2nd power */
                d__2 = abs(x[l + j * x_dim1]) / qraux[j];
                tt = GSL_MAX(1. - d__2 * d__2, 0);
                t = tt;
                tt = tt * .05 * gsl_pow_2(qraux[j] / work[j]) + 1.;
                if (tt == 1.)
                    goto L130;
                qraux[j] *= sqrt(t);
                continue;
                L130:
                i__3 = *n - l;
                qraux[j] = dnrm2(&i__3, &x[l + 1 + j * x_dim1]);
                work[j] = qraux[j];
            }
        /*              save the transformation. */
        qraux[l] = x[l + l * x_dim1];
        x[l + l * x_dim1] = -nrmxl;
    }
} /* dqrdc_ */

static void lowesb_(double *xx, double *yy, double *ww, double *diagl, double trl,
        integer *iv, integer *liv, integer * lv, double *wv) {
    static integer execnt = 0, setlf;
    --wv;
    --iv;

    ++execnt;
    if (! (iv[28] != 173))
        loess_error(174);
    if (iv[28] != 172 && !(iv[28] == 171))
	    loess_error(171);
    iv[28] = 173;
    setlf = iv[27] != iv[25];
    integer     i__1 = ifloor(iv[3] * wv[2]);
    ehg131_(xx, yy, ww, &trl, diagl, &iv[20], &iv[29], &iv[3], &iv[2], &iv[5], &iv[17], 
            &iv[4], &iv[6], &iv[14], &iv[19], &wv[1] , &iv[iv[7]], &iv[iv[8]], &iv[iv[9]], &iv[iv[10]], 
            &iv[iv[22]], & iv[iv[27]], &wv[iv[11]], &iv[iv[23]], &wv[iv[13]], &wv[iv[12]], & wv[iv[15]], 
            &wv[iv[16]], &wv[iv[18]], &i__1, &wv[3], &wv[iv[26]], &wv[iv[24]], &wv[4], &iv[30], &iv[33], 
            &iv[32], &iv[41], &iv[iv[25]], &wv[iv[34]], &setlf);
    if ((double) iv[14] < iv[6] + (double) iv[4] / 2.)
        ehg183_("Warning. k-d tree limited by memory; nvmax=", &iv[14], 1, 1);
    else if (iv[17] < iv[5] + 2)
	    ehg183_("Warning. k-d tree limited by memory. ncmax=", &iv[17], 1, 1);
} /* lowesb_ */

static void lowesd_(integer *iv, integer *liv, integer *lv, double *v, 
        integer d__, integer n, double f, integer ideg, integer *nvmax, logical *setlf) {
    static integer execnt = 0, i__, j, i1, i2, nf, vc, ncmax, bound;
    --iv;
    --v;

    ++execnt;
    iv[28] = 171;
    iv[2] = d__;
    iv[3] = n;
    vc = pow_ii(2, d__);
    iv[4] = vc;
    if (! (0. < f))
        loess_error(120);
/* Computing MIN */
    nf = GSL_MIN(n, ifloor(n * f));
    iv[19] = nf;
    iv[20] = 1;
    if (ideg == 0)
        i1 = 1;
    else {
        if (ideg == 1)
            i1 = d__ + 1;
        else if (ideg == 2)
            i1 = (integer) ((double) ((d__ + 2) * (d__ + 1)) / 2.);
    }
    iv[29] = i1;
    iv[21] = 1;
    iv[14] = *nvmax;
    ncmax = *nvmax;
    iv[17] = ncmax;
    iv[30] = 0;
    iv[32] = ideg;
    if (! (ideg >= 0) || (! (ideg <= 2)))
        loess_error(195);
    iv[33] = d__;
    for (i2 = 41; i2 <= 49; ++i2)
        iv[i2] = ideg;
    iv[7] = 50;
    iv[8] = iv[7] + ncmax;
    iv[9] = iv[8] + vc * ncmax;
    iv[10] = iv[9] + ncmax;
    iv[22] = iv[10] + ncmax;
/*     initialize permutation */
    j = iv[22] - 1;
    for (i__ = 1; i__ <= n; ++i__)
        iv[j + i__] = i__;
    iv[23] = iv[22] + n;
    iv[25] = iv[23] + *nvmax;
    if (*setlf)
        iv[27] = iv[25] + *nvmax * nf;
    else
        iv[27] = iv[25];
    bound = iv[27] + n;
    if (! (bound - 1 <= *liv))
        loess_error(102);
    iv[11] = 50;
    iv[13] = iv[11] + *nvmax * d__;
    iv[12] = iv[13] + (d__ + 1) * *nvmax;
    iv[15] = iv[12] + ncmax;
    iv[16] = iv[15] + n;
    iv[18] = iv[16] + nf;
    iv[24] = iv[18] + iv[29] * nf;
    iv[34] = iv[24] + (d__ + 1) * *nvmax;
    if (*setlf)
        iv[26] = iv[34] + (d__ + 1) * *nvmax * nf;
    else
        iv[26] = iv[34];
    bound = iv[26] + nf;
    if (! (bound - 1 <= *lv))
       loess_error(103);

    v[1] = f;
    v[2] = .05;
    v[3] = 0.;
    v[4] = 1.;
} /* lowesd_ */

static void lowese_(integer *iv, integer *liv, integer *lv, double *wv, integer m, double *z, double *s) {
    static integer execnt = 0;
    ++execnt;
    --iv;
    --wv;

    if (! (iv[28] != 172))
        loess_error(172);
    if (! (iv[28] == 173))
        loess_error(173);
    ehg133_(&iv[3], &iv[2], &iv[4], &iv[14], &iv[5], &iv[17], &iv[iv[7]], &iv[iv[8]], &iv[iv[9]],
            &iv[iv[10]], &wv[iv[11]], &wv[iv[13]], &wv[iv[12]], m, z, s);
}

static void lowesf_(double *xx, double *yy, double *ww, integer *iv, integer *liv, 
        integer *lv, double *wv, integer *m, double *z__, double *l, integer ihat, double *s) {
    static integer execnt = 0;
    integer l_dim1, l_offset, z_dim1, z_offset;
    static logical i1;
    --xx;
    --yy;
    --ww;
    --iv;
    --wv;
    l_dim1 = *m;
    l -= l_offset = 1 + l_dim1;
    z_dim1 = *m;
    z__ -= z_offset = 1 + z_dim1;

    ++execnt;
    i1 = (171 <= iv[28])
          ? iv[28] <= 174
	      : FALSE_;
    if (! i1)
        loess_error(171);
    iv[28] = 172;
    if (! (iv[14] >= iv[19]))
        loess_error(186);
    ehg136_(&z__[z_offset], m, m, &iv[3], &iv[2], &iv[19], &wv[1], &xx[1], &iv[iv[22]], &yy[1],
            &ww[1], &iv[20], &iv[29], &wv[iv[15]], &wv[iv[16]], &wv[iv[18]], &c__0, &l[l_offset],
            &ihat, &wv[iv[26]], &wv[4], &iv[30], &iv[33], &iv[32], &iv[41], s);
} /* lowesf_ */

static void lowesl_(integer *iv, integer *liv, integer *lv, double *wv, integer *m, double *z__, double *l) {
    static integer execnt = 0;
    integer l_dim1, l_offset, z_dim1, z_offset;

    --iv;
    --wv;
    l_dim1 = *m;
    l -= l_offset = 1 + l_dim1;
    z_dim1 = *m;
    z__ -= z_offset = 1 + z_dim1;

    ++execnt;
    if (! (iv[28] != 172))
        loess_error(172);
    if (! (iv[28] == 173))
        loess_error(173);
    if (! (iv[26] != iv[34]))
        loess_error(175);
    ehg191_(m, &z__[z_offset], &l[l_offset], &iv[2], &iv[3], &iv[19], &iv[6],
	    &iv[17], &iv[4], &iv[iv[7]], &wv[iv[12]], &iv[iv[10]], &iv[iv[9]],
	     &iv[iv[8]], &wv[iv[11]], &iv[14], &wv[iv[24]], &wv[iv[34]], &iv[iv[25]]);
} /* lowesl_ */

static void lowesw_(double *res, integer *n, double *rw, integer *pi) {
    static integer i1, nh, execnt = 0, identi;
    static double cmad, rsmall;
    --pi;
    --rw;
    --res;
    ++execnt;
/*     tranliterated from Devlin's ratfor */
/*     find median of absolute residuals */
    for (i1 = 1; i1 <= *n; ++i1)
        rw[i1] = abs(res[i1]);
    for (identi = 1; identi <= *n; ++identi)
        pi[identi] = identi;
    nh = ifloor((double) (*n) / 2.) + 1;
/*     partial sort to find 6*mad */
    find_kth_smallest(1, *n, nh, 1, &rw[1], &pi[1]);
    if (*n - nh + 1 < nh) {
        find_kth_smallest(1, 1, 2, 1, &rw[1], &pi[1]);
        cmad = (rw[pi[nh]] + rw[pi[nh - 1]]) * 3;
    } else
        cmad = rw[pi[nh]] * 6;
    rsmall = DBL_MIN;
    if (cmad < rsmall)
        for (i1 = 1; i1 <= *n; ++i1)
            rw[i1] = 1.;
    else
        for (int i__ = 1; i__ <= *n; ++i__)
            if (cmad * .999 < rw[i__])
                rw[i__] = 0.;
            else
                rw[i__] = (cmad * .001 < rw[i__])
                            ? gsl_pow_2(1 - gsl_pow_2(rw[i__] / cmad))
                            : 1.;
}

static void pseudovals(integer n, double *y, double *yhat, double *pwgts,  //formerly lowesp
                double *rwgts, integer *pi, double *ytilde) {
    static integer m, i5, identi, execnt = 0;
    static double i4, mad;

    --ytilde;
    --pi;
    --rwgts;
    --pwgts;
    --yhat;
    --y;

    ++execnt;
    /*     median absolute deviation */
    for (i5 = 1; i5 <= n; ++i5)
        ytilde[i5] = abs(y[i5] - yhat[i5]) * sqrt(pwgts[i5]);
    for (identi = 1; identi <= n; ++identi)
        pi[identi] = identi;
    m = ifloor((double) (n) / 2.) + 1;
    find_kth_smallest(1, n, m, 1, &ytilde[1], &pi[1]);
    if (n - m + 1 < m) {
        find_kth_smallest(1, m-1, m-1, 1, &ytilde[1], &pi[1]);
        mad = (ytilde[pi[m - 1]] + ytilde[pi[m]]) / 2;
    } else
        mad = ytilde[pi[m]];
    for (i5 = 1; i5 <= n; ++i5)
        ytilde[i5] = 1 - gsl_pow_2(y[i5] - yhat[i5]) * pwgts[i5] / (gsl_pow_2(mad * 6)/ 5);
    for (i5 = 1; i5 <= n; ++i5)
        ytilde[i5] *= sqrt(rwgts[i5]);
    if (n <= 0)
        i4 = 0.;
    else {
        double i1 = ytilde[n];
        for (integer i2 = n - 1; i2 >= 1; --i2)
            i1 += ytilde[i2];
        i4 = i1;
    }
    //     pseudovalues
    for (i5 = 1; i5 <= n; ++i5)
        ytilde[i5] = yhat[i5] + (n / i4) * rwgts[i5] * (y[i5] - yhat[i5]);
}

////// Back to loessc.c
static void loess_workspace(long D, long N, double	span, long degree,
			long *nonparametric, long *drop_square, long *sum_drop_sqr, long setLf){
	long tau0, nvmax, nf, i;
	nvmax = max(200, N);
        nf = min(N, floor(N * span));
        tau0 = (degree > 1) ? ((D + 2) * (D + 1) * 0.5) : (D + 1);
        tau = tau0 - (*sum_drop_sqr);
        lv = 50 + (3 * D + 3) * nvmax + N + (tau0 + 2) * nf;
	liv = 50 + ((long)pow((double)2, (double)D) + 4) * nvmax + 2 * N;
	if(setLf) {
		lv = lv + (D + 1) * nf * nvmax;
		liv = liv + nf * nvmax;	
	}
    iv = Calloc(liv, long);
    v = Calloc(lv, double);

    lowesd_(iv, &liv, &lv, v, D, N, span, degree, &nvmax, &setLf);
    iv[32] = *nonparametric;
    for(i = 0; i < D; i++)
        iv[i + 40] = drop_square[i];
}

static void loess_free() {
    free(v);
    free(iv);
}

static void loess_dfit( double	*y, double *x, double *x_evaluate, double *weights,
			double span, long degree, long *nonparametric, long *drop_square,
			long *sum_drop_sqr, long d, long n, long *m, double *fit) {
    loess_workspace(d, n, span, degree, nonparametric, drop_square, sum_drop_sqr, 0);
	lowesf_(x, y, weights, iv, &liv, &lv, v, m, x_evaluate, &doublepluszero, 0, fit);
	loess_free();
}

static void loess_dfitse( double	*y, double *x, double *x_evaluate, double *weights, double *robust,
        int	family, double span, long degree, long *nonparametric, long *drop_square,
         long *sum_drop_sqr, long d, long n, long *m, double *fit, double *L) {
    loess_workspace(d, n, span, degree, nonparametric, drop_square, sum_drop_sqr, 0);
	if(family == GAUSSIAN)
		lowesf_(x, y, weights, iv, &liv, &lv, v, m, x_evaluate, L, 2, fit);
	else if(family == SYMMETRIC) {
		lowesf_(x, y, weights, iv, &liv, &lv, v, m, x_evaluate, L, 2, fit);
		lowesf_(x, y, robust, iv, &liv, &lv, v, m, x_evaluate, &doublepluszero, 0, fit);
	}	
	loess_free();
}

static void loess_grow(long	const * restrict parameter,long const*restrict a,
                       double	const *restrict xi, double const *restrict vert, 
                       const double *restrict vval) {
	long	d, vc, nc, nv, a1, v1, xi1, vv1, i, k;
	d = parameter[0];
	vc = parameter[2];
	nc = parameter[3];
	nv = parameter[4];
	liv = parameter[5];
	lv = parameter[6];
	iv = Calloc(liv, long);
	v = Calloc(lv, double);

	iv[1] = d;
	iv[2] = parameter[1];
	iv[3] = vc;
	iv[5] = iv[13] = nv;
	iv[4] = iv[16] = nc;
	iv[6] = 50;
	iv[7] = iv[6] + nc;
	iv[8] = iv[7] + vc * nc;
	iv[9] = iv[8] + nc;
	iv[10] = 50;
	iv[12] = iv[10] + nv * d;
	iv[11] = iv[12] + (d + 1) * nv;
	iv[27] = 173;

	v1 = iv[10] - 1;
	xi1 = iv[11] - 1;
	a1 = iv[6] - 1;
	vv1 = iv[12] - 1;
	
    for(i = 0; i < d; i++) {
		k = nv * i;
		v[v1 + k] = vert[i];
		v[v1 + vc - 1 + k] = vert[i + d];
	}
    for(i = 0; i < nc; i++) {
            v[xi1 + i] = xi[i];
            iv[a1 + i] = a[i];
    }
	k = (d + 1) * nv;
	for(i = 0; i < k; i++)
		v[vv1 + i] = vval[i];
	ehg169_(d, &vc, &nc, &nc, &nv, nv, v+v1, iv+a1, v+xi1, iv+iv[7]-1, iv+iv[8]-1, iv+iv[9]-1);
}

static void loess_ifit(long const * restrict parameter, long const *restrict a, 
                double const *restrict xi, double const *restrict vert,
                 const double *restrict vval, long m, double *x_evaluate, double *fit) {
	loess_grow(parameter, a, xi, vert, vval);
	lowese_(iv, &liv, &lv, v, m, x_evaluate, fit);
	loess_free();
}

static void loess_ise( double	*y, double *x, double *x_evaluate, double *weights, double span, long degree,
             long int *nonparametric, long int *drop_square, long int *sum_drop_sqr, double *cell, long int d,
             long int n, long int *m, double *fit, double *L) {
    loess_workspace(d, n, span, degree, nonparametric, drop_square, sum_drop_sqr, 1);
	v[1] = *cell;
	lowesb_(x, y, weights, &doublepluszero, 0, iv, &liv, &lv, v);
	lowesl_(iv, &liv, &lv, v, m, x_evaluate, L);
	loess_free();
}

static void loess_prune( long	*parameter, long *a, double	*xi, double *vert, double *vval) {
	long	d, vc, a1, v1, xi1, vv1, nc, nv, nvmax, i, k;
	d = iv[1];
	vc = iv[3] - 1;
	nc = iv[4];
	nv = iv[5];
	a1 = iv[6] - 1;
	v1 = iv[10] - 1;
	xi1 = iv[11] - 1;
	vv1 = iv[12] - 1;
	nvmax = iv[13];

	for(i = 0; i < 5; i++)
		parameter[i] = iv[i + 1];
	parameter[5] = iv[21] - 1;
	parameter[6] = iv[14] - 1;

	for(i = 0; i < d; i++){
		k = nvmax * i;
		vert[i] = v[v1 + k];
		vert[i + d] = v[v1 + vc + k];
	}
	for(i = 0; i < nc; i++) {
		xi[i] = v[xi1 + i];
		a[i] = iv[a1 + i];
	}
	k = (d + 1) * nv;
	for(i = 0; i < k; i++)
		vval[i] = v[vv1 + i];
}

////// predict.c

struct pred_struct {
	double	*fit;           //The evaluated loess surface at eval.
	double	*se_fit;        //Estimates of the standard errors of the surface values.
	double  residual_scale; //Estimate of the scale of the residuals.
	double  df;             //The degrees of freedom of the t-distribution used to compute pointwise 
                            //   confidence intervals for the evaluated surface. 
};

void predict(double  *new_x, long M, struct loess_struct *lo, struct pred_struct *pre, int want_cov) {
	
    long D = lo->in.p;//Aliases for the purposes of merging some fn.s
    long N = lo->in.n;
            
	int     i, j, k, p;
	double x[N * D], x_tmp[N * D], x_evaluate[M * D];

	for(i = 0; i < D; i++) {
		k = i * M;
		for(j = 0; j < M; j++) {
			p = k + j;
			new_x[p] /= lo->out.divisor[i];
		}
	}
    memcpy(x_tmp, lo->in.x, N * D * sizeof(double));
	if(!strcmp(lo->control.surface, "direct") || want_cov)
        for(i = 0; i < D; i++) {
			k = i * N;
			for(j = 0; j < N; j++) {
                p = k + j;
                x_tmp[p] = lo->in.x[p] / lo->out.divisor[i];
            }
		}
	j = D - 1;
    long sum_drop_sqr = 0, sum_parametric = 0, nonparametric = 0;
    long order_parametric[D], order_drop_sqr[D];
	for(i = 0; i < D; i++) {
        sum_drop_sqr += lo->model.drop_square[i];
        sum_parametric += lo->model.parametric[i];
        if(lo->model.parametric[i])
            order_parametric[j--] = i;
        else
            order_parametric[nonparametric++] = i;
	}
    for(i = 0; i < D; i++) {
        order_drop_sqr[i] = 2 - lo->model.drop_square[order_parametric[i]];
        k = i * M;
        p = order_parametric[i] * M;
        for(j = 0; j < M; j++)
            x_evaluate[k + j] = new_x[p + j];
        k = i * N;
        p = order_parametric[i] * N;
        for(j = 0; j < N; j++)
            x[k + j] = x_tmp[p + j];
    }
	for(i = 0; i < N; i++)
		lo->out.robust[i] *= lo->in.weights[i];

    pre->fit = malloc(M * sizeof(double));
	pre->residual_scale = lo->out.s;
	pre->df = (lo->out.one_delta * lo->out.one_delta) / lo->out.two_delta;
    double L[N * M];
	if(!strcmp(lo->control.surface, "direct")) {
        if(want_cov)
            loess_dfitse(lo->in.y, x, x_evaluate, lo->in.weights, lo->out.robust, !strcmp(lo->model.family, "gaussian"), 
                lo->model.span, lo->model.degree, &nonparametric, order_drop_sqr, &sum_drop_sqr, D, N, &M, pre->fit, L);
        else
            loess_dfit(lo->in.y, x, x_evaluate, lo->out.robust, lo->model.span, lo->model.degree, &nonparametric,
                order_drop_sqr, &sum_drop_sqr, D, N, &M, pre->fit);
    } else {
        loess_ifit(lo->kd_tree.parameter, lo->kd_tree.a, lo->kd_tree.xi, lo->kd_tree.vert, 
                        lo->kd_tree.vval, M, x_evaluate, pre->fit);
        if(want_cov) {
            double new_cell = lo->model.span * lo->control.cell;
            double fit_tmp[M];
            loess_ise(lo->in.y, x, x_evaluate, lo->in.weights, lo->model.span, lo->model.degree, &nonparametric, 
                    order_drop_sqr, &sum_drop_sqr, &new_cell, D, N, &M, fit_tmp, L);
        }
    }
	if (want_cov) {
        pre->se_fit = malloc(M * sizeof(double));
        for(i = 0; i < N; i++) {
            k = i * M;
            for(j = 0; j < M; j++) {
                p = k + j;
                L[p] /= lo->in.weights[i];
                L[p] *= L[p]; //i.e., square
            }
		}
		for(i = 0; i < M; i++) {
            double tmp = 0;
			for(j = 0; j < N; j++)
                tmp += L[i + j * M];
			pre->se_fit[i] = lo->out.s * sqrt(tmp);
		}
	}
}

void pred_free_mem(struct	pred_struct	*pre){
	free(pre->fit);
	free(pre->se_fit);
}

 ///// loess.c
static  char    *surf_stat;

int comp(const void *d1_in, const void *d2_in) {
    const double *d1 = d1_in;
    const double *d2 = d1_in;
        if(*d1 < *d2)
                return(-1);
        else if(*d1 == *d2)
                return(0);
        else
                return(1);
}

static void condition(char	**surface, char *new_stat, char **trace_hat_in) {
	if(!strcmp(*surface, "interpolate")) {
		if(!strcmp(new_stat, "none"))
			surf_stat = "interpolate/none";
		else if(!strcmp(new_stat, "exact"))
			surf_stat = "interpolate/exact";
		else if(!strcmp(new_stat, "approximate"))
		{
			if(!strcmp(*trace_hat_in, "approximate"))
				surf_stat = "interpolate/2.approx";
			else if(!strcmp(*trace_hat_in, "exact"))
				surf_stat = "interpolate/1.approx";
		}
	}
	else if(!strcmp(*surface, "direct")) {
		if(!strcmp(new_stat, "none"))
			surf_stat = "direct/none";
		else if(!strcmp(new_stat, "exact"))
			surf_stat = "direct/exact";
		else if(!strcmp(new_stat, "approximate"))
			surf_stat = "direct/approximate";
	}
}

static void loess_raw( double	*y, double *x, double *weights, double *robust, long	*d, 
            long *n, double *span, long *degree, long *nonparametric, long *drop_square, 
            long *sum_drop_sqr, double *cell, char	**surf_stat, double *surface, long	*parameter, 
            long *a, double *xi, double *vert, double	*vval,double *diagonal, double*trL, 
            double*one_delta, double*two_delta, long *setLf) {
	long nsing, i, k;
	double	*hat_matrix, *LL;
	*trL = 0;
	loess_workspace(*d, *n, *span, *degree, nonparametric, drop_square, sum_drop_sqr, *setLf);
        v[1] = *cell;
	if(!strcmp(*surf_stat, "interpolate/none")) {
		lowesb_(x, y, robust, &doublepluszero, 0, iv, &liv, &lv, v);
		lowese_(iv, &liv, &lv, v, *n, x, surface);
		loess_prune(parameter, a, xi, vert, vval);
	}			
	else if (!strcmp(*surf_stat, "direct/none"))
		lowesf_(x, y, robust, iv, &liv, &lv, v, n, x, &doublepluszero, 0, surface);
	else if (!strcmp(*surf_stat, "interpolate/1.approx")) {
		lowesb_(x, y, weights, diagonal, 1, iv, &liv, &lv, v);
		lowese_(iv, &liv, &lv, v, *n, x, surface);
		nsing = iv[29];
		for(i = 0; i < *n; i++) *trL = *trL + diagonal[i];
		lowesa_(trL, n, d, &tau, &nsing, one_delta, two_delta);
		loess_prune(parameter, a, xi, vert, vval);
	}
    else if (!strcmp(*surf_stat, "interpolate/2.approx")) {
		lowesb_(x, y, robust, &doublepluszero, 0, iv, &liv, &lv, v);
		lowese_(iv, &liv, &lv, v, *n, x, surface);
		nsing = iv[29];
		ehg196_(tau, *d, *span, trL);
		lowesa_(trL, n, d, &tau, &nsing, one_delta, two_delta);
		loess_prune(parameter, a, xi, vert, vval);
	}
	else if (!strcmp(*surf_stat, "direct/approximate")) {
		lowesf_(x, y, weights, iv, &liv, &lv, v, n, x, diagonal, 1, surface);
		nsing = iv[29];
		for(i = 0; i < (*n); i++) *trL = *trL + diagonal[i];
		lowesa_(trL, n, d, &tau, &nsing, one_delta, two_delta);
	}
	else if (!strcmp(*surf_stat, "interpolate/exact")) {
		hat_matrix = Calloc((*n)*(*n), double);
		LL = Calloc((*n)*(*n), double);
		lowesb_(x, y, weights, diagonal, 1, iv, &liv, &lv, v);
		lowesl_(iv, &liv, &lv, v, n, x, hat_matrix);
		lowesc_(n, hat_matrix, LL, trL, one_delta, two_delta);
		lowese_(iv, &liv, &lv, v, *n, x, surface);
		loess_prune(parameter, a, xi, vert, vval);
		free(hat_matrix);
		free(LL);
	}
	else if (!strcmp(*surf_stat, "direct/exact")) {
		hat_matrix = Calloc((*n)*(*n), double);
		LL = Calloc((*n)*(*n), double);
		//lowesf_(x, y, weights, iv, liv, lv, v, n, x, hat_matrix, &two, surface);//seems wrong.
		lowesf_(x, y, weights, iv, &liv, &lv, v, n, x, hat_matrix, 2, surface);
		lowesc_(n, hat_matrix, LL, trL, one_delta, two_delta);
        k = (*n) + 1;
		for(i = 0; i < (*n); i++)
			diagonal[i] = hat_matrix[i * k];
		free(hat_matrix);
		free(LL);
	}
	loess_free();
}

static void loess_(double *y, double *x_, long *size_info, double *weights,
            double *span, long  *degree, long  *parametric, long *drop_square,
            long  *normalize, char	**statistics, char **surface, double *cell,
            char **trace_hat_in, long *iterations, double*fitted_values,
            double *fitted_residuals, double *enp, double *s, double *one_delta,
            double *two_delta, double *pseudovalues, double*trace_hat_out, double *diagonal,
            double *robust, double *divisor, long  *parameter, long  *a,
            double *xi, double *vert, double *vval){
	double new_cell, trL, delta1, delta2, sum_squares = 0, *pseudo_resid=NULL,
                trL_tmp = 0, d1_tmp = 0, d2_tmp = 0, sum, mean;
	long	i, j, k, p, N, D, sum_drop_sqr = 0, sum_parametric = 0, setLf,	
                nonparametric = 0, zero = 0, max_kd;
	char   *new_stat;

	D = size_info[0];
	N = size_info[1];
	max_kd = (N > 200 ? N : 200);
	*one_delta = *two_delta = *trace_hat_out = 0;

	double x[D * N], x_tmp[D*N], temp[N], xi_tmp[max_kd], vert_tmp[D * 2], 
                vval_tmp[(D + 1) * max_kd], diag_tmp[N];
	long a_tmp[max_kd], param_tmp[N], order_parametric[D], order_drop_sqr[D];
    integer int_temp[N];//original code sent double, but lowesw & lowesp want an int

    if((*iterations) > 0)
        pseudo_resid =  malloc(N * sizeof(double));

	new_cell = (*span) * (*cell);
	for(i = 0; i < N; i++)
		robust[i] = 1;
        for(i = 0; i < (N * D); i++)
                x_tmp[i] = x_[i];
	if((*normalize) && (D > 1)) {
		int cut = ceil(0.100000000000000000001 * N);
		for(i = 0; i < D; i++) {
			k = i * N;
			for(j = 0; j < N; j++)
				temp[j] = x_[k + j];
			qsort(temp, N, sizeof(double), comp);
			sum = 0;
			for(j = cut; j <= (N - cut - 1); j++)
			        sum = sum + temp[j];
			mean = sum / (N - 2 * cut);
			sum = 0;
			for(j = cut; j <= (N - cut - 1); j++) {
				temp[j] = temp[j] - mean;
				sum = sum + temp[j] * temp[j];
			}
			divisor[i] = sqrt(sum / (N - 2 * cut - 1));
			for(j = 0; j < N; j++) {
				p = k + j;
				x_tmp[p] = x_[p] / divisor[i];		
			}
		}
	}
	else
		for(i = 0; i < D; i++) divisor[i] = 1;
	j = D - 1;
	for(i = 0; i < D; i++) {
		sum_drop_sqr = sum_drop_sqr + drop_square[i];
		sum_parametric = sum_parametric + parametric[i];
		if(parametric[i])
			order_parametric[j--] = i;
		else
			order_parametric[nonparametric++] = i;
	}
    for(i = 0; i < D; i++) {
        order_drop_sqr[i] = 2 - drop_square[order_parametric[i]];
        k = i * N;
        p = order_parametric[i] * N;
        for(j = 0; j < N; j++)
            x[k + j] = x_tmp[p + j];
    }
    Apop_assert_n(!((*degree) == 1 && sum_drop_sqr), 
                "Specified the square of a factor predictor to be dropped when degree = 1");
	Apop_assert_n(!(D == 1 && sum_drop_sqr), 
                "Specified the square of a predictor to be dropped with only one numeric predictor");
	Apop_assert_n(sum_parametric != D, "Specified parametric for all predictors");
	for(j = 0; j <= (*iterations); j++) {
		new_stat = j ? "none" : *statistics;
		for(i = 0; i < N; i++)
			robust[i] = weights[i] * robust[i];
		condition(surface, new_stat, trace_hat_in);
		setLf = !strcmp(surf_stat, "interpolate/exact");
		loess_raw(y, x, weights, robust, &D, &N, span, degree, &nonparametric, order_drop_sqr, 
                &sum_drop_sqr, &new_cell, &surf_stat, fitted_values, parameter, a,
                xi, vert, vval, diagonal, &trL, &delta1, &delta2, &setLf);
		if(j == 0) {
			*trace_hat_out = trL;
			*one_delta = delta1;
			*two_delta = delta2;
		}
		for(i = 0; i < N; i++)
			fitted_residuals[i] = y[i] - fitted_values[i];
		if(j < (*iterations))
			lowesw_(fitted_residuals, &N, robust, int_temp);
	}
	if((*iterations) > 0) {
		pseudovals(N, y, fitted_values, weights, robust, int_temp, pseudovalues);
		
        //BK: I believe that temp here does not rely on prior temp
		loess_raw(pseudovalues, x, weights, weights, &D, &N, span, degree, &nonparametric, 
            order_drop_sqr, &sum_drop_sqr, &new_cell, &surf_stat, temp, param_tmp, a_tmp, 
            xi_tmp, vert_tmp, vval_tmp, diag_tmp, &trL_tmp, &d1_tmp, &d2_tmp, &zero);
		for(i = 0; i < N; i++)
			pseudo_resid[i] = pseudovalues[i] - temp[i];
	}
	if(*iterations == 0)
		for(i = 0; i < N; i++)
			sum_squares = sum_squares + weights[i] * fitted_residuals[i] * fitted_residuals[i];
	else
		for(i = 0; i < N; i++)
			sum_squares = sum_squares + weights[i] * pseudo_resid[i] * pseudo_resid[i];
	*enp = (*one_delta) + 2 * (*trace_hat_out) - N;
	*s = sqrt(sum_squares / (*one_delta));

    if((*iterations) > 0)
        free(pseudo_resid);
}

void loess( struct	loess_struct	*lo) {
	long size_info[2] = {lo->in.p, lo->in.n};
	long iterations = (!strcmp(lo->model.family, "gaussian"))
                        ? 0
                        : lo->control.iterations;		
    if(!strcmp(lo->control.trace_hat, "wait.to.decide")) {
        if(!strcmp(lo->control.surface, "interpolate"))
            lo->control.trace_hat = (lo->in.n < 500) ? "exact" : "approximate";
        else
            lo->control.trace_hat = "exact";
    }
	loess_(lo->in.y, lo->in.x, size_info, lo->in.weights,
		&lo->model.span,
		&lo->model.degree,
		lo->model.parametric,
		lo->model.drop_square,
		&lo->model.normalize,
		&lo->control.statistics,
		&lo->control.surface,
		&lo->control.cell,
		&lo->control.trace_hat,
		&iterations,
		lo->out.fitted_values,
		lo->out.fitted_residuals,
		&lo->out.enp,
		&lo->out.s,
		&lo->out.one_delta,
		&lo->out.two_delta,
		lo->out.pseudovalues,
		&lo->out.trace_hat,
		lo->out.diagonal,
		lo->out.robust,
		lo->out.divisor,
		lo->kd_tree.parameter,
		lo->kd_tree.a,
		lo->kd_tree.xi,
		lo->kd_tree.vert,
		lo->kd_tree.vval);
}	

void loess_free_mem(struct loess_struct *lo) {
    free(lo->in.x);
    free(lo->in.y);
    free(lo->in.weights);
    free(lo->out.fitted_values);
    free(lo->out.fitted_residuals);
    free(lo->out.pseudovalues);
    free(lo->out.diagonal);
    free(lo->out.robust);
    free(lo->out.divisor);
    free(lo->kd_tree.parameter);
    free(lo->kd_tree.a);
    free(lo->kd_tree.xi);
    free(lo->kd_tree.vert);
    free(lo->kd_tree.vval);
}

void loess_summary(struct loess_struct lo, FILE *ap) {
    fprintf(ap, "Number of Observations: %ld\n", lo.in.n);
    fprintf(ap, "Equivalent Number of Parameters: %.1f\n", lo.out.enp);
    if(!strcmp(lo.model.family, "gaussian"))
        fprintf(ap, "Residual Standard Error: ");
    else
        fprintf(ap, "Residual Scale Estimate: ");
    fprintf(ap, "%.4f\n", lo.out.s);
}

//misc.c ---anova and support fns

/* Incomplete beta function.
 * Reference:  Abramowitz and Stegun, 26.5.8.
 * Assumptions: 0 <= x <= 1; a,b > 0.
 */
#define DOUBLE_EPS      2.2204460492503131E-16
#define IBETA_LARGE     1.0e30
#define IBETA_SMALL     1.0e-30

double ibeta(double x, double a, double b) {
    int flipped = 0, i, k, count;
    double I, temp, pn[6], ak, bk, next, prev, factor, val;

    if (x <= 0)
        return(0);
    if (x >= 1)
        return(1);

    /* use ibeta(x,a,b) = 1-ibeta(1-x,b,a) */
    if ((a+b+1)*x > (a+1)) {
        flipped = 1;
        temp = a;
        a = b;
        b = temp;
        x = 1 - x;
    }

    pn[0] = 0.0;
    pn[2] = pn[3] = pn[1] = 1.0;
    count = 1;
    val = x/(1.0-x);
    bk = 1.0;
    next = 1.0;
    do {
        count++;
        k = count/2;
        prev = next;
        if (count%2 == 0)
            ak = -((a+k-1.0)*(b-k)*val)/ ((a+2.0*k-2.0)*(a+2.0*k-1.0));
        else
            ak = ((a+b+k-1.0)*k*val)/ ((a+2.0*k)*(a+2.0*k-1.0));
        pn[4] = bk*pn[2] + ak*pn[0];
        pn[5] = bk*pn[3] + ak*pn[1];
        next = pn[4] / pn[5];
        for (i=0; i<=3; i++)
            pn[i] = pn[i+2];
        if (fabs(pn[4]) >= IBETA_LARGE)
            for (i=0; i<=3; i++)
                    pn[i] /= IBETA_LARGE;
        if (fabs(pn[4]) <= IBETA_SMALL)
            for (i=0; i<=3; i++)
                    pn[i] /= IBETA_SMALL;
    } while (fabs(next-prev) > DOUBLE_EPS*prev);
    factor = a*log(x) + (b-1)*log(1-x);
    factor -= tgammal(a+1) + tgammal(b) - tgammal(a+b);
    I = exp(factor) * next;
    return(flipped ? 1-I : I);
}

/* For comparing two loess models, which is not as useful as one would hope...
double pf(double q, double df1, double df2) {
	return ibeta(q*df1/(df2+q*df1), df1/2, df2/2);
}

struct anova_struct {
	double	dfn;
	double	dfd;
	double  F_value;
	double  Pr_F;
};

void anova(struct loess_struct *one, struct loess_struct *two, struct anova_struct *out){
	double	one_d1, one_d2, one_s, two_d1, two_d2, two_s, rssdiff, d1diff, tmp;
	int     max_enp;
  
	one_d1 = one->out.one_delta;
	one_d2 = one->out.two_delta;
	one_s = one->out.s;
	two_d1 = two->out.one_delta;
	two_d2 = two->out.two_delta;
	two_s = two->out.s;

        rssdiff = fabs(one_s * one_s * one_d1 - two_s * two_s * two_d1);
        d1diff = fabs(one_d1 - two_d1);
        out->dfn = d1diff * d1diff / fabs(one_d2 - two_d2);
	max_enp = (one->out.enp > two->out.enp);
	tmp = max_enp ? one_d1 : two_d1;
        out->dfd = tmp * tmp / (max_enp ? one_d2 : two_d2);
	tmp = max_enp ? one_s : two_s;
        out->F_value = (rssdiff / d1diff) / (tmp * tmp);
        out->Pr_F = 1 - pf(out->F_value, out->dfn, out->dfd);
}
*/


/*
 * Rational approximation to inverse Gaussian distribution.
 * Absolute error is bounded by 4.5e-4.
 * Reference: Abramowitz and Stegun, page 933.
 * Assumption: 0 < p < 1.
 */

static double num[] = {
        2.515517,
        0.802853,
        0.010328
};

static double den[] = {
        1.000000,
        1.432788,
        0.189269,
        0.001308
};

double invigauss_quick(double p) {
  int lower;
  double t, n, d, q;
    if(p == 0.5)
        return(0);
    lower = p < 0.5;
    p = lower ? p : 1 - p;
    t = sqrt(-2 * log(p));
    n = (num[2]*t + num[1])*t + num[0];
    d = ((den[3]*t + den[2])*t + den[1])*t + den[0];
    q = lower ? n/d - t : t - n/d;
    return(q);
}

/*
 * Quick approximation to inverse incomplete beta function,
 * by matching first two moments with the Gaussian distribution.
 * Assumption: 0 < p < 1, a,b > 0.
 */

static double invibeta_quick(double p, double a, double b) {
  double x, m, s, invigauss_quick();
    x = a + b;
    m = a / x;
    s = sqrt((a*b) / (x*x*(x+1)));
    return(GSL_MAX(0.0, GSL_MAX(1.0, invigauss_quick(p)*s + m)));
}

/*
 * Inverse incomplete beta function.
 * Assumption: 0 <= p <= 1, a,b > 0.
 */

static double invibeta(double p,double  a,double  b) {
    int i;
    double ql, qr, qm, qdiff;
    double pl, pr, pm, pdiff;

/*    MEANINGFUL(qm);*/
	qm = 0;
    if(p == 0 || p == 1)
        return(p);

    /* initialize [ql,qr] containing the root */
    ql = qr = invibeta_quick(p, a, b);
    pl = pr = ibeta(ql, a, b);
    if(pl == p)
        return(ql);
    if(pl < p)
        while(1) {
            qr += 0.05;
            if(qr >= 1) {
                pr = qr = 1;
                break;
            }
            pr = ibeta(qr, a, b);
            if(pr == p)
                return(pr);
            if(pr > p)
                break;
        }
    else
        while(1) {
            ql -= 0.05;
            if(ql <= 0) {
                pl = ql = 0;
                break;
            }
            pl = ibeta(ql, a, b);
            if(pl == p)
                return(pl);
            if(pl < p)
                break;
        }

    /* a few steps of bisection */
    for(i = 0; i < 5; i++) {
        qm = (ql + qr) / 2;
        pm = ibeta(qm, a, b);
        qdiff = qr - ql;
        pdiff = pm - p;
        if(fabs(qdiff) < DOUBLE_EPS*qm || fabs(pdiff) < DOUBLE_EPS)
            return(qm);
        if(pdiff < 0) {
            ql = qm;
            pl = pm;
        } else {
            qr = qm;
            pr = pm;
        }
    }

    /* a few steps of secant */
    for(i = 0; i < 40; i++) {
        qm = ql + (p-pl)*(qr-ql)/(pr-pl);
        pm = ibeta(qm, a, b);
        qdiff = qr - ql;
        pdiff = pm - p;
        if(fabs(qdiff) < 2*DOUBLE_EPS*qm || fabs(pdiff) < 2*DOUBLE_EPS)
            return(qm);
        if(pdiff < 0) {
            ql = qm;
            pl = pm;
        } else {
            qr = qm;
            pr = pm;
        }
    }

    /* no convergence */
    return(qm);
}

static double qt(double p, double df) {
  double t = invibeta(fabs(2*p-1), 0.5, df/2);
    return((p>0.5?1:-1) * sqrt(t*df/(1-t)));
}



////The apop_model front end.

void loess(struct loess_struct *lo) ;
void loess_setup( double  *x, double *y, long n, long p, struct  loess_struct *lo) ;


Apop_settings_copy(apop_loess, )
Apop_settings_free(apop_loess, loess_free_mem(&(in->lo_s));)

void matrix_to_FORTRAN(gsl_matrix *inmatrix, double *outFORTRAN, int start_col){
    double *current_outcol = outFORTRAN; 
    for (int i=start_col; i< inmatrix->size2; i++){
        Apop_matrix_col(inmatrix, i, col);
        for (int j=0; j< col->size; j++)
            current_outcol[j]=gsl_vector_get(col,j);
        current_outcol += col->size;
    }
}

#define lo_set(var, dflt) .var = in.lo_s.var ? in.lo_s.var : dflt
#define apop_strcmp(a, b) (((a)&&(b) && !strcmp((a), (b))) || (!(a) && !(b)))

Apop_settings_init(apop_loess,
    Apop_assert(in.data, "I need a .data element to allocate apop_loess_settings.");
    int n = in.data->matrix->size1;
	int	max_kd = n > 200 ? n : 200;
    int p =  (in.data->vector)
                ? in.data->matrix->size2
                : in.data->matrix->size2-1;
    out->lo_s = (struct loess_struct){
        .in.n = n,
        .in.p = p,
	    .in.y =  malloc(n * sizeof(double)),
        .in.x =  malloc(n * p * sizeof(double)),

	    lo_set(model.degree , 2),
	    lo_set(model.normalize , 'y'),
        .model.family = apop_strcmp(in.lo_s.model.family , "symmetric") ? "symmetric": "gaussian",
        lo_set(model.span , 0.75),
        //.model.span = in.span ? in.span : 0.75,

        .control.surface = apop_strcmp(in.lo_s.control.surface , "direct") ? "direct" : "interpolate",
        .control.statistics = apop_strcmp(in.lo_s.control.statistics , "exact") ? "exact" : "approximate",
        lo_set(control.cell , 0.2),
        .control.trace_hat = apop_strcmp(in.lo_s.control.trace_hat , "exact") ? "exact" 
                        : apop_strcmp(in.lo_s.control.trace_hat , "approximate") ? "approximate" 
                        : "wait.to.decide",
        lo_set(control.iterations, (apop_strcmp(in.lo_s.model.family , "symmetric") ? 4 : 0)),

        .out.fitted_values =  malloc(n * sizeof(double)),
        .out.fitted_residuals =  malloc(n * sizeof(double)),
        .out.pseudovalues =  malloc(n * sizeof(double)),
        .out.diagonal =  malloc(n * sizeof(double)),
        .out.robust =  malloc(n * sizeof(double)),
        .out.divisor =  malloc(p * sizeof(double)),

        .kd_tree.parameter =  malloc(7 * sizeof(long)),
        .kd_tree.a =  malloc(max_kd * sizeof(long)),
        .kd_tree.xi =  malloc(max_kd * sizeof(double)),
        .kd_tree.vert =  malloc(p * 2 * sizeof(double)),
        .kd_tree.vval =  malloc((p + 1) * max_kd * sizeof(double)),
    };
    Apop_varad_set(ci_level, 0.95);
    struct loess_struct *lo = &(out->lo_s);
    if (in.data->weights)
        lo->in.weights = in.data->weights->data;
    else {
        lo->in.weights = malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) 
            lo->in.weights[i] = 1;
    }
    int startat = 0;
    if (in.data->vector) //OK, then that's the dependent var.
        memcpy(lo->in.y, in.data->vector->data, n*sizeof(double));
    else {  //use the first col as the dep. var.
        startat =1;
        Apop_col_v(in.data, 0, col);
        for (int j=0; j< col->size; j++)
            lo->in.y[j] = gsl_vector_get(col,j);
    }
    matrix_to_FORTRAN(in.data->matrix, lo->in.x, startat);
	for(int i = 0; i < 8; i++)
        lo->model.parametric[i] = lo->model.drop_square[i] = 0;
)

apop_data * loess_predict (apop_data *in, apop_model *m){
    //Massage inputs to FORTRAN's format
    double *eval_here = malloc(sizeof(double)*in->matrix->size1*(in->matrix->size2-1));
    matrix_to_FORTRAN(in->matrix, eval_here, 1);
    int want_cov = Apop_settings_get(m, apop_loess, want_predict_ci)=='y';
    struct pred_struct pred = (struct pred_struct){ };

    predict(eval_here, in->matrix->size1, &(Apop_settings_get(m, apop_loess, lo_s)), &pred, want_cov);

    //Massage FORTRAN's output to Apophenia's formats
    Apop_col_v(in, 0, firstcol)
    gsl_vector_view v = gsl_vector_view_array(pred.fit, firstcol->size);
    gsl_vector_memcpy(firstcol, &(v.vector));
    apop_data *ci =  apop_data_add_page(in, apop_data_alloc(in->matrix->size1, 3), "<Confidence>");
    apop_name_add(ci->names, "standard error", 'c');
    apop_name_add(ci->names, "lower CI", 'c');
    apop_name_add(ci->names, "upper CI", 'c');

    if (want_cov){
        //Find confidence intervals. Used to be in loess's pointwise().
        double coverage = Apop_settings_get(m, apop_loess, ci_level);
        double t_dist = qt(1 - (1 - coverage)/2, pred.df);
        for(int i = 0; i < in->matrix->size1; i++) {
            double limit = pred.se_fit[i] * t_dist;
            apop_data_set(ci, i, 0, limit);
            apop_data_set(ci, i, 1, pred.fit[i] + limit);
            apop_data_set(ci, i, 2, pred.fit[i] - limit);
        }	
    }

    free(eval_here);
    pred_free_mem(&pred);
    return in;
}

static double onerow(gsl_vector *v, void *sd){ 
    return log(gsl_ran_gaussian_pdf(v->data[2], *((double*)sd))); 
}

//Assumes one gaussian, unweighted.
static long double loess_ll(apop_data *d, apop_model *m){
    apop_data *exp = apop_data_get_page(d, "<Predicted>");
    Apop_col_tv(exp, "residual", residuals);
    double sd = sqrt(apop_vector_var(residuals));
    return apop_map_sum(exp, .param=&sd, .part='r', .fn_vp= onerow);
}

static void apop_loess_est(apop_data *d, apop_model *out){
    if (!Apop_settings_get_group(out, apop_loess))
        Apop_model_add_group(out, apop_loess, .data=d);
    out->data = d;
    loess(&Apop_settings_get(out, apop_loess, lo_s));

    //setup the expected matrix. In a perfect world, this wouldn't all be cut/pasted from apop_OLS.
    //Also, it wouldn't be 14 lines.
    apop_data *expect = apop_data_add_page(out->info, apop_data_alloc(d->matrix->size1, 3), "<Predicted>");
    if (!out->info) out->info = expect;
    apop_name_add(expect->names, (out->data->names->colct ? out->data->names->col[0] : "Expected"), 'c');
    apop_name_add(expect->names, "Predicted", 'c');
    apop_name_add(expect->names, "Residual", 'c');
    gsl_vector *v = gsl_vector_alloc(d->matrix->size1);
    double *holding = v->data;
    v->data = Apop_settings_get(out, apop_loess, lo_s.in.y);
    Apop_col_v(expect, 0, y_data)
    gsl_vector_memcpy(y_data, v);
    v->data = Apop_settings_get(out, apop_loess, lo_s.out.fitted_values);
    Apop_col_v(expect, 1, fitted)
    gsl_vector_memcpy(fitted, v);
    v->data = Apop_settings_get(out, apop_loess, lo_s.out.fitted_residuals);
    Apop_col_v(expect, 2, resid)
    gsl_vector_memcpy(resid, v);
    v->data = holding;
    gsl_vector_free(v);
}

static void apop_loess_print(apop_model *in, FILE *out){
    loess_summary(Apop_settings_get(in, apop_loess, lo_s), out);
}

static void loess_prep(apop_data *data, apop_model *params){
    apop_predict_vtable_add(loess_predict, apop_loess);
    apop_model_print_vtable_add(apop_loess_print, apop_loess);
    apop_model_clear(data, params);
}

apop_model *apop_loess = &(apop_model){.name="Loess smoothing", .vsize = -1, .dsize=1, 
        .estimate =apop_loess_est, .log_likelihood = loess_ll, .prep = loess_prep};
