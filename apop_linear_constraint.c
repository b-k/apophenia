/** \file apop_linear_constraint.c 
  \c apop_linear_constraint finds a point that meets a set of linear constraints. This takes a lot of machinery, so it gets its own file.

  (c) 2007 Ben Klemens. Licensed under the GNU GPL v2.
*/
#include <apop.h>

static double magnitude(gsl_vector *v){
 double out;
    gsl_blas_ddot(v, v, &out);
    return out;
}

static void find_nearest_point(gsl_vector *V, double k, gsl_vector *B, gsl_vector *out){
    /* Find X such that BX =K and there is an S such that X + SB=V. */
  double S=0; //S = (BV-K)/B'B.
    gsl_blas_ddot(B, V, &S);
    S   -= k;
assert(!gsl_isnan(S));
    S   /= magnitude(B);
assert(!gsl_isnan(S));
    gsl_vector_memcpy(out, B); //X = -SB +V
    gsl_vector_scale(out, -S);
    gsl_vector_add(out, V);
assert(!gsl_isnan(gsl_vector_get(out,0)));
}

static int binds(gsl_vector *v, int k, gsl_vector *b){
  double d;
    gsl_blas_ddot(v, b, &d);
    return d < k;
}

static double trig_bit(gsl_vector *dimv, gsl_vector *otherv, double off_by){
  double    theta, costheta, dot, out;
    gsl_blas_ddot(dimv, otherv, &dot);
    costheta = dot/(magnitude(dimv)*magnitude(otherv));
    theta   = acos(costheta);
    out     = off_by/gsl_pow_2(sin(theta)); 
    return out;
}

/* The hard part is when your candidate point does not satisfy other
   constraints, so you need to translate the point until it meets the new hypersurface.
   How far is that? Project beta onto the new surface, and find the
   distance between that projection and the original surface. Then
   translate beta toward the original surface by that amount. The
   projection of the translated beta onto the new surface now also touches the old
   surface.
   */
static void get_candiate(gsl_vector *beta, apop_data *constraint, int current, gsl_vector *candidate){
  double    k, ck, off_by, s;
  int       i;
  gsl_vector *pseudobeta        = NULL;
  gsl_vector *pseudocandidate   = NULL;
  gsl_vector *pseudocandidate2  = NULL;
  gsl_vector *fix               = NULL;
  APOP_ROW(constraint, current, cc);
    ck  =gsl_vector_get(constraint->vector, current);
    find_nearest_point(beta, ck, cc, candidate);
    for (i=0; i< constraint->vector->size; i++){
        if (i!=current){
            APOP_ROW(constraint, i, other);
            k   =apop_data_get(constraint, i, -1);
            if (binds(candidate, k, other)){
                if (!pseudobeta){
                    pseudobeta          = gsl_vector_alloc(beta->size);
                    gsl_vector_memcpy(pseudobeta, beta);
                    pseudocandidate     = gsl_vector_alloc(beta->size);
                    pseudocandidate2    = gsl_vector_alloc(beta->size);
                    fix                 = gsl_vector_alloc(beta->size);
                }
                find_nearest_point(pseudobeta, k, other, pseudocandidate);
                find_nearest_point(pseudocandidate, ck, cc, pseudocandidate2);
                off_by  = apop_vector_distance(pseudocandidate, pseudocandidate2);
                s       = trig_bit(cc, other, off_by);
                gsl_vector_memcpy(fix, cc);
                gsl_vector_scale(fix, magnitude(cc));
                gsl_vector_scale(fix, s);
                gsl_vector_add(pseudobeta, fix);
                find_nearest_point(pseudobeta, k, other, candidate);
                gsl_vector_memcpy(pseudobeta, candidate);
            } 
        }
    }
    if (fix){ 
        gsl_vector_free(fix); gsl_vector_free(pseudobeta);
        gsl_vector_free(pseudocandidate); gsl_vector_free(pseudocandidate2);
    }
}

/** This is designed to be called from within your own constraint
 function. Just write the constraint vector and this will do the rest.
 
 \param beta    The proposed vector about to be tested. 
 \param constraint  The constraints. See \ref apop_F_test on writing
 contrasts. To give a quick example, say your constraint is $3 < 2x +
 4y - 7z$; then the first row of your \c data->vector element would be 3, and the
 first row of the \c data->matrix element would be [2 4 -7].
 \param margin If zero, then this is a >= constraint, otherwise I will return a point this amount within the borders. You could try \c GSL_DBL_EPSILON, which is the smallest value a \c double can hold, or something like 1e-3.
 \param returned_beta If the constraint is not met, this is the closest point that meets the constraints (Euclidian distance)
 \return The penalty = the distance between beta and the closest point that meets the constraints.
\todo The apop_linear_constraint function doesn't check for odd cases like coplanar constraints.
 */
double  apop_linear_constraint(gsl_vector *beta, apop_data * constraint, double margin,  gsl_vector *returned_beta){
  static gsl_vector *closest_pt = NULL;
  static gsl_vector *candidate  = NULL;
  static gsl_vector *fix        = NULL;
  int               bindlist[beta->size];
  int               i, bound     = 0;
    /* For added efficiency, keep a scratch vector or two on hand. */
    if (closest_pt==NULL || closest_pt->size != constraint->matrix->size2){
        closest_pt  = gsl_vector_calloc(beta->size);
        candidate   = gsl_vector_alloc(beta->size);
        fix         = gsl_vector_alloc(beta->size);
        closest_pt->data[0] = GSL_NEGINF;
    }
    /* Do any constraints bind?*/
    memset(bindlist, 0, sizeof(int)*beta->size);
    for (i=0; i< constraint->matrix->size1; i++){
        APOP_ROW(constraint, i, c);
        bound           +=
        bindlist[i]      = binds(beta, apop_data_get(constraint, i, -1), c);
    }
    if (!bound)    //All constraints met.
        return 0;
    /* With only one constraint, it's easy. */
    if (constraint->vector->size==1){
        APOP_ROW(constraint, 0, c);
        find_nearest_point(beta, constraint->vector->data[0],c, returned_beta);
        return apop_vector_distance(beta,returned_beta);
    }
    /* Finally, multiple constraints, at least one binding.
       For each surface, pick a candidate point.
       Check whether the point meets the other constraints. 
            if not, translate to a new point that works.
            [Do this by maintaining a pseudopoint that translates by the
            necessary amount.]
        Once you have a candidate point, compare its distance to the
        current favorite; keep the best.
     */
    for (i=0; i< constraint->matrix->size2; i++){
        if (bindlist[i])
            get_candiate(beta, constraint, i, candidate);
        if(apop_vector_distance(beta, candidate) < apop_vector_distance(beta, closest_pt))
            gsl_vector_memcpy(closest_pt, candidate);
    }
    gsl_vector_memcpy(returned_beta, closest_pt);
    for (i=0; i< constraint->matrix->size1; i++){
        if(bindlist[i]){
            APOP_ROW(constraint, i, c);
            gsl_vector_memcpy(fix, c);
            gsl_vector_scale(fix, magnitude(fix));
            gsl_vector_scale(fix, margin);
            gsl_vector_add(returned_beta, fix);
        }
    }
    return apop_vector_distance(beta, returned_beta);
}
