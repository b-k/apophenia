
/** \file apop_linear_constraint.c 
  \c apop_linear_constraint finds a point that meets a set of linear constraints. This takes a lot of machinery, so it gets its own file.

Copyright (c) 2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
*/
#include "apop_internal.h"

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

static int binds(gsl_vector const *v, double k, gsl_vector const *b, double margin){
    double d;
    gsl_blas_ddot(v, b, &d);
    return d < k + margin;
}

static double trig_bit(gsl_vector *dimv, gsl_vector *otherv, double off_by){
    double theta, costheta, dot, out;
    gsl_blas_ddot(dimv, otherv, &dot);
    costheta = dot/(magnitude(dimv)*magnitude(otherv));
    theta = acos(costheta);
    out = off_by/gsl_pow_2(sin(theta)); 
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
static void get_candiate(gsl_vector *beta, apop_data *constraint, int current, gsl_vector *candidate, double margin){
    double k, ck, off_by, s;
    gsl_vector *pseudobeta        = NULL;
    gsl_vector *pseudocandidate   = NULL;
    gsl_vector *pseudocandidate2  = NULL;
    gsl_vector *fix               = NULL;
    Apop_row_v(constraint, current, cc);
    ck = gsl_vector_get(constraint->vector, current);
    find_nearest_point(beta, ck, cc, candidate);
    for (size_t i=0; i< constraint->vector->size; i++){
        if (i!=current){
            Apop_row_v(constraint, i, other);
            k   =apop_data_get(constraint, i, -1);
            if (binds(candidate, k, other, margin)){
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

/** This is designed to be called from within the constraint method of your \ref
apop_model. Just write the constraint vector+matrix and this will do the rest.
See the outline page for detailed discussion on setting contrasts. 
 
\param beta    The proposed vector about to be tested. No default, must not be \c NULL.

\param constraint  
A vector/matrix pair [v | m1 m2 ... mn] where each row is interpreted as a less-than inequality:
\f$v < m1x1+ m2x2 + ... + mnxn\f$.  For example, say your constraints are 
\f$3 < 2x + 4y - 7z\f$ and \f$y\f$ is positive, i.e. \f$0 < y\f$.
Allocate and fill the matrix representing these two constraints via:
\code
apop_data *constr = apop_data_falloc((2,2,3), 3,  2, 4, 7,
                                              0,  0, 1, 0);
\endcode
. Default: each elements is greater than zero. E.g., for three parameters:
\code
apop_data *constr = apop_data_falloc((3,3,3), 0,  1, 0, 0,
                                              0,  0, 1, 0,
                                              0,  0, 0, 1);
\endcode

\param margin If zero, then this is a >= constraint, otherwise I will return a point this amount within the borders. You could try \c GSL_DBL_EPSILON, which is the smallest value a \c double can hold, or something like 1e-3. Default = 0.

return The penalty = the distance between beta and the closest point that meets the constraints.
If the constraint is not met, this \c beta is shifted by \c margin (Euclidean distance) to meet the constraints. 

\li If your \ref apop_data is not just a vector, try \ref apop_data_pack to pack it into a vector. This is what \ref apop_maximum_likelihood does.

\li This function uses the \ref designated syntax for inputs.
todo The apop_linear_constraint function doesn't check for odd cases like coplanar constraints.
*/
#ifdef APOP_NO_VARIADIC
long double apop_linear_constraint(gsl_vector *beta, apop_data * constraint, double margin){
#else
apop_varad_head(long double, apop_linear_constraint){
    static threadlocal apop_data *default_constraint;
    gsl_vector * apop_varad_var(beta, NULL);
    double apop_varad_var(margin, 0);
    apop_data * apop_varad_var(constraint, NULL);
    Apop_assert(beta, "The vector to be checked is NULL.");
    if (!constraint){
        if (default_constraint && beta->size != default_constraint->vector->size){
            apop_data_free(default_constraint);
            default_constraint = NULL;
        }
        if (!default_constraint){
            default_constraint = apop_data_alloc(0,beta->size, beta->size);
            default_constraint->vector = gsl_vector_calloc(beta->size);
            gsl_matrix_set_identity(default_constraint->matrix);
        }
        constraint = default_constraint;
    }
    return apop_linear_constraint_base(beta, constraint, margin);
}

 long double apop_linear_constraint_base(gsl_vector *beta, apop_data * constraint, double margin){
#endif
    static threadlocal gsl_vector *closest_pt = NULL;
    static threadlocal gsl_vector *candidate  = NULL;
    static threadlocal gsl_vector *fix        = NULL;
    int constraint_ct = constraint->matrix->size1;
    int bindlist[constraint_ct];
    int i, bound = 0;
    /* For added efficiency, keep a scratch vector or two on hand. */
    if (closest_pt==NULL || closest_pt->size != constraint->matrix->size2){
        closest_pt  = gsl_vector_calloc(beta->size);
        candidate   = gsl_vector_alloc(beta->size);
        fix         = gsl_vector_alloc(beta->size);
        closest_pt->data[0] = GSL_NEGINF;
    }
    /* Do any constraints bind?*/
    memset(bindlist, 0, sizeof(int)*constraint_ct);
    for (i=0; i< constraint_ct; i++){
        Apop_row_v(constraint, i, c);
        bound           +=
        bindlist[i]      = binds(beta, apop_data_get(constraint, i, -1), c, margin);
    }
    if (!bound) return 0;   //All constraints met.
    gsl_vector *base_beta = apop_vector_copy(beta);
    /* With only one constraint, it's easy. */
    if (constraint->vector->size==1){
        Apop_row_v(constraint, 0, c);
        find_nearest_point(base_beta, constraint->vector->data[0], c, beta);
        goto add_margin;
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
    for (i=0; i< constraint_ct; i++)
        if (bindlist[i]){
            get_candiate(base_beta, constraint, i, candidate, margin);
            if(apop_vector_distance(base_beta, candidate) < apop_vector_distance(base_beta, closest_pt))
                gsl_vector_memcpy(closest_pt, candidate);
        }
    gsl_vector_memcpy(beta, closest_pt);
add_margin:
    for (i=0; i< constraint_ct; i++){
        if(bindlist[i]){
            Apop_row_v(constraint, i, c);
            gsl_vector_memcpy(fix, c);
            gsl_vector_scale(fix, magnitude(fix));
            gsl_vector_scale(fix, margin);
            gsl_vector_add(beta, fix);
        }
    }
    long double out = apop_vector_distance(base_beta, beta);
    gsl_vector_free(base_beta);
    return out;
}
