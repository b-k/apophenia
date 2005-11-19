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

#include "wnlib.h"


void wn_copy_mat(double **out_mat,double **in_mat,int len_i,int len_j)
{
  int i;

  for(i=0;i<len_i;++i)
  {
    wn_copy_vect(out_mat[i],in_mat[i],len_j);
  }
}


/**************************************************************************

wn_random_permutation(permutation,size)
wn_identity_permutation(permutation,size)

wn_permute_permutation(result,perm1,perm2,size)
wn_invert_permutation(inverse,permutation,size)

bool wn_is_valid_permutation(permutation,size)

**************************************************************************/

void wn_identity_permutation(int permutation[],int size)
{
  int i;

  for(i=0;i<size;i++)
  {
    permutation[i] = i;
  }
}

static gsl_rng *r;


void wn_random_permutation(int permutation[],int size)
{
  int *remaining_numbers,top_of_remaining_numbers,selection,count;

	r = gsl_rng_alloc(gsl_rng_taus);

  wn_gpmake("no_free");

  remaining_numbers = (int *)wn_zalloc(size*sizeof(int));
  wn_identity_permutation(remaining_numbers,size);

  for(count=0,top_of_remaining_numbers=size;
      count<size;
      count++,top_of_remaining_numbers--)
  {
    selection = (int) gsl_rng_uniform(r) * top_of_remaining_numbers;

    permutation[count] = remaining_numbers[selection];
    remaining_numbers[selection] = 
                          remaining_numbers[top_of_remaining_numbers-1];
  }

  wn_gpfree();
}


void wn_permute_permutation(int result[],int perm1[],int perm2[],int size)
{
  wn_permute_array(result,perm1,perm2,size);
}


void wn_invert_permutation(int inverse[],int permutation[],int size)
{
  int i;

  wn_assert(inverse != permutation);

  for(i=0;i<size;i++)
  {
    inverse[permutation[i]] = i;
  }
}


bool wn_is_valid_permutation(int permutation[],int size)
{
  int i,entry;
  int *counts;
  bool ret;

  ret = TRUE;

  wn_gpmake("no_free");

  counts = (int *)wn_zalloc(size*sizeof(bool));

  for(i=0;i<size;i++)
  {
    entry = permutation[i];

    if(!((0<=entry)&&(entry<size)))
    {
      ret = FALSE;

      break;
    }

    counts[permutation[i]]++;
  }

  if(ret == TRUE)
  {
    for(i=0;i<size;i++)
    {
      if(counts[i] != 1)
      {
        ret = FALSE;
  
        break;
      }
    }
  }

  wn_gpfree();

  return(ret);
}


/**********************************************************************

wn_zero_mat(mat,len_i,len_j)

wn_identity_mat(mat,len_i)
wn_hilbert_mat(mat,len_i)

**********************************************************************/
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



void wn_zero_mat(double **mat,int len_i,int len_j)
{
  int i;

  for(i=0;i<len_i;++i)
  {
    wn_zero_vect(mat[i],len_j);
  }
}


void wn_identity_mat(double **mat,int len)
{
  int i;

  wn_zero_mat(mat,len,len);

  for(i=0;i<len;++i)
  {
    mat[i][i] = 1.0;
  }
}


void wn_hilbert_mat(double **mat,int len)
{
  int i,j;

  for(i=0;i<len;++i)  
  for(j=0;j<len;++j)
  {
    mat[i][j] = 1.0/((double)(i+j+1));
  }
}



/**********************************************************************

wn_make_mat(&mat,len_i,len_j)
wn_free_mat(mat,len_i,len_j)

**********************************************************************/
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




void wn_make_mat(double ***pmat,int len_i,int len_j)
{
  int i;

  *pmat = (double **)wn_alloc(len_i*sizeof(double *));

  for(i=0;i<len_i;i++)
  {
    wn_make_vect(&((*pmat)[i]),len_j);
  }
}


void wn_free_mat(double **mat,int len_i,int len_j)
{
  int i;

  for(i=0;i<len_i;i++)
  {
    wn_free_vect(mat[i],len_j);
  }

  wn_free(mat);
}




/**********************************************************************

wn_mult_mat_by_vect(result_vect,mat,vect,len_i,len_j)

wn_mult_mats(result_mat,mat1,mat2,len_i,len_j,len_k)

**********************************************************************/
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




void wn_mult_mat_by_vect(double *result_vect,double **mat,double *vect,
			 int len_i,int len_j)
{
  int i;

  for(i=0;i<len_i;i++)
  {
    result_vect[i] = wn_dot_vects(mat[i],vect,len_j);
  }
}


void wn_mult_mats(double **result_mat,double **mat1,double **mat2,
		  int len_i,int len_j,int len_k)
{
  double *vect2;
  int i,j,k;

  wn_gpmake("no_free");
  wn_gplabel("matrix multiply group");

  wn_make_vect(&vect2,len_j);

  for(k=0;k<len_k;k++)
  {
    for(j=0;j<len_j;j++)
    {
      vect2[j] = mat2[j][k];
    }

    for(i=0;i<len_i;i++)
    {
      result_mat[i][k] = wn_dot_vects(mat1[i],vect2,len_j);
    }
  }

  wn_gpfree();
}




/**********************************************************************

wn_invert_mat(&code,mat,len_i)

**********************************************************************/



local int *row_permutation,*var_permutation,
          *row_permutation_inverse,*var_permutation_inverse;
local double **temp_mat;


local void set_up_arrays(double **mat,int len_i)
{
  int var;

  temp_mat = (double **)wn_zalloc(len_i*sizeof(double *));

  row_permutation = (int *)wn_zalloc(len_i*sizeof(int));
  var_permutation = (int *)wn_zalloc(len_i*sizeof(int));
  row_permutation_inverse = (int *)wn_zalloc(len_i*sizeof(int));
  var_permutation_inverse = (int *)wn_zalloc(len_i*sizeof(int));

  wn_random_permutation(row_permutation,len_i);
  memcpy(temp_mat,mat,len_i*sizeof(double *));
  wn_permute_array(mat,row_permutation,temp_mat,len_i);

  for(var=0;var<len_i;var++)
  {
    var_permutation[var] = var+len_i;
  }
}


local void get_best_var_from_pivot_row
(
  int *pcode,
  int *pvar,
  double **mat,
  int len_i,
  int pivot_row
)
{
  int var,best_var;
  double best_a,a;

  best_var = -1;
  best_a = 0.0;

  for(var=0;var<len_i;var++)
  {
    if(var_permutation[var] >= len_i)
    {
      a = mat[pivot_row][var];
      a = wn_abs(a);
  
      if(a > best_a)
      {
        best_a = a;
        best_var = var;
      }
    }
  }

  if(best_var == -1)
  {
    *pcode = WN_SINGULAR;

    return;
  }

  *pvar = best_var;

  *pcode = WN_SUCCESS;
}



/********************************************************************

  This is the slow part of the algorithm.

********************************************************************/
local void eliminate_var_from_all_other_rows
(
  double **mat,
  int len_i,
  int pivot_row,
  int pivot_var
)
{
  int row,col;
  double self_multiplier,multiplier;

  wn_assert(mat[pivot_row][pivot_var] != 0.0);

  self_multiplier = 1.0 / mat[pivot_row][pivot_var];
  
  mat[pivot_row][pivot_var] = 1.0;      /* for new inverted var */

  for(col=0;col<len_i;++col)
  {
    mat[pivot_row][col] *= self_multiplier;
  }

  for(row=0;row<len_i;row++)
  {
    if(row != pivot_row)
    {
      multiplier = -mat[row][pivot_var];

      mat[row][pivot_var] = 0.0;              /* for new inverted var */

      if(multiplier != 0.0)
      {
        wn_add_scaled_vect(mat[row],mat[pivot_row],multiplier,len_i); 
      }
    }
  }
}


local void in_place_invert(int *pcode,double **mat,int len_i)
{
  int pivot_row,pivot_var;

  for(pivot_row=0;pivot_row<len_i;pivot_row++)
  {
    get_best_var_from_pivot_row(pcode,&pivot_var,mat,len_i,pivot_row);
    if(*pcode != WN_SUCCESS)
    {
      break;
    }

    eliminate_var_from_all_other_rows(mat,len_i,pivot_row,pivot_var);

    wn_swap(row_permutation[pivot_row],var_permutation[pivot_var],int);
  }
}


local double *temp_row;

local void apply_permutation_to_all_rows(double **mat,int *permutation,
					 int len_i)
{
  double *row;
  int row_index;

  temp_row = (double *)wn_zalloc(len_i*sizeof(double));

  for(row_index=0;row_index<len_i;row_index++)
  {
    row = mat[row_index];

    memcpy(temp_row,row,len_i*sizeof(double));
    wn_permute_array(row,permutation,temp_row,len_i);
  }
}


local void unscramble_inverse(double **mat,int len_i)
{
  int row;

  for(row=0;row<len_i;row++)
  {
    row_permutation[row] -= len_i;
  }

  wn_invert_permutation(row_permutation_inverse,row_permutation,len_i);
  wn_invert_permutation(var_permutation_inverse,var_permutation,len_i);

  apply_permutation_to_all_rows(mat,var_permutation_inverse,len_i);

  memcpy(temp_mat,mat,len_i*sizeof(double *));
  wn_permute_array(mat,row_permutation_inverse,temp_mat,len_i);
}


void wn_invert_mat(int *pcode,double **mat,int len_i)
{
  wn_gpmake("no_free");
  wn_gplabel("invert matrix group");

  set_up_arrays(mat,len_i);

  in_place_invert(pcode,mat,len_i);

  if(*pcode == WN_SUCCESS)
  {
    unscramble_inverse(mat,len_i);
  }

  wn_gpfree();
}




void wn_transpose_mat(double **to_mat,double **from_mat,int len_i,int len_j)
{
  int i,j;

  for(i=0;i<len_i;++i)  
  for(j=0;j<len_j;++j)
  {
    to_mat[j][i] = from_mat[i][j];
  }
}




