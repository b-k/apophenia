/**********************************************************************

wn_make_vect(&vect,len)
wn_free_vect(vect,len)

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

#include "wnlib.h"


void wn_make_vect(double **pvect,int len)
{
  *pvect = (double *)wn_alloc(len*sizeof(double));

  wn_zero_vect(*pvect,len);
}


/*ARGSUSED*/ void wn_free_vect(double *vect,int len)
{
  wn_free((ptr)vect);
}



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

void wn_memzero(ptr out,int n)
{
	  memset(out,'\0',n);
}



void wn_copy_vect(double to_vect[],double from_vect[],int len)
{
  memcpy((ptr)to_vect,(ptr)from_vect,len*sizeof(double));
}



/*****************************************************************

  This code is the inner loop of various n^3 algorithms, so it must
  be fast.  That is why this code is so ugly.

*****************************************************************/
double wn_dot_vects
(
  register double *vect1,
  register double *vect2,
  int len
)
{
  register double result;

  if(len <= 0)
  {
    return(0.0);
  }

  result = ( (*(vect1))*(*(vect2)) );

small:
  switch(len)
  {
    case(16):  result += ( (*(++vect1))*(*(++vect2)) );
    case(15):  result += ( (*(++vect1))*(*(++vect2)) );
    case(14):  result += ( (*(++vect1))*(*(++vect2)) );
    case(13):  result += ( (*(++vect1))*(*(++vect2)) );
    case(12):  result += ( (*(++vect1))*(*(++vect2)) );
    case(11):  result += ( (*(++vect1))*(*(++vect2)) );
    case(10):  result += ( (*(++vect1))*(*(++vect2)) );
    case(9):   result += ( (*(++vect1))*(*(++vect2)) );
    case(8):   result += ( (*(++vect1))*(*(++vect2)) );
    case(7):   result += ( (*(++vect1))*(*(++vect2)) );
    case(6):   result += ( (*(++vect1))*(*(++vect2)) );
    case(5):   result += ( (*(++vect1))*(*(++vect2)) );
    case(4):   result += ( (*(++vect1))*(*(++vect2)) );
    case(3):   result += ( (*(++vect1))*(*(++vect2)) );
    case(2):   result += ( (*(++vect1))*(*(++vect2)) );
    case(1):   break;
    case(0):   result = 0.0;
	       break;

    default:
    {
      do
      {
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );
        result += ( (*(++vect1))*(*(++vect2)) );

	len -= 16;
      }
      while(len > 16);

      goto small;
    }
  }

  return(result);
}

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




/***************************************************************************

  This must be as fast as possible, because it is the inner loop for most
  n^3 matrix algebra algorithms.  That is why the code is so ugly.

***************************************************************************/
void wn_add_scaled_vect
(
  register double *to_vect,
  register double *from_vect,
  register double scalar,
  int len
)
{
small:
  switch(len)
  {
    case(16):   *(to_vect++) += (scalar*(*(from_vect++)));
    case(15):   *(to_vect++) += (scalar*(*(from_vect++)));
    case(14):   *(to_vect++) += (scalar*(*(from_vect++)));
    case(13):   *(to_vect++) += (scalar*(*(from_vect++)));
    case(12):   *(to_vect++) += (scalar*(*(from_vect++)));
    case(11):   *(to_vect++) += (scalar*(*(from_vect++)));
    case(10):   *(to_vect++) += (scalar*(*(from_vect++)));
    case(9):    *(to_vect++) += (scalar*(*(from_vect++)));
    case(8):    *(to_vect++) += (scalar*(*(from_vect++)));
    case(7):    *(to_vect++) += (scalar*(*(from_vect++)));
    case(6):    *(to_vect++) += (scalar*(*(from_vect++)));
    case(5):    *(to_vect++) += (scalar*(*(from_vect++)));
    case(4):    *(to_vect++) += (scalar*(*(from_vect++)));
    case(3):    *(to_vect++) += (scalar*(*(from_vect++)));
    case(2):    *(to_vect++) += (scalar*(*(from_vect++)));
    case(1):    *(to_vect++) += (scalar*(*(from_vect++)));
    case(0):    return;
    default:
    {
      do
      {
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));
        *(to_vect++) += (scalar*(*(from_vect++)));

	len -= 16;
      }
      while(len > 16);

      goto small;
    }
  }
}

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




/***************************************************************************

  This must be as fast as possible, because it is the inner loop for 
  conjugate gradient.  That is why the code is so ugly.

***************************************************************************/
void wn_add_vect_and_scaled_vect
(
  register double *out,
  register double *v1,
  register double *v2,
  register double s2,
  int len
)
{
small:
  switch(len)
  {
    case(16):   *(out++) = *(v1++) + (*(v2++))*s2;
    case(15):   *(out++) = *(v1++) + (*(v2++))*s2;
    case(14):   *(out++) = *(v1++) + (*(v2++))*s2;
    case(13):   *(out++) = *(v1++) + (*(v2++))*s2;
    case(12):   *(out++) = *(v1++) + (*(v2++))*s2;
    case(11):   *(out++) = *(v1++) + (*(v2++))*s2;
    case(10):   *(out++) = *(v1++) + (*(v2++))*s2;
    case(9):    *(out++) = *(v1++) + (*(v2++))*s2;
    case(8):    *(out++) = *(v1++) + (*(v2++))*s2;
    case(7):    *(out++) = *(v1++) + (*(v2++))*s2;
    case(6):    *(out++) = *(v1++) + (*(v2++))*s2;
    case(5):    *(out++) = *(v1++) + (*(v2++))*s2;
    case(4):    *(out++) = *(v1++) + (*(v2++))*s2;
    case(3):    *(out++) = *(v1++) + (*(v2++))*s2;
    case(2):    *(out++) = *(v1++) + (*(v2++))*s2;
    case(1):    *(out++) = *(v1++) + (*(v2++))*s2;
    case(0):    return;
    default:
    {
      do
      {
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;
        *(out++) = *(v1++) + (*(v2++))*s2;

	len -= 16;
      }
      while(len > 16);

      goto small;
    }
  }
}


/***************************************************************************

  This must be as fast as possible, because it is the inner loop for 
  conjugate gradient.  That is why the code is so ugly.

***************************************************************************/
void wn_add_scaled_vects
(
  register double *out,
  register double *v1,
  register double s1,
  register double *v2,
  register double s2,
  int len
)
{
small:
  switch(len)
  {
    case(16):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(15):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(14):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(13):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(12):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(11):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(10):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(9):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(8):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(7):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(6):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(5):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(4):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(3):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(2):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(1):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
    case(0):    return;
    default:
    {
      do
      {
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2;

	len -= 16;
      }
      while(len > 16);

      goto small;
    }
  }
}


/***************************************************************************

  This must be as fast as possible, because it is the inner loop for 
  conjugate gradient.  That is why the code is so ugly.

***************************************************************************/
void wn_add_3_scaled_vects
(
  register double *out,
  register double *v1,
  register double s1,
  register double *v2,
  register double s2,
  register double *v3,
  register double s3,
  int len
)
{
small:
  switch(len)
  {
    case(16):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(15):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(14):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(13):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(12):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(11):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(10):   *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(9):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(8):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(7):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(6):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(5):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(4):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(3):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(2):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(1):    *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
    case(0):    return;
    default:
    {
      do
      {
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;
        *(out++) = (*(v1++))*s1 + (*(v2++))*s2 + (*(v3++))*s3;

	len -= 16;
      }
      while(len > 16);

      goto small;
    }
  }
}
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

/*****************************************************************

  This code is the inner loop of various crucial algorithms, so it must
  be fast.  That is why this code is so ugly.

*****************************************************************/
double wn_norm2_vect
(
  register double *vect,
  int len
)
{
  register double result;
  register double num;

  if(len <= 0)
  {
    return(0.0);
  }

  num = *vect;  result = num*num;

small:
  switch(len)
  {
    case(16):  num = *(++vect);  result += num*num;
    case(15):  num = *(++vect);  result += num*num;
    case(14):  num = *(++vect);  result += num*num;
    case(13):  num = *(++vect);  result += num*num;
    case(12):  num = *(++vect);  result += num*num;
    case(11):  num = *(++vect);  result += num*num;
    case(10):  num = *(++vect);  result += num*num;
    case( 9):  num = *(++vect);  result += num*num;
    case( 8):  num = *(++vect);  result += num*num;
    case( 7):  num = *(++vect);  result += num*num;
    case( 6):  num = *(++vect);  result += num*num;
    case( 5):  num = *(++vect);  result += num*num;
    case( 4):  num = *(++vect);  result += num*num;
    case( 3):  num = *(++vect);  result += num*num;
    case( 2):  num = *(++vect);  result += num*num;
    case( 1):   break;
    case( 0):   result = 0.0;
	       break;

    default:
    {
      do
      {
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;
        num = *(++vect);  result += num*num;

	len -= 16;
      }
      while(len > 16);

      goto small;
    }
  }

  return(result);
}


double wn_norm_vect(double *vect,int len)
{
  return pow(wn_norm2_vect(vect,len),0.5);
}


void wn_unit_vect(double vect[],int len)
{
  double scale;
  int i;

  scale = wn_norm_vect(vect,len);
  if(scale == 0.0)
  {
    return;
  }

  scale = 1.0/scale;

  for(i=0;i<len;++i)
  {
    vect[i] *= scale;
  }
}

void wn_zero_vect(register double vect[],register int len)
{
  double test_val;

  /* make sure that memset produces floating-point zeros before using it */
  test_val = 12.345;
  memset(&test_val, 0, sizeof(double));
  if (test_val == 0)
  {
    memset(vect, 0, sizeof(double) * len);
  }
  else
  {
    register int i;

    for(i=0;i<len;++i)
    {
      vect[i] = 0.0;
    }
  }
}

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




/***************************************************************************

  This must be as fast as possible.
  That is why the code is so ugly.

***************************************************************************/
void wn_scale_vect
(
  register double *vect,
  register double scalar,
  int len
)
{
small:
  switch(len)
  {
    case(16):   *(vect++) *= scalar;
    case(15):   *(vect++) *= scalar;
    case(14):   *(vect++) *= scalar;
    case(13):   *(vect++) *= scalar;
    case(12):   *(vect++) *= scalar;
    case(11):   *(vect++) *= scalar;
    case(10):   *(vect++) *= scalar;
    case(9):    *(vect++) *= scalar;
    case(8):    *(vect++) *= scalar;
    case(7):    *(vect++) *= scalar;
    case(6):    *(vect++) *= scalar;
    case(5):    *(vect++) *= scalar;
    case(4):    *(vect++) *= scalar;
    case(3):    *(vect++) *= scalar;
    case(2):    *(vect++) *= scalar;
    case(1):    *(vect++) *= scalar;
    case(0):    return;
    default:
    {
      do
      {
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;
        *(vect++) *= scalar;

	len -= 16;
      }
      while(len > 16);

      goto small;
    }
  }
}

void
wn_multiply_vect_by_vect(double *v1, double *v2, int len)
{
  int i;

  for(i = 0; i < len; i++)
    v1[i] *= v2[i];
}



void wn_print_vect(double vect[],int len)
{
  int i;

  printf("[ ");

  for(i=0;i<len;i++)
  {
    printf("%lg ",vect[i]);
  }

  printf("]\n");
}


void wn_mult_polys(double out[],double in1[],int len1,double in2[],int len2)
{
  int i1,i2;

  wn_zero_vect(out,len1+len2-1);

  for(i1=0;i1<len1;++i1)
  for(i2=0;i2<len2;++i2)
  {
    out[i1+i2] += in1[i1]*in2[i2];
  }
}


