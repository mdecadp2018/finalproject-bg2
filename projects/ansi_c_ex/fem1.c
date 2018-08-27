// FEM1D_BVP_LINEAR is a C program which applies the finite element method, with piecewise linear elements, to a two point boundary value problem in one spatial dimension.

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

double *fem1d_bvp_linear ( int n, double a ( double x ), double c ( double x ), 
  double f ( double x ), double x[] );
int *i4vec_zero_new ( int n );
double r8_abs ( double x );
double *r8mat_solve2 ( int n, double a[], double b[], int *ierror );
double *r8mat_zero_new ( int m, int n );
double *r8vec_even ( int n, double alo, double ahi );
double *r8vec_zero_new ( int n );
void timestamp ( );

/******************************************************************************/

double *fem1d_bvp_linear ( int n, double a ( double x ), double c ( double x ), 
  double f ( double x ), double x[] )

/******************************************************************************/
/*
  Purpose:

    FEM1D_BVP_LINEAR solves a two point boundary value problem.

  Discussion:

    The program uses the finite element method, with piecewise linear basis
    functions to solve a boundary value problem in one dimension.

    The problem is defined on the region 0 <= x <= 1.

    The following differential equation is imposed between 0 and 1:

      - d/dx a(x) du/dx + c(x) * u(x) = f(x)

    where a(x), c(x), and f(x) are given functions.

    At the boundaries, the following conditions are applied:

      u(0.0) = 0.0
      u(1.0) = 0.0

    A set of N equally spaced nodes is defined on this
    interval, with 0 = X(1) < X(2) < ... < X(N) = 1.0.

    At each node I, we associate a piecewise linear basis function V(I,X),
    which is 0 at all nodes except node I.  This implies that V(I,X) is
    everywhere 0 except that

    for X(I-1) <= X <= X(I):

      V(I,X) = ( X - X(I-1) ) / ( X(I) - X(I-1) ) 

    for X(I) <= X <= X(I+1):

      V(I,X) = ( X(I+1) - X ) / ( X(I+1) - X(I) )

    We now assume that the solution U(X) can be written as a linear
    sum of these basis functions:

      U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X)

    where U(X) on the left is the function of X, but on the right,
    is meant to indicate the coefficients of the basis functions.

    To determine the coefficient U(J), we multiply the original
    differential equation by the basis function V(J,X), and use
    integration by parts, to arrive at the I-th finite element equation:

        Integral A(X) * U'(X) * V'(I,X) + C(X) * U(X) * V(I,X) dx 
      = Integral F(X) * V(I,X) dx

    We note that the functions U(X) and U'(X) can be replaced by
    the finite element form involving the linear sum of basis functions,
    but we also note that the resulting integrand will only be nonzero
    for terms where J = I - 1, I, or I + 1.

    By writing this equation for basis functions I = 2 through N - 1,
    and using the boundary conditions, we have N linear equations
    for the N unknown coefficients U(1) through U(N), which can
    be easily solved.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of nodes.

    Input, double A ( double X ), evaluates a(x);

    Input, double C ( double X ), evaluates c(x);

    Input, double F ( double X ), evaluates f(x);

    Input, double X[N], the mesh points.

    Output, double FEM1D_BVP_LINEAR[N], the finite element coefficients, 
    which are also the value of the computed solution at the mesh points.
*/
{
# define QUAD_NUM 2

  double abscissa[QUAD_NUM] = {
    -0.577350269189625764509148780502,
    +0.577350269189625764509148780502 };
  double al;
  double am;
  double ar;
  double *amat;
  double axq;
  double *b;
  double bm;
  double cxq;
  double fxq;
  double h;
  int i;
  int ierror;
  int q;
  int quad_num = QUAD_NUM;
  double *u;
  double weight[QUAD_NUM] = { 1.0, 1.0 };
  double wq;
  double vl;
  double vlp;
  double vm;
  double vmp;
  double vr;
  double vrp;
  double xl;
  double xm;
  double xq;
  double xr;
/*
  Zero out the matrix and right hand side.
*/
  amat = r8mat_zero_new ( n, n );
  b = r8vec_zero_new ( n );
/*
  Equation 1 is the left boundary condition, U(0.0) = 0.0;
*/
  amat[0+0*n] = 1.0;
  b[0] = 0.0;
/*
  Equation I involves the basis function at node I.
  This basis function is nonzero from X(I-1) to X(I+1).
  Equation I looks like this:

    Integral A(X) U'(X) V'(I,X) 
           + C(X) * U(X) V(I,X) dx 
  = Integral F(X) V(I,X) dx

  Then, we realize that U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X), 
  (U(X) means the function; U(J) is the coefficient of V(J,X) ).

  The only V functions that are nonzero when V(I,X) is nonzero are
  V(I-1,X) and V(I+1,X). 

  Let's use the shorthand 

    VL(X) = V(I-1,X)
    VM(X) = V(I,X)
    VR(X) = V(I+1,X)

  So our equation becomes

    Integral A(X) [ VL'(X) U(I-1) + VM'(X) U(I) + VR'(X) U(I+1) ] * VM'(X)
           + C(X) [ VL(X)  U(I-1) + VM(X)  U(I) + VR(X)  U(I+1) ] * VM(X) dx
  = Integral F(X) VM(X) dx.

  

  This is actually a set of N-2 linear equations for the N coefficients U.

  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1), 
  and so on.
*/
  for ( i = 1; i < n - 1; i++ )
  {
/*
  Get the left, right and middle coordinates.
*/
    xl = x[i-1];
    xm = x[i];
    xr = x[i+1];
/*
  Make temporary variables for A(I,I-1), A(I,I), A(I,I+1) and B(I).
*/
    al = 0.0;
    am = 0.0;
    ar = 0.0;
    bm = 0.0;
/*
  We approximate the integrals by using a weighted sum of
  the integrand values at quadrature points.
*/
    for ( q = 0; q < quad_num; q++ )
    {
/*
  Integrate over the LEFT interval, between XL and XM, where:

  VL(X) = ( XM - X       ) / ( XM - XL )
  VM(X) = (      X  - XL ) / ( XM - XL )
  VR(X) = 0

  VL'(X) =             - 1 / ( XM - XL )
  VM'(X) =             + 1 / ( XM - XL ) 
  VR'(X) = 0
*/
      xq = ( ( 1.0 - abscissa[q] ) * xl 
           + ( 1.0 + abscissa[q] ) * xm ) 
           /   2.0;

      wq = weight[q] * ( xm - xl ) / 2.0;

      vl =  ( xm - xq ) / ( xm - xl );
      vlp =      - 1.0  / ( xm - xl );

      vm =  ( xq - xl ) / ( xm - xl );
      vmp =      + 1.0  / ( xm - xl );

      vr =  0.0;
      vrp = 0.0;

      axq = a ( xq );
      cxq = c ( xq );
      fxq = f ( xq );

      al = al + wq * ( axq * vlp * vmp + cxq * vl * vm );
      am = am + wq * ( axq * vmp * vmp + cxq * vm * vm );
      ar = ar + wq * ( axq * vrp * vmp + cxq * vr * vm );
      bm = bm + wq * ( fxq * vm );
/*
  Integrate over the RIGHT interval, between XM and XR, where:

  VL(X) = 0
  VM(X) = ( XR - X       ) / ( XR - XM )
  VR(X) = (      X  - XM ) / ( XR - XM )

  VL'(X) = 0
  VM'(X) =             - 1 / ( XR - XM )
  VR'(X) =             + 1 / ( XR - XM ) 
*/
      xq = ( ( 1.0 - abscissa[q] ) * xm 
           + ( 1.0 + abscissa[q] ) * xr ) 
           /   2.0;

      wq = weight[q] * ( xr - xm ) / 2.0;

      vl = 0.0;
      vlp = 0.0;

      vm = ( xr - xq ) / ( xr - xm );
      vmp =     - 1.0  / ( xr - xm );

      vr = ( xq - xm ) / ( xr - xm );
      vrp =      1.0   / ( xr - xm );

      axq = a ( xq );
      cxq = c ( xq );
      fxq = f ( xq );

      al = al + wq * ( axq * vlp * vmp + cxq * vl * vm );
      am = am + wq * ( axq * vmp * vmp + cxq * vm * vm );
      ar = ar + wq * ( axq * vrp * vmp + cxq * vr * vm );
      bm = bm + wq * ( fxq * vm );
    }
    amat[i+(i-1)*n] = al;
    amat[i+ i   *n] = am;
    amat[i+(i+1)*n] = ar;

    b[i] = bm;
  }
/*
  Equation N is the right boundary condition, U(1.0) = 0.0;
*/
  amat[n-1+(n-1)*n] = 1.0;
  b[n-1] = 0.0;
/*
  Solve the linear system.
*/
  u = r8mat_solve2 ( n, amat, b, &ierror );

  free ( amat );
  free ( b );

  return u;
# undef QUAD_NUM
}
/******************************************************************************/

int *i4vec_zero_new ( int n )

/******************************************************************************/
/*
  Purpose:

    I4VEC_ZERO_NEW creates and zeroes an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
*/
{
  int *a;
  int i;

  a = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double *r8mat_solve2 ( int n, double a[], double b[], int *ierror )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SOLVE2 computes the solution of an N by N linear system.

  Discussion: 							    

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
    in column-major order.

    The linear system may be represented as

      A*X = B

    If the linear system is singular, but consistent, then the routine will
    still produce a solution.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of equations.

    Input/output, double A[N*N].
    On input, A is the coefficient matrix to be inverted.
    On output, A has been overwritten.

    Input/output, double B[N].
    On input, B is the right hand side of the system.
    On output, B has been overwritten.

    Output, double R8MAT_SOLVE2[N], the solution of the linear system.

    Output, int *IERROR.
    0, no error detected.
    1, consistent singularity.
    2, inconsistent singularity.
*/
{
  double amax;
  int i;
  int imax;
  int j;
  int k;
  int *piv;
  double *x;

  *ierror = 0;

  piv = i4vec_zero_new ( n );
  x = r8vec_zero_new ( n );
/*
  Process the matrix.
*/
  for ( k = 1; k <= n; k++ )
  {
/*
  In column K:
    Seek the row IMAX with the properties that:
      IMAX has not already been used as a pivot;
      A(IMAX,K) is larger in magnitude than any other candidate.
*/
    amax = 0.0;
    imax = 0;
    for ( i = 1; i <= n; i++ )
    {
      if ( piv[i-1] == 0 )
      {
        if ( amax < r8_abs ( a[i-1+(k-1)*n] ) )
        {
          imax = i;
          amax = r8_abs ( a[i-1+(k-1)*n] );
        }
      }
    }
/*
  If you found a pivot row IMAX, then,
    eliminate the K-th entry in all rows that have not been used for pivoting.
*/
    if ( imax != 0 )
    {
      piv[imax-1] = k;
      for ( j = k+1; j <= n; j++ )
      {
        a[imax-1+(j-1)*n] = a[imax-1+(j-1)*n] / a[imax-1+(k-1)*n];
      }
      b[imax-1] = b[imax-1] / a[imax-1+(k-1)*n];
      a[imax-1+(k-1)*n] = 1.0;

      for ( i = 1; i <= n; i++ )
      {
        if ( piv[i-1] == 0 )
        {
          for ( j = k+1; j <= n; j++ )
          {
            a[i-1+(j-1)*n] = a[i-1+(j-1)*n] - a[i-1+(k-1)*n] * a[imax-1+(j-1)*n];
          }
          b[i-1] = b[i-1] - a[i-1+(k-1)*n] * b[imax-1];
          a[i-1+(k-1)*n] = 0.0;
        }
      }
    }
  }
/*
  Now, every row with nonzero PIV begins with a 1, and
  all other rows are all zero.  Begin solution.
*/
  for ( j = n; 1 <= j; j-- )
  {
    imax = 0;
    for ( k = 1; k <= n; k++ )
    {
      if ( piv[k-1] == j )
      {
        imax = k;
      }
    }

    if ( imax == 0 )
    {
      x[j-1] = 0.0;

      if ( b[j-1] == 0.0 )
      {
        *ierror = 1;
        printf ( "\n" );
        printf ( "R8MAT_SOLVE2 - Warning:\n" );
        printf ( "  Consistent singularity, equation = %d\n", j );
      }
      else
      {
        *ierror = 2;
        printf ( "\n" );
        printf ( "R8MAT_SOLVE2 - Warning:\n" );
        printf ( "  Inconsistent singularity, equation = %d\n", j );
      }
    }
    else
    {
      x[j-1] = b[imax-1];

      for ( i = 1; i <= n; i++ )
      {
        if ( i != imax )
        {
          b[i-1] = b[i-1] - a[i-1+(j-1)*n] * x[j-1];
        }
      }
    }
  }

  free ( piv );

  return x;
}
/******************************************************************************/

double *r8mat_zero_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ZERO_NEW returns a new zeroed R8MAT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, double R8MAT_ZERO[M*N], the new zeroed matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
/******************************************************************************/

double *r8vec_even ( int n, double alo, double ahi )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EVEN returns N real values, evenly spaced between ALO and AHI.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 February 2004

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values.

    Input, double ALO, AHI, the low and high values.

    Output, double R8VEC_EVEN[N], N evenly spaced values.
*/
{
  double *a;
  int i;

  a = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 1; i <= n; i++ )
    {
      a[i-1] = ( ( double ) ( n - i     ) * alo 
               + ( double ) (     i - 1 ) * ahi ) 
               / ( double ) ( n     - 1 );
    }
  }

  return a;
}
/******************************************************************************/

double *r8vec_zero_new ( int n )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ZERO_NEW creates and zeroes an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
*/
{
  double *a;
  int i;

  a = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return a;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

int main ( void );
void fem1d_bvp_linear_test01 ( void );
void fem1d_bvp_linear_test02 ( void );
void fem1d_bvp_linear_test03 ( void );
void fem1d_bvp_linear_test04 ( void );
void fem1d_bvp_linear_test05 ( void );
double a1 ( double x );
double a2 ( double x );
double a3 ( double x );
double c1 ( double x );
double c2 ( double x );
double c3 ( double x );
double f1 ( double x );
double f2 ( double x );
double f3 ( double x );
double f4 ( double x );
double f5 ( double x );
double exact1 ( double x );
double exact2 ( double x );
double exact3 ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    FEM1D_BVP_LINEAR_PRB tests the routines in FEM1D_BVP_LINEAR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FEM1D_BVP_LINEAR_PRB\n" );
  printf ( "  C version\n" );

  fem1d_bvp_linear_test01 ( );
  fem1d_bvp_linear_test02 ( );
  fem1d_bvp_linear_test03 ( );
  fem1d_bvp_linear_test04 ( );
  fem1d_bvp_linear_test05 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM1D_BVP_LINEAR_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void fem1d_bvp_linear_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    FEM1D_BVP_LINEAR_TEST01 carries out test case #1.

  Discussion:

    Use A1, C1, F1, EXACT1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n = 11;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "FEM1D_BVP_LINEAR_TEST01\n" );
  printf ( "  A1(X)  = 1.0\n" );
  printf ( "  C1(X)  = 0.0\n" );
  printf ( "  F1(X)  = X * ( X + 3 ) * exp ( X )\n" );
  printf ( "  U1(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a1, c1, f1, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact1 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, r8_abs ( u[i] - uexact ) );
  }

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

void fem1d_bvp_linear_test02 ( void )

/******************************************************************************/
/*
  Purpose:

    FEM1D_BVP_LINEAR_TEST02 carries out test case #2.

  Discussion:

    Use A1, C2, F2, EXACT1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n = 11;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "FEM1D_BVP_LINEAR_TEST02\n" );
  printf ( "  A1(X)  = 1.0\n" );
  printf ( "  C2(X)  = 2.0\n" );
  printf ( "  F2(X)  = X * ( 5 - X ) * exp ( X )\n" );
  printf ( "  U1(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a1, c2, f2, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact1 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, r8_abs ( u[i] - uexact ) );
  }

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

void fem1d_bvp_linear_test03 ( void )

/******************************************************************************/
/*
  Purpose:

    FEM1D_BVP_LINEAR_TEST03 carries out test case #3.

  Discussion:

    Use A1, C3, F3, EXACT1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n = 11;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "FEM1D_BVP_LINEAR_TEST03\n" );
  printf ( "  A1(X)  = 1.0\n" );
  printf ( "  C3(X)  = 2.0 * X\n" );
  printf ( "  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )\n" );
  printf ( "  U1(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a1, c3, f3, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact1 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, r8_abs ( u[i] - uexact ) );
  }

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

void fem1d_bvp_linear_test04 ( void )

/******************************************************************************/
/*
  Purpose:

    FEM1D_BVP_LINEAR_TEST04 carries out test case #4.

  Discussion:

    Use A2, C1, F4, EXACT1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n = 11;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "FEM1D_BVP_LINEAR_TEST04\n" );
  printf ( "  A2(X)  = 1.0 + X * X\n" );
  printf ( "  C1(X)  = 0.0\n" );
  printf ( "  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )\n" );
  printf ( "  U1(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a2, c1, f4, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact1 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, r8_abs ( u[i] - uexact ) );
  }

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

void fem1d_bvp_linear_test05 ( )

/******************************************************************************/
/*
  Purpose:

    FEM1D_BVP_LINEAR_TEST05 carries out test case #5.

  Discussion:

    Use A3, C1, F5, EXACT1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n = 11;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "FEM1D_BVP_LINEAR_TEST05\n" );
  printf ( "  A3(X)  = 1.0 + X * X for X <= 1/3\n" );
  printf ( "         = 7/9 + X     for      1/3 < X\n" );
  printf ( "  C1(X)  = 0.0\n" );
  printf ( "  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )\n" );
  printf ( "                       for X <= 1/3\n" );
  printf ( "         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) .* exp ( X )\n" );
  printf ( "                       for      1/3 <= X\n" );
  printf ( "  U1(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a3, c1, f5, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact1 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, r8_abs ( u[i] - uexact ) );
  }

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

double a1 ( double x )

/******************************************************************************/
/*
  Purpose:

    A1 evaluates A function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A1, the value of A(X).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double a2 ( double x )

/******************************************************************************/
/*
  Purpose:

    A2 evaluates A function #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A2, the value of A(X).
*/
{
  double value;

  value = 1.0 + x * x;

  return value;
}
/******************************************************************************/

double a3 ( double x )

/******************************************************************************/
/*
  Purpose:

    A3 evaluates A function #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A3, the value of A(X).
*/
{
  double value;

  if ( x <= 1.0 / 3.0 )
  {
    value = 1.0 + x * x;
  }
  else
  {
    value = x + 7.0 / 9.0;
  }

  return value;
}
/******************************************************************************/

double c1 ( double x )

/******************************************************************************/
/*
  Purpose:

    C1 evaluates C function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C1, the value of C(X).
*/
{
  double value;

  value = 0.0;

  return value;
}
/******************************************************************************/

double c2 ( double x )

/******************************************************************************/
/*
  Purpose:

    C2 evaluates C function #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C2, the value of C(X).
*/
{
  double value;

  value = 2.0;

  return value;
}
/******************************************************************************/

double c3 ( double x )

/******************************************************************************/
/*
  Purpose:

    C3 evaluates C function #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C3, the value of C(X).
*/
{
  double value;

  value = 2.0 * x;

  return value;
}
/******************************************************************************/

double f1 ( double x )

/******************************************************************************/
/*
  Purpose:

    F1 evaluates right hand side function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F1, the value of F(X).
*/
{
  double value;

  value = x * ( x + 3.0 ) * exp ( x );

  return value;
}
/******************************************************************************/

double f2 ( double x )

/******************************************************************************/
/*
  Purpose:

    F2 evaluates right hand side function #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F2, the value of F(X).
*/
{
  double value;

  value = x * ( 5.0 - x ) * exp ( x );

  return value;
}
/******************************************************************************/

double f3 ( double x )

/******************************************************************************/
/*
  Purpose:

    F3 evaluates right hand side function #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F3, the value of F(X).
*/
{
  double value;

  value = - x * ( 2.0 * x * x - 3.0 * x - 3.0 ) * exp ( x );

  return value;
}
/******************************************************************************/

double f4 ( double x )

/******************************************************************************/
/*
  Purpose:

    F4 evaluates right hand side function #4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F4, the value of F(X).
*/
{
  double value;

  value = ( x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x ) * exp ( x );

  return value;
}
/******************************************************************************/

double f5 ( double x )

/******************************************************************************/
/*
  Purpose:

    F5 evaluates right hand side function #5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F5, the value of F(X).
*/
{
  double value;

  if ( x <= 1.0 / 3.0 )
  {
    value = ( x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x ) * exp ( x );
  }
  else
  {
    value = ( - 1.0 + ( 10.0 / 3.0 ) * x 
      + ( 43.0 / 9.0 ) * x * x + x * x * x ) * exp ( x );
  }

  return value;
}
/******************************************************************************/

double exact1 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT1 evaluates exact solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT1, the value of U(X).
*/
{
  double value;

  value = x * ( 1.0 - x ) * exp ( x );

  return value;
}
/******************************************************************************/

double exact2 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT2 returns exact solution #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT2, the value of U(X).
*/
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value =  x * ( 1.0 - x ) * exp ( x );
  }
  else
  {
    value = x * ( 1.0 - x )  * exp ( 2.0 / 3.0 );
  }
  return value;
}
/******************************************************************************/

double exact3 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT3 returns exact solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT3, the value of U(X).
*/
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value = x * ( 1.0 - x ) * exp ( x );
  }
  else
  {
    value = x * ( 1.0 - x );
  }
  return value;
}