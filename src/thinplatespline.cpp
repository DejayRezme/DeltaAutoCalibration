///*
// *  Thin Plate Spline demo/example in C++
// *
// *  - a simple TPS editor, using the Boost uBlas library for large
// *    matrix operations and OpenGL + GLUT for 2D function visualization
// *    (curved plane) and user interface
// *
// *  Copyright (C) 2003,2005 by Jarno Elonen
// *
// *  TPSDemo is Free Software / Open Source with a very permissive
// *  license:
// *
// *  Permission to use, copy, modify, distribute and sell this software
// *  and its documentation for any purpose is hereby granted without fee,
// *  provided that the above copyright notice appear in all copies and
// *  that both that copyright notice and this permission notice appear
// *  in supporting documentation.  The authors make no representations
// *  about the suitability of this software for any purpose.
// *  It is provided "as is" without express or implied warranty.
// *
// *  TODO:
// *    - implement TPS approximation 3 as suggested in paper
// *      Gianluca Donato and Serge Belongie, 2002: "Approximation
// *      Methods for Thin Plate Spline Mappings and Principal Warps"
// */
//
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
//
////#include "linalg3d.h"
////#include "ludecomposition.h"
//
//#include <vector>
//#include <cmath>
//#include <cstdio>
//#include <cstring>
//
//using namespace boost::numeric::ublas;
//
//// ========= BEGIN INTERESTING STUFF  =========
//#define EPSILON 0.00001f
//#define PI 3.1415926
//#define Deg2Rad(Ang) ((float)( Ang * PI / 180.0 ))
//#define Rad2Deg(Ang) ((float)( Ang * 180.0 / PI ))
//
//class Vec
//{
//public:
//
//  // Position
//  float x, y, z;
//
//  // Default constructor
//  Vec()
//  : x( 0 ), y( 0 ), z( 0 ) {}
//
//  // Element constructor
//  Vec( float x, float y, float z )
//  : x( x ), y( y ), z( z ) {}
//
//  // Copy constructor
//  Vec( const Vec& a )
//  : x( a.x ), y( a.y ), z( a.z ) {}
//
//  // Norm (len^2)
//  inline float norm() const { return x*x + y*y + z*z; }
//
//  // Length of the vector
//  inline float len() const { return (float)sqrt(norm()); }
//
//  Vec &operator += ( const Vec &src ) { x += src.x; y += src.y; z += src.z; return *this; }
//  Vec operator + ( const Vec &src ) const { Vec tmp( *this ); return ( tmp += src ); }
//  Vec &operator -= ( const Vec &src ) { x -= src.x; y -= src.y; z -= src.z; return *this; }
//  Vec operator - ( const Vec &src ) const { Vec tmp( *this ); return ( tmp -= src ); }
//
//  Vec operator - () const { return Vec(-x,-y,-z); }
//
//  Vec &operator *= ( const float src ) { x *= src; y *= src; z *= src;  return *this; }
//  Vec operator * ( const float src ) const { Vec tmp( *this ); return ( tmp *= src ); }
//  Vec &operator /= ( const float src ) { x /= src; y /= src; z /= src; return *this; }
//  Vec operator / ( const float src ) const { Vec tmp( *this ); return ( tmp /= src ); }
//
//  bool operator == ( const Vec& b) const { return ((*this)-b).norm() < EPSILON; }
//  //bool operator == ( const Vec& b) const { return x==b.x && y==b.y && z==b.z; }
//};
//
//
//// Left hand float multplication
//inline Vec operator * ( const float src, const Vec& v ) { Vec tmp( v ); return ( tmp *= src ); }
//
//// Dot product
//inline float dot( const Vec& a, const Vec& b )
//{ return a.x*b.x + a.y*b.y + a.z*b.z; }
//
//// Cross product
//inline Vec cross( const Vec &a, const Vec &b )
//{ return Vec( a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x ); }
//
//
//
//
//
//
//// Solve a linear equation system a*x=b using inplace LU decomposition.
////
//// Stores x in 'b' and overwrites 'a' (with a pivotted LUD).
////
//// Matrix 'b' may have any (>0) number of columns but
//// must contain as many rows as 'a'.
////
//// Possible return values:
////  0=success
////  1=singular matrix
////  2=a.rows != b.rows
//template <typename T> int LU_Solve(
//  boost::numeric::ublas::matrix<T>& a,
//  boost::numeric::ublas::matrix<T>& b )
//{
//  // This routine is originally based on the public domain draft for JAMA,
//  // Java matrix package available at http://math.nist.gov/javanumerics/jama/
//
//  typedef boost::numeric::ublas::matrix<T> Matrix;
//  typedef boost::numeric::ublas::matrix_row<Matrix> Matrix_Row;
//  typedef boost::numeric::ublas::matrix_column<Matrix> Matrix_Col;
//
//  if (a.size1() != b.size1())
//    return 2;
//
//  int m = a.size1(), n = a.size2();
//  int pivsign = 0;
//  int* piv = (int*)alloca( sizeof(int) * m);
//
//  // PART 1: DECOMPOSITION
//  //
//  // For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
//  // unit lower triangular matrix L, an n-by-n upper triangular matrix U,
//  // and a permutation vector piv of length m so that A(piv,:) = L*U.
//  // If m < n, then L is m-by-m and U is m-by-n.
//  {
//    // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
//    for (int i = 0; i < m; ++i)
//      piv[i] = i;
//    pivsign = 1;
//
//    // Outer loop.
//    for (int j=0; j<n; ++j)
//    {
//      // Make a copy of the j-th column to localize references.
//      Matrix_Col LUcolj(a,j);
//
//      // Apply previous transformations.
//      for (int i = 0; i < m; ++i)
//      {
//          Matrix_Row LUrowi(a,i);
//
//          // This dot product is very expensive.
//          // Optimize for SSE2?
//          int kmax = (i<=j)?i:j;
//          typename Matrix_Row::const_iterator ri_ite( LUrowi.begin());
//          typename Matrix_Col::const_iterator cj_ite( LUcolj.begin());
//          typename Matrix::value_type sum = 0.0;
//          while( kmax-- > 0 )
//            sum += (*(ri_ite++)) * (*(cj_ite++));
//          LUrowi[j] = LUcolj[i] -= sum;
//      }
//
//      // Find pivot and exchange if necessary.
//      //
//      // Slightly optimized version of:
//      //  for (int i = j+1; i < m; ++i)
//      //    if ( fabs(LUcolj[i]) > fabs(LUcolj[p]) )
//      //      p = i;
//      int p = j;
//      typename Matrix::value_type coljp_abs = fabs(LUcolj[p]);
//      for ( typename Matrix_Col::const_iterator
//              beg = LUcolj.begin(),
//              ite = beg + j+1,
//              end = LUcolj.end();
//            ite < end;
//            ++ite )
//      {
//        if (fabs(*ite) > coljp_abs)
//        {
//          p = ite-beg;
//          coljp_abs = fabs(LUcolj[p]);
//        }
//      }
//
//      if (p != j)
//      {
//          Matrix_Row raj(a,j);
//          Matrix_Row(a,p).swap(raj);
//
//          int tmp = piv[p];
//          piv[p] = piv[j];
//          piv[j] = tmp;
//          pivsign = -pivsign;
//      }
//
//      // Compute multipliers.
//      if (j < m && a(j,j) != 0.0)
//          for (int i = j+1; i < m; ++i)
//            LUcolj[i] /= LUcolj[j];
//    }
//  }
//
//  // PART 2: SOLVE
//
//  // Check singluarity
//  for (int j = 0; j < n; ++j)
//    if (a(j,j) == 0)
//      return 1;
//
//  // Reorder b according to pivotting
//  for (int i=0; i<m; ++i)
//  {
//    if ( piv[i] != i )
//    {
//      Matrix_Row b_ri( b, i );
//      Matrix_Row( b, piv[i] ).swap( b_ri );
//      for ( int j=i; j<m; ++j )
//        if ( piv[j] == i )
//        {
//          piv[j] = piv[i];
//          break;
//        }
//    }
//  }
//
//  // Solve L*Y = B(piv,:)
//  for (int k=0; k<n; ++k)
//  {
//    const Matrix_Row& b_rk = Matrix_Row( b, k );
//    for (int i = k+1; i < n; ++i)
//    {
//      const typename Matrix_Row::value_type aik = a(i,k);
//      Matrix_Row( b, i ) -= b_rk * aik;
//    }
//  }
//
//  // Solve U*X = Y;
//  for (int k=n-1; k>=0; --k)
//  {
//    Matrix_Row(b,k) *= 1.0/a(k,k);
//
//    const Matrix_Row& b_rk = Matrix_Row(b, k );
//    for (int i=0; i<k; ++i)
//    {
//      const typename Matrix_Row::value_type aik = a(i,k);
//      Matrix_Row(b,i) -= b_rk * aik;
//    }
//  }
//
//  return 0;
//}
//
//
//
//
//
//
//
//
//
//
//
//#define GRID_W 10
//#define GRID_H 10
//static float grid[GRID_W][GRID_H];
//
//
//std::vector< Vec > control_points;
//
//int selected_cp = -1;
//
//double regularization = 0.0;
//double bending_energy = 0.0;
//
//static double tps_base_func(double r)
//{
//  if ( r == 0.0 )
//    return 0.0;
//  else
//    return r*r * log(r);
//}
//
//
///*
// *  Calculate Thin Plate Spline (TPS) weights from
// *  control points and build a new height grid by
// *  interpolating with them.
// */
//static void calc_tps()
//{
//  // You We need at least 3 points to define a plane
//  if ( control_points.size() < 3 )
//    return;
//
//  unsigned p = control_points.size();
//
//  // Allocate the matrix and vector
//  matrix<double> mtx_l(p+3, p+3);
//  matrix<double> mtx_v(p+3, 1);
//  matrix<double> mtx_orig_k(p, p);
//
//  // Fill K (p x p, upper left of L) and calculate
//  // mean edge length from control points
//  //
//  // K is symmetrical so we really have to
//  // calculate only about half of the coefficients.
//  double a = 0.0;
//  for ( unsigned i=0; i<p; ++i )
//  {
//    for ( unsigned j=i+1; j<p; ++j )
//    {
//      Vec pt_i = control_points[i];
//      Vec pt_j = control_points[j];
//      pt_i.y = pt_j.y = 0;
//      double elen = (pt_i - pt_j).len();
//      mtx_l(i,j) = mtx_l(j,i) =
//        mtx_orig_k(i,j) = mtx_orig_k(j,i) =
//          tps_base_func(elen);
//      a += elen * 2; // same for upper & lower tri
//    }
//  }
//  a /= (double)(p*p);
//
//  // Fill the rest of L
//  for ( unsigned i=0; i<p; ++i )
//  {
//    // diagonal: reqularization parameters (lambda * a^2)
//    mtx_l(i,i) = mtx_orig_k(i,i) =
//      regularization * (a*a);
//
//    // P (p x 3, upper right)
//    mtx_l(i, p+0) = 1.0;
//    mtx_l(i, p+1) = control_points[i].x;
//    mtx_l(i, p+2) = control_points[i].z;
//
//    // P transposed (3 x p, bottom left)
//    mtx_l(p+0, i) = 1.0;
//    mtx_l(p+1, i) = control_points[i].x;
//    mtx_l(p+2, i) = control_points[i].z;
//  }
//  // O (3 x 3, lower right)
//  for ( unsigned i=p; i<p+3; ++i )
//    for ( unsigned j=p; j<p+3; ++j )
//      mtx_l(i,j) = 0.0;
//
//
//  // Fill the right hand vector V
//  for ( unsigned i=0; i<p; ++i )
//    mtx_v(i,0) = control_points[i].y;
//  mtx_v(p+0, 0) = mtx_v(p+1, 0) = mtx_v(p+2, 0) = 0.0;
//
//  // Solve the linear system "inplace"
//  if (0 != LU_Solve(mtx_l, mtx_v))
//  {
//    puts( "Singular matrix! Aborting." );
//    exit(1);
//  }
//
//  // Interpolate grid heights
//  for ( int x=-GRID_W/2; x<GRID_W/2; ++x )
//  {
//    for ( int z=-GRID_H/2; z<GRID_H/2; ++z )
//    {
//      double h = mtx_v(p+0, 0) + mtx_v(p+1, 0)*x + mtx_v(p+2, 0)*z;
//      Vec pt_i, pt_cur(x,0,z);
//      for ( unsigned i=0; i<p; ++i )
//      {
//        pt_i = control_points[i];
//        pt_i.y = 0;
//        h += mtx_v(i,0) * tps_base_func( ( pt_i - pt_cur ).len());
//      }
//      grid[x+GRID_W/2][z+GRID_H/2] = h;
//      printf("%2f ", h);
//    }
//    printf("\n");
//  }
//
//  // Calc bending energy
//  matrix<double> w( p, 1 );
//  for (unsigned int i=0; i<p; ++i )
//    w(i,0) = mtx_v(i,0);
////  matrix<double> be = prod( prod<matrix<double> >( trans(w), mtx_orig_k ), w );
////  bending_energy = be(0,0);
//}
//
//
//static void clear_grid()
//{
//  for (int x=0; x<GRID_W; ++x)
//    for (int z=0; z<GRID_H; ++z)
//      grid[x][z] = 0;
//}
//
//
//// Startup
//int main( int argc, char *argv[] )
//{
//	printf("hello\n");
//  clear_grid();
//
//  control_points.push_back(Vec(0, 1, 0));
//  control_points.push_back(Vec(1, 0, 0));
//  control_points.push_back(Vec(0, 0, 1));
//  control_points.push_back(Vec(-1, 0, 0));
//  control_points.push_back(Vec(0, 0, -1));
//  control_points.push_back(Vec(2, 0, 0));
//  control_points.push_back(Vec(0, 0, 2));
//  control_points.push_back(Vec(-2, 0, 0));
//  control_points.push_back(Vec(0, 0, -2));
//
//  calc_tps();
//
//  return 0;
//}
