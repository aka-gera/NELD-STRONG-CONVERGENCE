#ifndef MATRIX_HPP
#define MATRIX_HPP

//------------------------------------------ 
//          Simple matrix class 
//-------------------------------------------

//--- required libraries ---
#include <limits>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std; 

//------- structure --------
class Matrix
{
public:

  //--- fields ---
  int row, col;
  double * coord;

  //--- constructors, destructors ---
  Matrix();
  ~Matrix();
  Matrix(int,int);
  double & operator() (int,int);
  void set_size(int,int);
  void zeros();
  void identity();
  Matrix inverse();  //Only for 2 by 2
  Matrix subslice(int, int, int, int);
  Matrix transpose();
  Matrix crosses();
  double column_norm(int);
  double norm();
  double contraction_product(Matrix&);
  Matrix exp_diag(); // Called on a vector, returns a diagonal matrix.
  
  // overloaded operations
  Matrix & operator= (const Matrix &);
  Matrix operator+ (const Matrix &);
  Matrix operator- (const Matrix &);
  Matrix operator* (const Matrix &);
  Matrix operator* (const double); 
  Matrix& operator+= (const Matrix &);
  Matrix& operator-= (const Matrix &);
  
  friend ostream& operator<<(ostream&, const Matrix&);
  friend istream& operator>>(istream&, Matrix&);
};

#endif
