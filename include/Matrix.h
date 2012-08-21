/********************************************************************
 * Matrix.h - header file of the Matrix library which defines 
 *           a new class "Matrix" and the associated member functions
 *
 * Note: 
 *   This is a simple C++ library for matrix processing. 
 *   The purpose is not for high performance, but to show how 
 *   the algorithm works through programming.
 *
 * Copyright (C) hqi@utk.edu, Jan. 2002
 *
 * Other contributors:
 *   - Ryan Kerekes: determinant()
 *   - Xiaoling Wang: jacobi() based on "Numerical Recipe"
 *   - Xiaoling Wang: eigsrt() based on "Numerical Recipe"
 *   - Xiaoling Wang: feature extraction related functions
 *
 * Modifications:
 *   - 01/12/08: separate jacobi() and eigsrt(), it does have a reason!
 *   - 01/11/08: compress the Matrix class member function,
 *               make most C-style functions
 *   - 05/03/05: combined jacobi() and eigsrt() 
 *               added distance() and insertsort()
 *   - 04/07/05: fix jacobi() that does not alter matrix (Tom Karnowski)
 *   - 04/07/05: fix readPPM and writePPM using C style file I/O
 *   - 01/13/05: delete the Data class, only keep Matrix
 *   - 01/13/05: change all the overloading operators to pass-by-reference
 *   - 01/13/05: add overloading for "<<" 
 *
 ********************************************************************/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
using namespace std;

/** 
 * Matrix is an m x n matrix, where
 *     m - # of rows (or # of samples)
 *     n - # of columns (or # or features)
 *         (might have one more column that indicates the sample label/class)
 *     channel - # of channels
 **/   
class Matrix {
  friend ostream & operator<<(ostream &, Matrix &);
  friend Matrix operator/(Matrix &, double);  // matrix divided by a scalar
  friend Matrix operator*(Matrix &, double);  // matrix multiplied by a scalar
  friend Matrix operator+(Matrix &, double);  // matrix added by a scalar
  friend Matrix operator-(Matrix &, double);  // matrix subtracted by a scalar
 
 public:
  // constructors and destructor
  Matrix();                          // default constructor 
  Matrix(int,                        // constructor with row
         int);                       // column
  Matrix(const Matrix &);            // copy constructor 
  ~Matrix();                         // destructor 

  // create a matrix
  void createMatrix(int,             // row 
                    int c=1);        // column (default 1, a column vector)
  void initMatrix(float init=0.0);   // initiate the matrix component
                                     // the default is 0.0

  // get and set functions
  int getRow() const;                // get row number 
  int getCol() const;                // get column number
  void setRow(int);                  // set row number 
  void setCol(int);                  // set column number 


  // operator overloading functions
  double & operator()(int,                  // row index
		      int) const;           // column index
  const Matrix operator=(const Matrix &);   // = operator overloading
  Matrix operator+(const Matrix &) const;   // overloading + operator
  Matrix & operator+=(const Matrix &);      // overloading += operator
  Matrix operator-(const Matrix &) const;   // overloading - operator
  Matrix & operator-=(const Matrix &);      // overloading -= operator
  Matrix operator/(const Matrix &) const;   // overloading / operator
                                            // (element-wised division)
  Matrix & operator/=(const Matrix &);      // overloading /= operator
  Matrix operator*(const Matrix &) const;   // overloading * operator 
                                            // (element-wised multiplication)
  Matrix & operator*=(const Matrix &);      // overloading *= operator
  Matrix operator->*(const Matrix &) const; // overloading ->* operator 
                                            // (matrix multiplication) 
  //  Matrix & operator->*=(const Matrix &);    // overloading ->*= operator

 private:
  int row;                           // number of rows / samples 
  int col;                           // number of columns / features
  double *matrix;                    // matrix buffer
};


////////////////////////////////////
// matrix manipulation
Matrix transpose(Matrix &);        // matrix transpose
Matrix inverse(Matrix &);          // matrix inverse
Matrix mean(Matrix &, int nf);     // return a column vector with
                                   // mean value of each column
                                   // nf is the number of features
Matrix cov(Matrix &, int nf);      // return the covariance of a matrix

double det(Matrix &);              // return the determinant of a matrix
double detLU(Matrix &);            // determinant of a matrix using LU decomposition
double deteig(Matrix &);           // det of a matrix using eigenvalue decomposition
Matrix ludcmp(Matrix &, double &); // LU decomposition

Matrix subMatrix(Matrix &,         // crop a matrix
               int,                // starting row index
               int,                // starting column index
               int,                // ending row index
               int);               // ending column index
Matrix getType(Matrix &,           // get samples of a certain class
	       int);               // the data type
                                   // can only be called when the last col 
                                   // is the data type

// matrix sorting (only for 1-D row vectors)
void insertsort(Matrix &,          // the vector
		Matrix &,          // the sorted vector (row vector), 
		Matrix &);         // the index of the sorted matrix

// calculate the eigensystem
void jacobi(Matrix &,              // the input matrix
	    Matrix &,              // generated eigenvalues in vector format
	    Matrix &);             // generated eigenvectors
// sort the eigenvalues
void eigsrt(Matrix &,              // the eigenvalue vector
	    Matrix &);             // the eigenvector matrix

#endif











