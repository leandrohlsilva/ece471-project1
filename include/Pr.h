/*
 * pr.h - header file of the pattern recognition library
 *
 * Author: Hairong Qi, ECE, University of Tennessee
 *
 * Date: 01/25/04
 *
 * Please send all your comments to hqi@utk.edu 
 * 
 * Modified:
 *   - 04/26/05: reorganized for the Spring 2005 classs
 */

#ifndef _PR_H_
#define _PR_H_

#include "Matrix.h"


/////////////////////////  
// file I/O
Matrix readData(char *,            // the file name
		int);              // the number of columns of the matrix
Matrix readData(char *,            // the file name
		int,               // the number of columns
		int);              // the number of rows (or samples)
Matrix readData(char *);           // read data file to a matrix with 1 row
void writeData(Matrix &, char *);  // write data to a file
Matrix readImage(char *, int *, int *); 
                                   // read the image into a matrix (row,3col)
void writeImage(char *, Matrix &, int, int);    
                                   // write the image to a file


#endif






