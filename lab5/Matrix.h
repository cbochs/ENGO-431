#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

class Matrix
{
	public:
	Matrix();										
	Matrix(unsigned int rows,unsigned int cols);
	Matrix(unsigned int rows,unsigned int cols,double initial);
	void resize(unsigned int rows, unsigned int cols);
	//First at function is used for both reading and assigning a element of the matrix
	double& at(unsigned int rows, unsigned int cols);	
	//The extra at function is required to access members of a const Matrix
	//note that the second at function does not a return a reference, so it can't be used 
	//to change a value, it can only return a value.
	double at(unsigned int rows, unsigned int cols) const;

	//Assignment operator
	const Matrix& operator= (const Matrix& mat1);
	
	//operator to return a reference to the "row_th" row of the private data
	vector<double>& operator[](unsigned int row);

	//overload that works on const matrix objects (it returns the value, not the reference
	vector<double> operator[](unsigned int row) const;

	unsigned int getrows() const;  //const is required so it will work on a const Matrix as well as a Matrix
	unsigned int getcols() const;

	Matrix trans();	       // Return the transpose of the matrix
	Matrix inv();	       // Return the inverse of the matrix
	double det();	       // Return the determinant of the a positive definite square matrix
	double maxAbsElem();   // Returns the maximum absolute value in the matrix
	double maxElem();      // Returns the maximum value in the matrix
	double minElem();      // Returns the minimum value in the matrix
			
	// Print the matrix values to a specified stream 
	// mode:
	// ios_base::out [default] --> to overwrite the file if it exists
	// ios_base::app           --> to append to the file if it exists
	void print(const char *filename,int prec = 3, int width = 12,const char *header = "",
		const char *footer = "",ios_base::open_mode mode = ios_base::out );

	void print(int prec = 3, int width = 12,const char *header = "",
        const char *footer = ""); //Print matrix to the screen
	
	void read(const char* filename);  //Read a file into the matrix where the file is given by filename
	void read(const char* filename, unsigned int rows , unsigned int cols);  //Read a file into the matrix where the file is given by filename
	

	void sortCol(unsigned int n=0); // Sort rows in ascending order according to the nth column [n= 0 to noCol-1 ]
	void clear (); // Zero all elements

	friend const Matrix operator* (const Matrix& mat1, const Matrix& mat2);
	friend const Matrix operator* (int a, const Matrix& mat2);
	friend const Matrix operator* (const Matrix& mat1, int b);
	friend const Matrix operator* (double a, const Matrix& mat2);
	friend const Matrix operator* (const Matrix& mat1, double b);
    
	friend const Matrix operator+ (const Matrix& mat1, const Matrix& mat2);
	friend const Matrix operator+ (int a, const Matrix& mat2);
	friend const Matrix operator+ (const Matrix& mat1, int b);
	friend const Matrix operator+ (double a, const Matrix& mat2);
	friend const Matrix operator+ (const Matrix& mat1, double b);
    
	friend const Matrix operator- (const Matrix& mat1, const Matrix& mat2);
	friend const Matrix operator- (int a, const Matrix& mat2);
	friend const Matrix operator- (const Matrix& mat1, int b);
	friend const Matrix operator- (double a, const Matrix& mat2);
	friend const Matrix operator- (const Matrix& mat1, double b);

	friend const Matrix operator/ (const Matrix& mat1, double b);
	friend const Matrix operator/ (const Matrix& mat1, int b);
	
	friend const Matrix operator/ (const Matrix& mat1, const Matrix& mat2);

	Matrix exclude(int row, int col);
	double determinant();
	double determinant(Matrix *mat);
	Matrix inverse();

private:
	unsigned int n_rows;
	unsigned int n_cols;
	vector<vector<double> > data;

};


#endif