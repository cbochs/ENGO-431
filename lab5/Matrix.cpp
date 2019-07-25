//Implemenation of a matrix class, Lab 4 version;
//Author Kyle O'Keefe

#include "Matrix.h"

//Constructor: Used to set size to zero
Matrix::Matrix()
{
	n_rows = 0;
	n_cols = 0;
	data.resize(0);  //this line not necessary as the vector constructor will do this.
}

//Overload of constructor, to declare a matrix with dimensions
Matrix::Matrix(unsigned int rows,unsigned int cols)
{
	n_rows = rows;
	n_cols = cols;
	data.resize(rows);
	for (unsigned int i = 0; i < n_rows; i++)
		data[i].resize(n_cols);
	
	return;
}

Matrix::Matrix(unsigned int rows,unsigned int cols, double initial)
{
	n_rows = rows;
	n_cols = cols;
	data.resize(rows);
	for (unsigned int i = 0; i < n_rows; i++){
		data[i].resize(n_cols);
		data[i].assign(data[i].size(),initial);}
	
	return;
}


//Resize function, uses the resize function of vector
void Matrix::resize(unsigned int rows, unsigned int cols)
{
	n_rows = rows;
	n_cols = cols;
	data.resize(rows);
	for (unsigned int i = 0; i < n_rows; i++)
		data[i].resize(n_cols);
	return;
}

double& Matrix::at(unsigned int rows, unsigned int cols)
{
	if(rows < n_rows && cols < n_cols)
		return data.at(rows).at(cols);
	else
	{
		cout << "Error: Matrix::at Index requested exceeds matrix dimensions" << endl;
		exit(1);
	}
}

double Matrix::at(unsigned int rows, unsigned int cols) const
{
	if(rows < n_rows && cols < n_cols)
		return data.at(rows).at(cols);
	else
	{
		cout << "Error: Matrix::at Index requested exceeds matrix dimensions" << endl;
		exit(1);
	}
}


vector<double>& Matrix::operator[](unsigned int row)
{
	if(row < n_rows)
		return data.at(row);
	else
	{
		cout << "Error: Matrix::[] Index requested exceeds number of rows" << endl;
		exit(1);
	}
}
vector<double> Matrix::operator[](unsigned int row) const
{
	if(row < n_rows)
		return data.at(row);
	else
	{
		cout << "Error: Matrix::[] Index requested exceeds number of rows" << endl;
		exit(1);
	}
}

const Matrix& Matrix::operator= (const Matrix& mat1)
	{
		if(this == &mat1)
			return *this;

		this->resize(mat1.n_rows,mat1.n_cols);
		
		for (unsigned int i = 0; i < mat1.n_rows; i++)
               for (unsigned int j = 0; j < mat1.n_cols; j++)
                       this->at(i,j) = mat1.at(i,j);

		return *this;
	}
unsigned int Matrix::getrows() const
{
	return this->n_rows;
}
unsigned int Matrix::getcols() const
{
	return this->n_cols;
}
void Matrix::sortCol(unsigned int n){
	
	if ( n >= n_cols ){
		cout << "Error: Matrix::sortCol Index requested exceeds matrix dimensions" << endl;
		exit(1);
	}
    int flag = 1;     // set flag to 1 to begin initial pass
    double Temp = 0;  // holding variable
    
	for(unsigned int i = 0 ; (i < n_rows) && flag; i++)
     {
          flag = 0;
          for (unsigned int j=0; j < (n_rows -1); j++)
         {
               if (this->at(j+1,n) < this->at(j,n))      // descending order simply changes to >
              { 
                       Temp = this->at(j,n);             // swap elements
                       this->at(j,n) = this->at(j+1,n);
                       this->at(j+1,n) = Temp;

					   for( unsigned int k=0; k <n_cols ;k++)
						   if( k == n)
							   continue;
						   else{
							   Temp = this->at(j,k);             // swap elements
							   this->at(j,k) = this->at(j+1,k);
							   this->at(j+1,k) = Temp;

						   }
    
					   
					   flag = 1;               // indicates that a swap occurred.
               }
		  }
	}
 }
Matrix Matrix::trans()
{
	Matrix temp(n_cols,n_rows);
	for (unsigned int i = 0; i < n_rows; i++)
		for (unsigned int j = 0; j < n_cols; j++)
			temp.at(j,i) = data.at(i).at(j);
			//temp.at(j,i) = this->at(i,j);  //alternative way to access data using Matrix::at
	return temp;
}

Matrix Matrix::inv()
{
	const double EPSILON = 1.0E-12;

       if (n_rows == 0 || n_cols == 0)
	   {
		   cout << "Error: Matrix::inv Inverse called on matrix with zero rows and/or columns" << endl;
		   exit(1);
	   }

       if (n_rows != n_cols)
	   {
		   cout << "Error: Matrix::inv Non-square matrix in Cholesky Decomposition" << endl;
		   exit(1);
	   }

	  
	   // temp object to hold the result
       Matrix A (n_rows, n_cols); // store result

	   for (unsigned int i = 0; i < this->n_rows; i++)
               for (unsigned int j = 0; j < this->n_cols; j++)
                       A.at(i,j) = this->at(i,j);
	  
	   // check for positive definiteness
       for (unsigned int i = 0; i < n_rows; i++)
       {
               if (A.at(i,i) < 0.0)
			   {
				   cout << "Error: Matrix::inv Matrix not positive definite in Cholesky Decomposition" << endl;
				   exit(1);
			   }

               if (fabs (A.at(i,i)) < EPSILON)
			   {
				   cout << "Error: Matrix::inv Singular matrix in Cholesky Decomposition" << endl;
				   exit(1);
			   }
       }


       // Perform Choleski decomposition
       for (unsigned int j = 0; j < n_rows; j++)
       {
               for(unsigned int k = 0; k < j; k++)
                       A.at(j,j) -= A.at(j,k) * A.at(j,k);

               if (A.at(j,j) < 0.0)
			   {
				   cout << "Error: Matrix::inv Square root of negative number in Cholesky Decomposition" << endl;
				   exit(1);
			   }

               A.at(j,j) = sqrt (A.at(j,j));
			  
               for(unsigned int i = j + 1; i < n_rows; i++)
               {
                       for (unsigned int k = 0; k < j; k++)
                               A.at(i,j) -= A.at(i,k) * A.at(j,k);

                       if (fabs(A.at(j,j)) < EPSILON)
					   {
						   cout << "Error: Matrix::inv Division by zero in Cholesky Decomposition" << endl;
						   exit(1);
					   }

                       A.at(i,j) /= A.at(j,j);

               }
       }

       // inversion of lower triangular matrix
       for (unsigned int j = 0; j < n_rows; j++)
       {
               A.at(j,j) = 1.0 / A.at(j,j);

               for (unsigned int i = j + 1; i < n_rows; i++)
               {
                       A.at(i,j) *= -A.at(j,j) / A.at(i,i);

                       for (unsigned int k = j + 1; k < i; k++)
                               A.at(i,j) -= A.at(i,k) * A.at(k,j) / A.at(i,i);
               }
       }

       // construction of lower triangular inverse matrix
       for (unsigned int j = 0; j < n_rows; j++)
       {
               for (unsigned int i = j; i < n_rows; i++)
               {
                       A.at(i,j) *= A.at(i,i);

                       for (unsigned int k = i + 1; k < n_rows; k++)
                               A.at(i,j) += A.at(k,i) * A.at(k,j);
               }
       }

       // fill upper diagonal
       for (unsigned int i = 1; i < n_rows; i++)
       {
               for (unsigned int j = 0; j < i; j++)
                       A.at(j,i) = A.at(i,j);
       }
       return A; // return a copy
}

double Matrix::maxAbsElem()
{
	if (n_rows == 0 || n_cols == 0)
	{
		cout <<"MaxAbsElem called on matrix with zero rows and/or columns"<<endl;
		exit(1);
	}

       double max = data[0][0];
       for (unsigned int i = 0; i < n_rows; i++)
       {
               for (unsigned int j = 0; j < n_cols; j++)
               {
                       if ( fabs( data[i][j] ) > max )
                               max = fabs(data[i][j]);
               }
       }
       return max;
}

double Matrix::maxElem(){

	if (n_rows == 0 || n_cols == 0)
	{
		cout <<"MaxElem called on matrix with zero rows and/or columns"<<endl;
		exit(1);
	}

       double max = data[0][0];
       for (unsigned int i = 0; i < n_rows; i++)
       {
               for (unsigned int j = 0; j < n_cols; j++)
               {
                       if (  data[i][j]  > max )
                               max = data[i][j];
               }
       }
       return max;

}

double Matrix::minElem(){

	cout<<n_rows<<endl<<n_cols<<endl;
	if (n_rows == 0 || n_cols == 0)
	{
		cout <<"MaxElem called on matrix with zero rows and/or columns"<<endl;
		exit(1);
	}

       double min = data[0][0];
       for (unsigned int i = 0; i < n_rows; i++)
       {
               for (unsigned int j = 0; j < n_cols; j++)
               {
                       if (  data[i][j]  < min )
                               min = data[i][j];
               }
       }
       return min;
}




double Matrix::det()
{
	const double EPSILON = 1.0E-12;

       if (n_rows == 0 || n_cols == 0)
	   {
		   cout << "Error: Matrix::det Inverse called on matrix with zero rows and/or columns" << endl;
		   exit(1);
	   }

       if (n_rows != n_cols)
	   {
		   cout << "Error: Matrix::det Non-square matrix in Cholesky Decomposition" << endl;
		   exit(1);
	   }

	  
	   // temp object to hold the result
       Matrix A (n_rows, n_cols); // store result

	   for (unsigned int i = 0; i < this->n_rows; i++)
               for (unsigned int j = 0; j < this->n_cols; j++)
                       A.at(i,j) = this->at(i,j);
	  
	   // check for positive definiteness
       for (unsigned int i = 0; i < n_rows; i++)
       {
               if (A.at(i,i) < 0.0)
			   {
				   cout << "Error: Matrix::inv Matrix not positive definite in Cholesky Decomposition" << endl;
				   exit(1);
			   }

               if (fabs (A.at(i,i)) < EPSILON)
			   {
				   cout << "Error: Matrix::inv Singular matrix in Cholesky Decomposition" << endl;
				   exit(1);
			   }
       }

	   double det = 1.0;
       // Perform Choleski decomposition
       for (unsigned int j = 0; j < n_rows; j++)
       {
               for(unsigned int k = 0; k < j; k++)
                       A.at(j,j) -= A.at(j,k) * A.at(j,k);

               if (A.at(j,j) < 0.0)
			   {
				   cout << "Error: Matrix::inv Square root of negative number in Cholesky Decomposition" << endl;
				   exit(1);
			   }

               A.at(j,j) = sqrt (A.at(j,j));
			   det *= A.at(j,j) * A.at(j,j);

               for(unsigned int i = j + 1; i < n_rows; i++)
               {
                       for (unsigned int k = 0; k < j; k++)
                               A.at(i,j) -= A.at(i,k) * A.at(j,k);

                       if (fabs(A.at(j,j)) < EPSILON)
					   {
						   cout << "Error: Matrix::inv Division by zero in Cholesky Decomposition" << endl;
						   exit(1);
					   }

                       A.at(i,j) /= A.at(j,j);
               }
       }

	   return det;
}

void Matrix::print(int prec, int width, const char *header,const char *footer )
{

	   // set output format
       cout.precision (prec);
       cout.setf (std::ios::showpoint);
       cout.setf (std::ios::fixed, std::ios::floatfield);

       // overflow in fixed format can cause runtime error
       // switch to scientific format if req'd
       if (maxAbsElem() > 1.0E+15) 
		   cout.setf (std::ios::scientific, std::ios::floatfield);

	   cout<< header <<endl ;

	for (unsigned int i = 0; i < n_rows; i++)
	{
		for (unsigned int j = 0; j < n_cols; j++)
			cout << setw(width)<<data[i][j] << flush;
		if ( i != n_rows-1)
			cout << endl; //carriage return after the end of the ith row
	}
	
	cout << "\t"<<footer << endl <<endl <<endl ;

	cout.unsetf ( ios::fixed );
	cout.setf(0,ios::floatfield);
	return;
}

void Matrix::print(const char *filename,int prec, int width,const char *header,const char *footer ,ios_base::open_mode mode )
{
	// mode can be on of the following:
	// ios_base::app  --> Append to the file if it exists
	// ios_base::out  --> Overwrite the file if it exists 
	


	ofstream outfile;
	outfile.open  (filename, mode);

	if (!outfile.is_open())
	{
		cout << "Error: print - Error opening file" << endl;
		exit(1);
	}

//	 set output format
       outfile.precision (prec);
       outfile.setf (std::ios::showpoint);
       outfile.setf (std::ios::fixed, std::ios::floatfield);

//	 overflow in fixed format can cause runtime error
 //   // switch to scientific format if req'd
	if (maxAbsElem() > 1.0E+15) 
		cout.setf (std::ios::scientific, std::ios::floatfield);
	
	outfile << header << endl ;

	for (unsigned int i = 0; i < n_rows; i++)
	{
			for (unsigned int j = 0; j < n_cols; j++)
				outfile << setw(width)<<data[i][j];// << flush;
	
			if ( i != n_rows-1)
				outfile << endl; //carriage return after the end of the ith row
							
   }
	
	outfile << "\t"<<footer << endl <<endl <<endl ;

	outfile.unsetf ( ios::fixed );
	outfile.setf(0,ios::floatfield);

	return; 

}

void Matrix::read(const char *filename)
{
	//Define and open the input file stream
	ifstream infile;                    
	infile.open  (filename, ifstream::in);

	if (!infile.is_open())
	{
		cout << "Error: read - Error opening file" << endl;
		exit(1);
	}

	// Read noRows and noCols from the text file
	infile >> n_rows >> n_cols;

	this->resize(n_rows,n_cols);

	// Define matrix to store data in
	Matrix temp(n_rows,n_cols);

	//Read in the values from the text file and store them into Read Matrix
	for(unsigned int i=0; i<n_rows;i++)
		for(unsigned int j=0;j<n_cols;j++)
			infile>> this->at(i,j);
		
    // Close file
	infile.close();
}


void Matrix::read(const char *filename, unsigned int rows , unsigned int cols)
{
	//Define and open the input file stream
	ifstream infile;                    
	infile.open  (filename, ifstream::in);

	if (!infile.is_open())
	{
		cout << "Error: read - Error opening file" << endl;
		exit(1);
	}

	// Read noRows and noCols from the text file
	n_rows = rows;
	n_cols = cols;

	this->resize(n_rows,n_cols);

	// Define matrix to store data in
	Matrix temp(n_rows,n_cols);

	//Read in the values from the text file and store them into Read Matrix
	for(unsigned int i=0; i<n_rows;i++)
		for(unsigned int j=0;j<n_cols;j++)
			infile>> this->at(i,j);
		
    // Close file
	infile.close();
}

	
// set all matrix elements to zero	
void Matrix::clear()
{
       for (unsigned int i = 0; i < n_rows; i++)
               for (unsigned int j = 0; j < n_cols; j++)
                       data[i][j] = 0.0;

       return;
}

const Matrix operator+(const Matrix& mat1, const Matrix& mat2)
{
	if (mat1.n_rows != mat2.n_rows || mat1.n_cols != mat2.n_cols)
	{
		cout << "Error: operator+ Matrices not conformal for addition" << endl;
		exit(1);
	}

	Matrix temp(mat1.n_rows, mat1.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = mat1.at(i,j) + mat2.at(i,j);

	return temp;
}

const Matrix operator+ (int a, const Matrix& mat2)
{
	Matrix temp(mat2.n_rows, mat2.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = a + mat2.at(i,j);

	return temp;
}

const Matrix operator+ (const Matrix& mat1, int b)
{
	return b + mat1;
}

const Matrix operator+ (double a, const Matrix& mat2)
{
	Matrix temp(mat2.n_rows, mat2.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = a + mat2.at(i,j);

	return temp;
}

const Matrix operator+ (const Matrix& mat1, double b)
{
	return b + mat1;
}

const Matrix operator-(const Matrix& mat1, const Matrix& mat2)
{
	if (mat1.n_rows != mat2.n_rows && mat1.n_cols != mat2.n_cols)
	{
		cout << "Error: operator- Matrices not conformal for addition" << endl;
		exit(1);
	}

	else if(mat1.n_rows == mat2.n_rows && mat1.n_cols == mat2.n_cols ){

	Matrix temp(mat1.n_rows, mat1.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = mat1.at(i,j) - mat2.at(i,j);
	return temp;}

	else if(mat1.n_rows == mat2.n_rows && mat2.n_cols == 1 ){
		Matrix temp(mat1.n_rows, mat1.n_cols);
		for (unsigned int i = 0; i < temp.n_rows; i++)
			for (unsigned int j = 0; j < temp.n_cols; j++)
				temp.at(i,j) = mat1.at(i,j) - mat2.at(i,0);
	return temp;}

	else if(mat1.n_rows == mat2.n_rows && mat1.n_cols == 1 ){
		Matrix temp(mat2.n_rows, mat2.n_cols);
		for (unsigned int i = 0; i < temp.n_rows; i++)
			for (unsigned int j = 0; j < temp.n_cols; j++)
				temp.at(i,j) = mat1.at(i,0) - mat2.at(i,j);
	return temp;}

	else if(mat1.n_cols == mat2.n_cols && mat2.n_rows == 1 ){
		Matrix temp(mat1.n_rows, mat1.n_cols);
		for (unsigned int i = 0; i < temp.n_rows; i++)
			for (unsigned int j = 0; j < temp.n_cols; j++)
				temp.at(i,j) = mat1.at(i,j) - mat2.at(0,j);
	return temp;}

	else if(mat1.n_cols == mat2.n_cols && mat1.n_rows == 1 ){
		Matrix temp(mat2.n_rows, mat2.n_cols);
		for (unsigned int i = 0; i < temp.n_rows; i++)
			for (unsigned int j = 0; j < temp.n_cols; j++)
				temp.at(i,j) = mat1.at(0,j) - mat2.at(0,j);
	return temp;}

	// temp .. Check why I did the above
	// the exit is needed here just to get rid of the warning:
	// warning C4715: 'operator-' : not all control paths return a value

	exit(1);
	
}

const Matrix operator- (int a, const Matrix& mat2)
{
	Matrix temp(mat2.n_rows, mat2.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = a - mat2.at(i,j);

	return temp;
}

const Matrix operator- (const Matrix& mat1, int b)
{
	Matrix temp(mat1.n_rows, mat1.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = mat1.at(i,j) - b;

	return temp;
}

const Matrix operator- (double a, const Matrix& mat2)
{
	Matrix temp(mat2.n_rows, mat2.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = a - mat2.at(i,j);

	return temp;
}

const Matrix operator- (const Matrix& mat1, double b)
{
	Matrix temp(mat1.n_rows, mat1.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = mat1.at(i,j) - b;

	return temp;
}


const Matrix operator* (const Matrix& mat1, const Matrix& mat2)
{
	if (mat1.n_cols == 0 || mat1.n_rows == 0 ||
               mat2.n_cols == 0 || mat2.n_rows == 0)
       {
		   cout << "Error: operator* Multiplication called on matrix with zero rows and/or columns" << endl;
		   exit(1);
               
       }

	if (mat1.n_cols != mat2.n_rows)
	{
		 cout << "Error: operator* Matrices not conformal for multiplication" << endl;
		 exit(1);
	}

	Matrix temp(mat1.n_rows, mat2.n_cols);
	for (unsigned int i = 0; i < mat1.n_rows; i++)
		for (unsigned int j = 0; j < mat2.n_cols; j++)
			for (unsigned int k = 0;k < mat1.n_cols; k++)
				temp[i][j] += mat1[i][k] * mat2[k][j];

	return temp;
}

const Matrix operator* (int a, const Matrix& mat2)
{
	Matrix temp(mat2.n_rows, mat2.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = a * mat2.at(i,j);

	return temp;
}

const Matrix operator* (const Matrix& mat1, int b)
{
	Matrix temp(mat1.n_rows, mat1.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = mat1.at(i,j) * b;

	return temp;
}

const Matrix operator* (double a, const Matrix& mat2)
{
	Matrix temp(mat2.n_rows, mat2.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = a * mat2.at(i,j);

	return temp;
}

const Matrix operator* (const Matrix& mat1, double b)
{
	Matrix temp(mat1.n_rows, mat1.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = mat1.at(i,j) * b;

	return temp;
}

const Matrix operator/ (const Matrix& mat1, int b)
{
	Matrix temp(mat1.n_rows, mat1.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = mat1.at(i,j) / b;

	return temp;
}

const Matrix operator/ (const Matrix& mat1, double b)
{
	Matrix temp(mat1.n_rows, mat1.n_cols);
	for (unsigned int i = 0; i < temp.n_rows; i++)
		for (unsigned int j = 0; j < temp.n_cols; j++)
			temp.at(i,j) = mat1.at(i,j) / b;

	return temp;
}

const Matrix operator/ (const Matrix& mat1, const Matrix& mat2)
{

	if (mat1.n_rows !=  mat2.n_rows )
	{
		cout << "Error: operator/ LSQR : The dimensions of the design matrix and observation vector do not match" << endl;
		   exit(1);
               
    }

	Matrix A = mat1;
	Matrix y = mat2;
	Matrix x(A.getcols(),1);

	Matrix N = A.trans() * A;
	Matrix N_inv = N.inv();

	double temp = 0.0;
	unsigned int count = 0;
	Matrix d;
	do{
		count = count + 1;
		Matrix U = A.trans() * (A*x - y);
		d = -1* N_inv * U ;
		x = x + d ;
		temp = d.maxAbsElem();} while(temp > 1e-9 && count <20 );

	
	return x;

}

Matrix Matrix::exclude(int row, int col)
{
	Matrix A(this->n_rows - 1, this->n_cols - 1);
	for(int y = 0; y < this->n_rows; y++) 
	{
		int shift_x = 0;
		int shift_y = 0;
		for(int x = 0; x < n_cols; x++)
		{
			if(y == row || x == col) continue;
			if(y > row) shift_y = -1;
			if(x > col) shift_x = -1;
				
			A.at(y + shift_y, x + shift_x) = this->at(y,x);
		}
	}
	
	return A;
}

double Matrix::determinant()
{
	return determinant(this);
}

double Matrix::determinant(Matrix *mat)
{
	int rows = mat->getrows();
	int cols = mat->getcols();
	if(rows != cols)
	{
		cout << "Not Square" << endl;
		exit(1);
	} 
	else if(rows == 2)
	{
		return mat->at(0,0)*mat->at(1,1)-mat->at(0,1)*mat->at(1,0);
	}

	double det = 0;
	for(int i = 0; i < mat->getrows(); i++)
	{
		det += pow(-1., i) * mat->at(0,i) * determinant(&(mat->exclude(i,0)));
	}

	return det;
}

Matrix Matrix::inverse()
{
	double det = this->determinant();

	if(det == 0)
	{
		cout << "Cannot compute inverse" << endl;
		exit(1);
	}

	Matrix adjugate(this->n_rows, this->n_cols);

	for(int y = 0; y < this->n_rows; y++)
	{
		for(int x = 0; x < this->n_cols; x++)
		{
			adjugate.at(y,x) = exclude(y,x).determinant();
		}
	}

	return (1/det) * adjugate;
}