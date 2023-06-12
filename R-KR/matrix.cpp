//------------------------------------------
//          Simple matrix class
//-------------------------------------------

#include "matrix.hpp"
#include "hamiltonian.hpp" 

//---------------------------------------------------
//             Creations, destructions
//---------------------------------------------------

//------------ creator ------------
Matrix::Matrix()
{
	row = 0;
	col = 0;
	coord = NULL;
}

//------------ destructor ------------
Matrix::~Matrix()
{
	if (coord != NULL)
		delete[] coord;
	coord = NULL;
	row = 0;
	col = 0;
}

//------------ initialization ------------
Matrix::Matrix(int Row, int Col)
{
	row = Row;
	col = Col;
	coord = new double[Row * Col];
}

void Matrix::set_size(int Row, int Col)
{
	row = Row;
	col = Col;
	if (coord != NULL)
		delete[] coord;
	coord = new double[col * row];
}

//------------ creation same matrice ------------
Matrix &Matrix::operator=(const Matrix &a)
{
	if (this != &a)
	{
		row = a.row;
		col = a.col;
		if (coord != NULL)
			delete[] coord;
		if (a.coord != NULL)
		{
			coord = new double[row * col];
			for (int i = 0; i < row * col; i++)
				coord[i] = a.coord[i];
		}
		else
		{
			coord = NULL;
		}
	}
	return *this;
}

//------------ set to zero ------------
void Matrix::zeros()
{
	if (coord != NULL)
	{
		for (int i = 0; i < col * row; i++)
			coord[i] = 0.0;
	}
}

//------------ set to zero ------------
void Matrix::identity()
{
	if (coord != NULL)
	{
		zeros();
		if (row == col)
		{
			for (int i = 0; i < col; i++)
			{
				(*this)(i, i) = 1.0;
			}
		}
		else
		{
			cout << "ERROR: identity() only works for a square matrix.  "
				 << "Called on a " << row << " by " << col << " matrix." << endl;
			exit(-1);
		}
	}
}

// Output a square submatrix with indices [i1:i2, j1:j2].
Matrix Matrix::subslice(int i1, int i2, int j1, int j2)
{
	Matrix ret;
	if (i2 >= i1 && j2 >= j1)
	{
		ret.set_size(i2 - i1 + 1, j2 - j1 + 1);
		for (int i = i1; i <= i2; i++)
		{
			for (int j = j1; j <= j2; j++)
			{
				ret(i - i1, j - j1) = (*this)(i, j);
			}
		}
	}
	else
	{
		cout << "For subslice, must have i1<=i2 j1<=j2.  (i1,i2,j1,j2): ("
			 << i1 << ", " << i2 << ", " << j1 << ", " << j2 << ")" << endl;
		exit(1);
	}
	return ret;
}

//------------ access an element ------------
double &Matrix::operator()(int a, int b)
{
	if (coord != NULL)
	{
		if (a >= row || a < 0 || b >= col || b < 0)
		{
			cerr << "Array Out of Bounds! " << endl;
			exit(1);
		}
		return coord[a * col + b];
	}
	else
	{
		cerr << " ### ERROR ### requiring the component of an empty matrix! Requesting (" << a << ", " << b << ")" << endl;
		exit(1);
	}
}

//---------------------------------------------------
//        Overloaded operations
//---------------------------------------------------

//------------ addition ------------
Matrix Matrix::operator+(const Matrix &a)
{
	Matrix ret;
	if ((row != a.row) && (col != a.col))
	{
		cerr << " ### ERROR ### adding matrices of differents sizes! " << endl;
		exit(1);
	}
	else
	{
		if ((coord != NULL) && (a.coord != NULL))
		{
			ret.row = a.row;
			ret.col = a.col;
			ret.coord = new double[col * row];
			for (int i = 0; i < col * row; i++)
				ret.coord[i] = coord[i] + a.coord[i];
		}
	}
	return ret;
}

Matrix &Matrix::operator+=(const Matrix &a)
{
	Matrix ret;
	if ((row != a.row) && (col != a.col))
	{
		cerr << " ### ERROR ### adding matrices of differents sizes! " << endl;
		exit(1);
	}
	else
	{
		if ((coord != NULL) && (a.coord != NULL))
		{
			for (int i = 0; i < col * row; i++)
				coord[i] += a.coord[i];
		}
	}
	return *this;
}

Matrix &Matrix::operator-=(const Matrix &a)
{
	Matrix ret;
	if ((row != a.row) && (col != a.col))
	{
		cerr << " ### ERROR ### adding matrices of differents sizes! " << endl;
		exit(1);
	}
	else
	{
		if ((coord != NULL) && (a.coord != NULL))
		{
			for (int i = 0; i < col * row; i++)
				coord[i] -= a.coord[i];
		}
	}
	return *this;
}

//------------ subtraction ------------
Matrix Matrix::operator-(const Matrix &a)
{
	Matrix ret;
	if ((row != a.row) && (col != a.col))
	{
		cerr << " ### ERROR ### substracting matrices of differents sizes! " << endl;
		exit(1);
	}
	else
	{
		if ((coord != NULL) && (a.coord != NULL))
		{
			ret.row = a.row;
			ret.col = a.col;
			ret.coord = new double[col * row];
			for (int i = 0; i < col * row; i++)
				ret.coord[i] = coord[i] - a.coord[i];
		}
	}
	return ret;
}

//------------ multiplication by a constant ------------
Matrix Matrix::operator*(const double c)
{
	Matrix ret;
	if (coord != NULL)
	{
		ret.row = row;
		ret.col = col;
		ret.coord = new double[col * row];
		for (int i = 0; i < col * row; i++)
			ret.coord[i] = coord[i] * c;
	}
	return ret;
}

//------------ multiplication by matrix a ------------
Matrix Matrix::operator*(const Matrix &a)
{
	Matrix ret;
	if (a.row != col)
	{
		cout << "Incompatible matrix multiplication "
			 << row << "x" << col << " * " << a.row << "x" << a.col << endl;
		exit(-1);
	}
	ret.set_size(row, a.col);
	ret.zeros();
	for (int i = 0; i < row; i++)
		for (int j = 0; j < a.col; j++)
			for (int k = 0; k < col; k++)
				ret(i, j) = ret(i, j) + coord[i * col + k] * a.coord[k * a.col + j];
	return ret;
}

Matrix Matrix::transpose()
{
	Matrix ret;
	if (coord != NULL)
	{
		ret.set_size(col, row);
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
			{
				ret(j, i) = (*this)(i, j);
			}
		}
	}
	return ret;
}

double Matrix::contraction_product(Matrix &m)
{
	double tot = 0;
	if (col == m.col && row == m.row)
	{
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
			{
				tot += (*this)(i, j) * m(i, j);
			}
		}
	}
	else
	{
		cout << "ERROR contraction_product() matrices need to be same size.  "
			 << " this " << row << " by " << col
			 << " that " << row << " by " << col << endl;
		exit(-1);
	}
	return tot;
}

ostream &operator<<(ostream &output, const Matrix &m)
{ 
	if (m.coord != NULL)
	{
		for (int i = 0; i < m.row; i++)
		{
			for (int j = 0; j < m.col; j++)
			{
				output << setw(14) << m.coord[i * m.col + j];
				output << "  ";
			} 
		}
	}
	else
	{
	} 
	return output; // for multiple << operators.
}

istream &operator>>(istream &input, Matrix &m)
{
	if (m.coord != NULL)
		delete[] m.coord;
	input >> m.row >> m.col;
	m.coord = new double[m.col * m.row];
	for (int i = 0; i < m.row; i++)
	{
		for (int j = 0; j < m.col; j++)
		{
			input >> m.coord[i * m.col + j];
		}
	}
	return input; // for multiple >> operators.
}

//  Computing inverses is ill-posed.  Why do I need to do it?
Matrix Matrix::inverse()
{
	Matrix ret;
	if (col == row)
	{
		ret.set_size(row, col);
		if (col == 1)
		{
			ret(0, 0) = 1 / coord[0];
		}
		if (col == 2)
		{
			double det = coord[0] * coord[3] - coord[1] * coord[2];
			if (fabs(det) > numeric_limits<double>::epsilon())
			{
				ret.coord[0] = coord[3] / det;
				ret.coord[1] = -coord[1] / det;
				ret.coord[2] = -coord[2] / det;
				ret.coord[3] = coord[0] / det;
			}
			else
			{
				cout << "Error: inverse called on non-invertible matrix."
					 << "  Determinant " << det << endl;
				exit(1);
			}
		}
		else
		{
			cout << "Error: inverse not implemented beyond 2 by 2." << endl;
			exit(1);
		}
	}
	else
	{
		cout << "Error: inverse called for non-square matrix." << endl;
		exit(1);
	}
	return ret;
}

// For a vector, returns a diagonal matrix that is the
// exponential of the vector.  Does not work for general cases.
Matrix Matrix::exp_diag()
{
	Matrix ret;
	if (row == 1 || col == 1)
	{
		if (col > row)
		{
			ret = this->transpose().exp_diag();
		}
		else
		{
			ret.set_size(row, row);
			ret.zeros();
			for (int i = 0; i < row; i++)
			{
				ret(i, i) = exp((*this)(i, 0));
			}
		}
	}
	else
	{
		cout << "exp_diag() called for " << row << " by " << col << ".  Only "
			 << "implemented for vectors." << endl;
		exit(-1);
	}
	return ret;
}

Matrix Matrix::crosses()
{
	Matrix ret;
	if (coord != NULL && col == 3 && row == 3)
	{
		ret.set_size(3, 3);
		int k, l;
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
			{
				ret(i, j) = (*this)((i + 1) % 3, (j + 1) % 3) * (*this)((i + 2) % 3, (j + 2) % 3) - (*this)((i + 1) % 3, (j + 2) % 3) * (*this)((i + 2) % 3, (j + 1) % 3);
			}
		}
	}
	else
	{
		cout << "crosses() computes cross products for a 3 by 3 matrix only.  "
			 << "It was called with " << row << " by " << col
			 << " matrix " << endl;
		exit(-1);
	}
	return ret;
}

double Matrix::column_norm(int j)
{
	double tot = 0;
	if (coord != NULL)
	{
		for (int i = 0; i < row; i++)
		{
			tot += pow(coord[i * col + j], 2);
		}
	}
	else
	{
		tot = -1;
	}
	return sqrt(tot);
}

double Matrix::norm()
{
	double tot = 0;
	if (coord != NULL)
	{
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
			{
				tot += pow(coord[i * col + j], 2);
			}
		}
	}
	else
	{
		tot = -1;
	}
	return sqrt(tot);
}
