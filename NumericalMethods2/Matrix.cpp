#include <iostream>
#include <iomanip>
#include "Matrix.h"
#include "MathOper.h"

// конструктор по умолчанию
Matrix::Matrix(): Matrix(3, 3) {}

// конструктор с параметрами
Matrix::Matrix(int lines, int columns) : n(lines), m(columns)
{
	for (size_t i = 0; i < n; ++i)
	{
		std::vector<double> temp;
		for (size_t j = 0; j < m; ++j)
			temp.push_back(0);
		A.push_back(temp);
	}
}

// конструктор копирования
Matrix::Matrix(const Matrix& matrixCopy): n(matrixCopy.GetLines()), m(matrixCopy.GetColumns())
{
	for (size_t i = 0; i < n; ++i)
	{
		A.push_back(matrixCopy[i]);
	}
}

// конструктор для двумерного массива
Matrix::Matrix(std::vector<std::vector<double>> twoDimMass): n(twoDimMass.size()), m(twoDimMass[0].size())
{
	for (size_t i = 0; i < n; ++i)
	{
		A.push_back(twoDimMass[i]);
	}
}

Matrix Matrix::Ones(int n)
{
	Matrix result(n, n);
	for (int i = 0; i < n; ++i)
	{
		result[i][i] = 1;
	}
	return result;
}

Matrix Matrix::Ort(int n)
{
	Matrix result(n, 1);
	result[0][0] = 1;
	return result;
}

Matrix Matrix::FormAForTest5(int n, int N, double e)
{
	Matrix result(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (i < j) result[i][j] = -1 - e * N;
			else if (i > j) result[i][j] = e * N;
			else result[i][j] = 1 + e * N;
		}
	}
	return result;
}

Matrix Matrix::FormbForTest5(int n)
{
	Matrix result(n, 1);
	for (int i = 0; i < n-1; ++i)
	{
		result[i][0] = -1;
	}
	result[n-1][0] = 1;
	return result;
}

// перегрузка оператора индексации для обращения к строке
std::vector<double> Matrix::operator[] (int i) const
{
	if (i < 0 || i >= n)
	{
		std::cerr << "\nIndex error: " << i << std::endl;
		exit(2); // завершить работу программы, неправильный индекс
	}
	return A[i];
}

std::vector<double>& Matrix::operator[] (int i)
{
	if (i < 0 || i >= n)
	{
		std::cerr << "\nIndex error: " << i << std::endl;
		exit(2); // завершить работу программы, неправильный индекс
	}
	return A[i];
}

// сложение матриц
Matrix Matrix::operator+(const Matrix& right) const 
{
	if (this->GetColumns()!=right.GetColumns() || this->GetLines() != right.GetLines())
	{
		std::cerr << "\nIt is impossible to perform addition. Arrays are different sizes" << std::endl;
		exit(3); // завершить работу программы
	}
	Matrix result(this->GetLines(), this->GetColumns());

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			result[i][j] = A[i][j] + right[i][j];
	return result;
}

// умножение матрицы на число
Matrix Matrix::operator*(double multiplier) const
{
	Matrix result(this->GetLines(), this->GetColumns());
	double num = 0;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			result[i][j] = A[i][j] * multiplier;
	return result;
}

// вычитание матриц
Matrix Matrix::operator-(const Matrix& right) const
{
	if (this->GetColumns() != right.GetColumns() || this->GetLines() != right.GetLines())
	{
		std::cerr << "\nIt is impossible to perform addition. Arrays are different sizes" << std::endl;
		exit(3); // завершить работу программы
	}
	Matrix result(this->GetLines(), this->GetColumns());

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			result[i][j] = A[i][j] - right[i][j];
	return result;
}

// деление матрицы на число
Matrix Matrix::operator/(const double divider) const
{
	Matrix result(this->GetLines(), this->GetColumns());

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			result[i][j] = A[i][j] / divider;
	return result;
}

Matrix Matrix::operator*(const Matrix& right) const
{
	if (this->GetColumns() != right.GetLines())
	{
		std::cerr << "\nIt is impossible to perform multiplication. The number of columns of the first multiplier is not equal to the number of rows of the second multiplier" << std::endl;
		exit(3); // завершить работу программы
	}
	Matrix result(this->GetLines(), right.GetColumns());

	int leftColumns = this->GetColumns();
	int leftLines = this->GetLines();
	int rightColumns = right.GetColumns();

	for (int i = 0; i < leftLines; ++i)
	{
		for (int j = 0; j < rightColumns; ++j)
		{
			for (int k = 0; k < leftColumns; ++k)
				result[i][j] += A[i][k] * right[k][j];
		}
	}
	return result;
}

// оператор присваивания
const Matrix& Matrix::operator= (const Matrix& right)
{
	if (&right != this)
	{
		if (this->GetColumns() != right.GetColumns() || this->GetLines() != right.GetLines())
		{
			A.clear();
			int rightLines = right.GetLines();
			for (size_t i = 0; i < rightLines; ++i)
			{
				A.push_back(right[i]);
			}
		}
		else
		{
			int rightLines = right.GetLines();
			int rightColumns = right.GetColumns();
			for (int i = 0; i < rightLines; ++i)
				for (int j = 0; j < rightColumns; ++j)
					A[i][j] = right[i][j];
		}
	}
	return *this;
}

Matrix Matrix::Transposition() const
{
	Matrix result(this->GetColumns(), this->GetLines());
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			result[j][i] = A[i][j];
	return result;
}

bool Matrix::DiagonalPredominance() const
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			if (A[i][i] < A[i][j]) return false;
	return true;
}

void Matrix::TransferLines(int a, int b)
{
	swap(A[a], A[b]);
}

int Matrix::LeadingElement(int j)
{
	double max = MathOper::Module(this->A[j][j]);
	int i = j;
	for (int counter =j;counter<this->GetLines();++counter)
		if (MathOper::Module(this->A[counter][j]) > max)
		{
			max = MathOper::Module(this->A[counter][j]);
			i = counter;
		}
	return i;
}

Matrix Matrix::MatrixEdition(int n)
{
	if (this->GetLines() == n)
		return *this;
	else
	{
		Matrix result = Matrix::Ones(n);
		int matrixSize = this->GetLines();
		for (int i = n - matrixSize, k = 0; i < n; ++k, ++i)
			for (int j = n - matrixSize, l = 0; j < n; ++l, ++j)
				result[i][j] = A[k][l];
		return result;
	}
}

Matrix Matrix::Delta(const Matrix& right)
{
	if (this->GetColumns() != right.GetColumns() || this->GetLines() != right.GetLines())
	{
		std::cerr << "\nIt is impossible to search Delta. Arrays are different sizes" << std::endl;
		exit(3); // завершить работу программы
	}
	Matrix result(*this - right);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			result[i][j] = MathOper::Module(result[i][j]);
	return result;
}

std::ostream& operator<<(std::ostream& out, const Matrix& x)
{
	int n = x.GetLines();
	int m = x.GetColumns();
	if (n != 1)
	{
		std::cerr << "\nYou can't write a matrix with more than 1 rows, use the 'Write' method" << std::endl;
		exit(4); // завершить работу программы
	}
	std::string matr = "";
	for (int i = 0; i < n; ++i)
	{
		matr += "(";
		for (int j = 0; j < m; ++j)
		{
			matr += std::to_string(x[i][j]);
			if (j != m - 1) matr += ", ";
		}
		matr += ")";
	}
	std::cout << matr;
	return out;
}

void Matrix::Write()
{
	for (int i = 0; i < n; ++i)
	{
		std::cout << "(";
		for (int j = 0; j < m; ++j)
		{
			std::cout << A[i][j];
			if (j != m-1) std::cout << ", ";
		}
		std::cout << ")";
		if (n != 1) std::cout << std::endl;
	}
}

void Matrix::Write() const
{
	for (int i = 0; i < n; ++i)
	{
		std::cout << "(";
		for (int j = 0; j < m; ++j)
		{
			std::cout << A[i][j];
			if (j != m - 1) std::cout << ", ";
		}
		std::cout << ")";
		if (n != 1) std::cout << std::endl;
	}
}