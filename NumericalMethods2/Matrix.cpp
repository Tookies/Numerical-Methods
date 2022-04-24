#include <iostream>
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
	for (int i = 0; i < n; i++)
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
			for (int k = 0; k < leftColumns; k++)
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