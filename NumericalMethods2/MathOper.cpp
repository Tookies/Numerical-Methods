#include <iostream>
#include "MathOper.h"


// функция модуля
double MathOper::Module(double numb) 
{
	return (numb >= 0) ? numb: -numb;
}

// функция квадратного корня
double MathOper::SqrtHer(double numb) 
{
	if (numb < 0)
	{
		std::cerr << "\nError root of a negative number: " << numb << std::endl;
		exit(1); // завершить работу программы, корень из отрицательного числа
	}
	double n = 1;
	double e = 0.00000000001;
	double sqrt = (1.0 / 2) * (n + numb / n);
	while (Module(sqrt - n) > e)
	{
		n = sqrt;
		sqrt = (1.0 / 2) * (n + numb / n);
	}
	return sqrt;
}

// евкидова норма
double MathOper::EuclideanNorm(const Matrix& M)
{
	double result = 0;
	int MColumns = M.GetColumns();
	int MLines = M.GetLines();
	for (int i = 0; i < MLines; ++i)
		for (int j = 0; j < MColumns; ++j)
			result += M[i][j] * M[i][j];
	return MathOper::SqrtHer(result);
}

// матричная норма подчинённая векторной норме ||x|| = SUM{i in 1,n}(|x[i]|)
double MathOper::ColumnNorm(const Matrix& M)
{
	double result = -1;
	double sum = 0;
	int MColumns = M.GetColumns();
	int MLines = M.GetLines();
	for (int j = 0; j < MColumns; ++j)
	{
		for (int i = 0; i < MLines; ++i)
		{
			sum += MathOper::Module(M[i][j]);
		}
		if (sum > result) result = sum;
		sum = 0;
	}
	return result;
}

// матричная норма подчинённая векторной норме ||x|| = MAX{i in 1,n}(|x[i]|)
double MathOper::LineNorm(const Matrix& M)
{
	double result = -1;
	double sum = 0;
	int MColumns = M.GetColumns();
	int MLines = M.GetLines();
	for (int i = 0; i < MLines; ++i)
	{
		for (int j = 0; j < MColumns; ++j)
		{

			sum += MathOper::Module(M[i][j]);
		}
		if (sum > result) result = sum;
		sum = 0;
	}
	return result;
}