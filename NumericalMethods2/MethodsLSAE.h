#pragma once
#ifndef METHODSLSAE_H
#define METHODSLSAE_H
#include "MathOper.h"
#include "Matrix.h"

class MethodsLSAE
{	
public:
	const double eps = 0.00001;
	Matrix GaussMethod(const Matrix& A, const Matrix& b); // метод Гаусса
	Matrix DownInverseGauss(const Matrix& A, const Matrix& b); // обратный ход Гаусса для нижней треугольной матрицы
	Matrix UpInverseGauss(const Matrix& A, const Matrix& b); // обратный ход Гаусса для верхней треугольной матрицы
	Matrix NouseholderMethod(const Matrix& A, const Matrix& b); // метод Хаусхолдера
	Matrix SimpleIterationMethod(const Matrix& A, const Matrix& b); // метод простой итерации
	Matrix SeidelMethod(const Matrix& A, const Matrix& b); // метод Зейделя
};
#endif