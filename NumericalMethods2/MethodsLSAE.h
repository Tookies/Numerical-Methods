#pragma once
#ifndef METHODSLSAE_H
#define METHODSLSAE_H
#include "MathOper.h"
#include "Matrix.h"

class MethodsLSAE
{	
public:
	static Matrix GaussMethod(const Matrix& A, const Matrix& b); // метод Гаусса
	static Matrix DownInverseGauss(const Matrix& A, const Matrix& b); // обратный ход Гаусса для нижней треугольной матрицы
	static Matrix UpInverseGauss(const Matrix& A, const Matrix& b); // обратный ход Гаусса для верхней треугольной матрицы
	static Matrix NouseholderMethod(const Matrix& A, const Matrix& b); // метод Хаусхолдера
	static Matrix SimpleIterationMethod(const Matrix& A, const Matrix& b, int& counter, double eps); // метод простой итерации
	static Matrix SeidelMethod(const Matrix& An, const Matrix& b, int& counter, double eps); // метод Зейделя

	static double Norm(int i ,const Matrix& M);
};

#endif