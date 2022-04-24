#pragma once
#ifndef MATHOPER_H
#define MATHOPER_H
#include "Matrix.h"

class MathOper
{
public:
	static double Module(double numb); // функция модуля
	static double SqrtHer(double numb); // функция квадратного корня
	static double EuclideanNorm(const Matrix& M); // евклидова норма
	static double ColumnNorm(const Matrix& M); // матричная норма подчинённая векторной норме ||x|| = SUM{i in 1,n}(|x[i]|)
	static double LineNorm(const Matrix& M); // матричная норма подчинённая векторной норме ||x|| = MAX{i in 1,n}(|x[i]|)
};
#endif