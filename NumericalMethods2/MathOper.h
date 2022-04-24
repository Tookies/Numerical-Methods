#pragma once
#ifndef MATHOPER_H
#define MATHOPER_H
#include "Matrix.h"

class MathOper
{
public:
	static double Module(double numb); // ������� ������
	static double SqrtHer(double numb); // ������� ����������� �����
	static double EuclideanNorm(const Matrix& M); // ��������� �����
	static double ColumnNorm(const Matrix& M); // ��������� ����� ���������� ��������� ����� ||x|| = SUM{i in 1,n}(|x[i]|)
	static double LineNorm(const Matrix& M); // ��������� ����� ���������� ��������� ����� ||x|| = MAX{i in 1,n}(|x[i]|)
};
#endif