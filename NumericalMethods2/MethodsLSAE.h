#pragma once
#ifndef METHODSLSAE_H
#define METHODSLSAE_H
#include "MathOper.h"
#include "Matrix.h"

class MethodsLSAE
{	
public:
	static const double eps;
	static Matrix GaussMethod(const Matrix& A, const Matrix& b); // ����� ������
	static Matrix DownInverseGauss(const Matrix& A, const Matrix& b); // �������� ��� ������ ��� ������ ����������� �������
	static Matrix UpInverseGauss(const Matrix& A, const Matrix& b); // �������� ��� ������ ��� ������� ����������� �������
	static Matrix NouseholderMethod(const Matrix& A, const Matrix& b); // ����� �����������
	static Matrix SimpleIterationMethod(const Matrix& A, const Matrix& b); // ����� ������� ��������
	static Matrix SeidelMethod(const Matrix& A, const Matrix& b); // ����� �������
};

#endif