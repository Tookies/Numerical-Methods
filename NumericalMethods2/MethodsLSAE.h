#pragma once
#ifndef METHODSLSAE_H
#define METHODSLSAE_H
#include "MathOper.h"
#include "Matrix.h"

class MethodsLSAE
{	
public:
	const double eps = 0.00001;
	Matrix GaussMethod(const Matrix& A, const Matrix& b); // ����� ������
	Matrix DownInverseGauss(const Matrix& A, const Matrix& b); // �������� ��� ������ ��� ������ ����������� �������
	Matrix UpInverseGauss(const Matrix& A, const Matrix& b); // �������� ��� ������ ��� ������� ����������� �������
	Matrix NouseholderMethod(const Matrix& A, const Matrix& b); // ����� �����������
	Matrix SimpleIterationMethod(const Matrix& A, const Matrix& b); // ����� ������� ��������
	Matrix SeidelMethod(const Matrix& A, const Matrix& b); // ����� �������
};
#endif