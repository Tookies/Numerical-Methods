#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <vector>


class Matrix
{
private:
	std::vector<std::vector<double>> A; // ������ �������� � �������� �������
	int n; // ���������� �����
	int m; // ���������� ��������
public:
	Matrix(); // ����������� �� ���������
	Matrix(int lines, int columns); // ����������� � �����������
	Matrix(const Matrix& matrixCopy); // ����������� �����������
	Matrix(std::vector<std::vector<double>> twoDimMass); // ����������� ��� ���������� �������

	static Matrix Ones(int n);
	static Matrix Ort(int n);

	int GetLines() const { return n; }
	int GetColumns() const { return m; }

	std::vector<double> operator[] (int i) const;
	std::vector<double>& operator[] (int i);
	Matrix operator+(const Matrix& right) const; // �������� ������
	Matrix operator-(const Matrix& right) const; // ��������� ������
	Matrix operator*(const double multiplier) const; // ��������� ������� �� �����
	Matrix operator/(const double divider) const; // ������� ������� �� �����
	Matrix operator*(const Matrix& right) const; // ��������� ������� �� �������
	const Matrix& operator= (const Matrix& right); // �������� ������������

	Matrix Transposition() const; // ���������������� �������
	void TransferLines(int a, int b); // ����������� ����� � �������
	int LeadingElement(int j); // ����������� �������� 
	Matrix MatrixEdition(int n);
};
#endif