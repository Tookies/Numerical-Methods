#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <iostream>
#include <string>

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
	static Matrix FormAForTest5(int n, int N, double e);
	static Matrix FormbForTest5(int n);

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
	bool DiagonalPredominance() const; // ������� ����� �� ������� ���������
	void TransferLines(int a, int b); // ����������� ����� � �������
	int LeadingElement(int j); // ����������� �������� 
	Matrix MatrixEdition(int n);
	Matrix Delta(const Matrix& right);
	void Write();
	void Write() const;
	friend std::ostream& operator<< (std::ostream& out, const Matrix& x);
};
#endif