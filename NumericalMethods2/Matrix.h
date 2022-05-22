#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <iostream>
#include <string>

class Matrix
{
private:
	std::vector<std::vector<double>> A; // вектор векторов в качестве матрицы
	int n; // количество строк
	int m; // количество столбцов
public:
	Matrix(); // конструктор по умолчанию
	Matrix(int lines, int columns); // конструктор с параметрами
	Matrix(const Matrix& matrixCopy); // конструктор копирования
	Matrix(std::vector<std::vector<double>> twoDimMass); // конструктор для двумерного массива

	static Matrix Ones(int n);
	static Matrix Ort(int n);
	static Matrix FormAForTest5(int n, int N, double e);
	static Matrix FormbForTest5(int n);

	int GetLines() const { return n; }
	int GetColumns() const { return m; }	

	std::vector<double> operator[] (int i) const;
	std::vector<double>& operator[] (int i);
	Matrix operator+(const Matrix& right) const; // сложение матриц
	Matrix operator-(const Matrix& right) const; // вычитание матриц
	Matrix operator*(const double multiplier) const; // умножение матрицы на число
	Matrix operator/(const double divider) const; // деление матрицы на число
	Matrix operator*(const Matrix& right) const; // умножение матрицы на матрицу
	const Matrix& operator= (const Matrix& right); // оператор присваивания

	Matrix Transposition() const; // транспонирование матрицы
	bool DiagonalPredominance() const; // наличие нулей на главной диагонале
	void TransferLines(int a, int b); // перемещение строк в матрице
	int LeadingElement(int j); // определение ведущего 
	Matrix MatrixEdition(int n);
	Matrix Delta(const Matrix& right);
	void Write();
	void Write() const;
	friend std::ostream& operator<< (std::ostream& out, const Matrix& x);
};
#endif