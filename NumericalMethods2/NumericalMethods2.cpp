#include <iostream>
#include <iomanip>
#include <sstream>
#include "MathOper.h"
#include "Matrix.h"
#include "MethodsLSAE.h"

void PrintTable1(int testNumber, const Matrix& A, const Matrix b, const Matrix& answer)
{
	std::cout << std::endl << "Number of test: " << testNumber << std::endl;
	std::cout << std::endl << "The exact solution: " << answer.Transposition() << std::endl;

	std::cout << " ______________________________________________________________________________________________________________________________________________________" << std::endl;
	std::cout << "|" << std::setw(12) << " " << "|" << std::setw(68) << "Simple iteration method" << "|" << std::setw(68) << "The Seidel method" << "|" << std::endl;
	std::cout << "|" << std::setw(12) << "eps" <<  "|____________________________________________________________________|____________________________________________________________________|" << std::endl;
	std::cout << "|" << std::setw(12) << " " << "|" << std::setw(48) << "x" << "|" << std::setw(12) << "delta" << "|" << std::setw(6) << "k" << "|" << std::setw(48) << "x" << "|" << std::setw(12) << "delta" << "|" << std::setw(6) << "k" << "|" << std::endl;
	std::cout << " ______________________________________________________________________________________________________________________________________________________" << std::endl;
	double eps = 0.1;
	for (int x = 2; x <= 6; ++x)
	{
		eps /= 10;
		int counter1 = 0;
		Matrix xSIM = MethodsLSAE::SimpleIterationMethod(A, b, counter1, eps);
		int counter2 = 0;
		Matrix xSM = MethodsLSAE::SeidelMethod(A, b, counter2, eps);
		std::cout << " ______________________________________________________________________________________________________________________________________________________" << std::endl;
		std::cout << "|" << std::setw(12) << eps << "|"  << std::setw(48) << xSIM.Transposition() << "|" << std::setw(12) << MathOper::EuclideanNorm(answer - xSIM) << "|" << std::setw(6) << counter1 << "|" << std::setw(48) << xSM.Transposition() << "|" << std::setw(12) << MathOper::EuclideanNorm(answer - xSM) << "|" << std::setw(6) << counter2<< "|" << std::endl;
	}
	std::cout << std::endl << std::endl;
	std::cout << " ___________________________________________________________________________________________________________________________" << std::endl;
	std::cout << "|" << std::setw(61) << "The Gauss method" << "|" << std::setw(61) << "The Nouseholder Method" << "|" << std::endl;
	std::cout << "|_____________________________________________________________|_____________________________________________________________|" << std::endl;
	std::cout << "|" << std::setw(48) << "x" << "|" << std::setw(12) << "delta" << "|" << std::setw(48) << "x" << "|" << std::setw(12) << "delta" << "|" << std::endl;

	Matrix xG = MethodsLSAE::GaussMethod(A, b);
	Matrix xH = MethodsLSAE::NouseholderMethod(A, b);
	std::cout << " ___________________________________________________________________________________________________________________________" << std::endl;
	std::cout << "|" << std::setw(48) << xG.Transposition() << "|" << std::setw(12) << MathOper::EuclideanNorm(answer - xG) << "|" << std::setw(48) << xH.Transposition() << "|" << std::setw(12) << MathOper::EuclideanNorm(answer - xH) << "|" << std::endl;
}

std::string _String(int n)
{
	std::stringstream ss;
	for (int i = 0; i <= n; ++i)
	{
		ss << "_";
	}
	return ss.str();
}

void PrintTable2(int N)
{
	std::cout << std::endl << "Number of test: " << 5 << std::endl;
	for (int i = 4; i <= 6; ++i)
	{
		double e = 1;
		for (int j = 0; j<2; ++j)
		{
			e /= 1000;
			Matrix A = Matrix::FormAForTest5(i, N, e);
			Matrix b = Matrix::FormbForTest5(i);
			Matrix xG = MethodsLSAE::GaussMethod(A, b);
			Matrix xH = MethodsLSAE::NouseholderMethod(A, b);
			std::cout << std::endl << "The exact solution: " << xG.Transposition() << std::endl;
			std::cout << std::endl << "e value: " << e << std::endl;

			std::cout << " _________________________" << _String(10 * i + 8 + 12 + 6 + 2)<< _String(10 * i + 8 + 12 + 6 + 2) << std::endl;
			std::cout << "|" << std::setw(12) << " " << "|" << std::setw(12) << " " << "|" << std::setw(10 * i + 8 + 12 + 6 + 2) << "Simple iteration method" << "|" << std::setw(10 * i + 8 + 12 + 6 + 2) << "The Seidel method" << "|" << std::endl;
			std::cout << "|" << std::setw(12) << "e" << "|" << std::setw(12) << "eps" << "|"<< _String(10 * i + 8 + 12 + 6+1)<<"|"<< _String(10 * i + 8 + 12 + 6+1)<< "|" << std::endl;
			std::cout << "|" << std::setw(12) << " " << "|" << std::setw(12) << " " << "|" << std::setw(10 * i + 8) << "x" << "|" << std::setw(12) << "delta" << "|" << std::setw(6) << "k" << "|" << std::setw(10 * i + 8) << "x" << "|" << std::setw(12) << "delta" << "|" << std::setw(6) << "k" << "|" << std::endl;
			std::cout << " _________________________" << _String(10 * i + 8 + 12 + 6 + 2) << _String(10 * i + 8 + 12 + 6 + 2) << std::endl;
			double eps = 0.1;
			for (int x = 2; x <= 6; ++x)
			{
				eps /= 10;
				int counter1 = 0;
				Matrix xSIM = MethodsLSAE::SimpleIterationMethod(A, b, counter1, eps);
				int counter2 = 0;
				Matrix xSM = MethodsLSAE::SeidelMethod(A, b, counter2, eps);
				std::cout << " _________________________" << _String(10 * i + 8 + 12 + 6 + 2) << _String(10 * i + 8 + 12 + 6 + 2) << std::endl;
				std::cout << "|" << std::setw(12) << e << "|" << std::setw(12) << eps << "|" << std::setw(10 * i + 8) << xSIM.Transposition() << "|" << std::setw(12) << MathOper::EuclideanNorm(xG - xSIM) << "|" << std::setw(6) << counter1 << "|" << std::setw(10 * i + 8) << xSM.Transposition() << "|" << std::setw(12) << MathOper::EuclideanNorm(xG - xSM) << "|" << std::setw(6) << counter2 << "|" << std::endl;
			}
			std::cout << std::endl << std::endl;
			std::cout << " ______________" << _String(10 * i + 8 + 12) << _String(10 * i + 8 + 12) << std::endl;
			std::cout << "|" << std::setw(12) << " " << "|" << std::setw(10 * i + 8 + 12 + 1) << "The Gauss method" << "|" << std::setw(10 * i + 8 + 12 + 1) << "The Nouseholder Method" << "|" << std::endl;
			std::cout << "|" << std::setw(12) << "e" << "|"<< _String(10 * i + 8 + 12) <<"|"<< _String(10 * i + 8 + 12) <<"|" << std::endl;
			std::cout << "|" << std::setw(12) << " " << "|" << std::setw(10 * i + 8) << "x" << "|" << std::setw(12) << "delta" << "|" << std::setw(10 * i + 8) << "x" << "|" << std::setw(12) << "delta" << "|" << std::endl;


			std::cout << " ______________" << _String(10 * i + 8 + 12) << _String(10 * i + 8 + 12) << std::endl;
			std::cout << "|" << std::setw(12) << e << "|" << std::setw(10 * i + 8) << xG.Transposition() << "|" << std::setw(12) << MathOper::EuclideanNorm(xG - xG) << "|" << std::setw(10 * i + 8) << xH.Transposition() << "|" << std::setw(12) << MathOper::EuclideanNorm(xG - xH) << "|" << std::endl;
		}
	}
}

int main()
{
	/*std::cout << std::fixed;
	std::cout.precision(4);*/
	std::vector<std::vector<double> > vect;
	vect.push_back({ 1, 2, 3 });
	const Matrix test0Answer((Matrix (vect)).Transposition());
	vect.pop_back();
	vect.push_back({ 1, 1, 1 });
	const Matrix test1Answer((Matrix(vect)).Transposition());
	vect.pop_back();
	vect.push_back({ 1403.0/1139, 1379.0/1139, 1359.0/1139 });
	const Matrix test2Answer((Matrix(vect)).Transposition());
	vect.pop_back();
	vect.push_back({ 21529.0 / 17726, 9643.0 / 8863, 18837.0 / 17726 });
	const Matrix test3Answer((Matrix(vect)).Transposition());
	vect.pop_back();
	vect.push_back({ -27.0 / 143, 259.0 / 429, 109.0 / 143 });
	const Matrix test4Answer((Matrix(vect)).Transposition());

	std::vector<std::vector<double> > vect1;
	vect1.push_back({ 0, 2, 3 });
	vect1.push_back({ 1, 2, 4 });
	vect1.push_back({ 4, 5, 6 });
	std::vector<std::vector<double> > vect2;
	vect2.push_back({ 13 });
	vect2.push_back({ 17 });
	vect2.push_back({ 32 });
	const Matrix A0(vect1);
	const Matrix b0(vect2);
	double N = 17;
	vect1.erase(vect1.begin(), vect1.end());
	vect2.erase(vect2.begin(), vect2.end());
	vect1.push_back({ N+2, 1, 1 });
	vect1.push_back({ 1, N+4, 1 });
	vect1.push_back({ 1, 1, N+6 });
	vect2.push_back({ N+4 });
	vect2.push_back({ N+6 });
	vect2.push_back({ N+8 });
	const Matrix A1(vect1);
	const Matrix b1(vect2);
	vect1.erase(vect1.begin(), vect1.end());
	vect2.erase(vect2.begin(), vect2.end());
	vect1.push_back({ -(N+2), 1, 1 });
	vect1.push_back({ 1, -(N+4), 1 });
	vect1.push_back({ 1, 1, -(N+6) });
	vect2.push_back({ -(N+4) });
	vect2.push_back({ -(N+6) });
	vect2.push_back({ -(N+8) });
	const Matrix A2(vect1);
	const Matrix b2(vect2);
	vect1.erase(vect1.begin(), vect1.end());
	vect2.erase(vect2.begin(), vect2.end());
	vect1.push_back({ -(N + 2), N+3, N+4 });
	vect1.push_back({ N+5, -(N + 4), N+1 });
	vect1.push_back({ N+4, N+5, -(N + 6) });
	vect2.push_back({ N + 4 });
	vect2.push_back({ N + 6 });
	vect2.push_back({ N + 8 });
	const Matrix A3(vect1);
	const Matrix b3(vect2);
	vect1.erase(vect1.begin(), vect1.end());
	vect2.erase(vect2.begin(), vect2.end());
	vect1.push_back({ N + 2, N + 1, N + 1 });
	vect1.push_back({ N + 1, N + 4, N + 1 });
	vect1.push_back({ N + 1, N + 1, N + 6 });
	vect2.push_back({ N + 4 });
	vect2.push_back({ N + 6 });
	vect2.push_back({ N + 8 });
	const Matrix A4(vect1);
	const Matrix b4(vect2);

	PrintTable1(0, A0, b0, test0Answer);
	PrintTable1(1, A1, b1, test1Answer);
	PrintTable1(2, A2, b2, test2Answer);
	PrintTable1(3, A3, b3, test3Answer);
	PrintTable1(4, A4, b4, test4Answer);

	PrintTable2(17);
}
