#include <iostream>
#include "MathOper.h"
#include "Matrix.h"


int main()
{
	std::vector<std::vector<double> > vect1;
	vect1.push_back({ 1,2,3 });
	vect1.push_back({ 1,0,-1 });
	std::vector<std::vector<double> > vect2;
	vect2.push_back({ 3,4,5 });
	vect2.push_back({ 6,0,-2 });
	vect2.push_back({ 7,1,8 });
	vect2.push_back({ 1,0,-3 });
	Matrix A(vect1);
	Matrix B(vect2);
	double euclideanNorm = MathOper::EuclideanNorm(B);
	double columnNorm = MathOper::ColumnNorm(B);
	double lineNorm = MathOper::LineNorm(B);
	
    std::cout << "Ok";
}
