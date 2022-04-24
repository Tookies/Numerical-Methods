#include <iostream>
#include "MathOper.h"
#include "Matrix.h"
#include "MethodsLSAE.h"


int main()
{
	std::vector<std::vector<double> > vect1;
	vect1.push_back({ 3, 1, 1 });
	vect1.push_back({ 1, 5, 1 });
	vect1.push_back({ 1, 1, 7 });
	std::vector<std::vector<double> > vect2;
	vect2.push_back({ 5 });
	vect2.push_back({ 7 });
	vect2.push_back({ 9 });
	Matrix A(vect1);
	Matrix b(vect2);
	Matrix x = MethodsLSAE::SeidelMethod(A, b);
	
    std::cout << "\u2552\u2550\u2550\u2564\u2550\u2550\u2555";
}
