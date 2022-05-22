#include "MethodsLSAE.h"

double MethodsLSAE::Norm(int i, const Matrix& M)
{
	if (i == 0) return MathOper::EuclideanNorm(M);
	else if (i == 1) return MathOper::ColumnNorm(M);
	else if (i == 2) return MathOper::LineNorm(M);
	else if (i == 3) return MathOper::MaxElemNorm(M);
	else
	{
		std::cerr << "\nYou are trying to trigger a norm that is not implemented" << std::endl;
		exit(5);
	}
}

Matrix MethodsLSAE::GaussMethod(const Matrix& A, const Matrix& b)
{
	Matrix M(A);
	Matrix P = Matrix::Ones(A.GetColumns());
	int n = A.GetLines();
	for (int j = 0; j < n - 1; ++j)
	{
		int leadingLine = M.LeadingElement(j);
		M.TransferLines(j, leadingLine);
		P.TransferLines(j, leadingLine);
		for (int i = j + 1; i < n; ++i)
		{
			M[i][j] = M[i][j] / M[j][j];
			for (int k = j + 1; k < n; ++k)
			{
				M[i][k] = M[i][k] - M[i][j] * M[j][k];
			}
		}
	}
	Matrix L = Matrix::Ones(n);
	M = M + L;
	for (int j = 0; j < n - 1; ++j)
	{
		for (int i = j + 1; i < n; ++i)
		{
			L[i][j] = M[i][j];
		}
	}
	M = M - L;
	Matrix y = DownInverseGauss(L, P*b);
	Matrix x = UpInverseGauss(M, y);
	return x;
}

Matrix MethodsLSAE::UpInverseGauss(const Matrix& A, const Matrix& b)
{
	Matrix x(b.GetLines(), b.GetColumns());
	int n = b.GetLines();
	for (int k = n-1; k >= 0; --k) 
	{
		double d = 0;
		for (int j = k + 1; j < n; ++j)
		{
			double s = A[k][j] * x[j][0]; 
			d = d + s; 
		}
		x[k][0] = (b[k][0] - d) / A[k][k]; 
	}
	return x;
}

Matrix MethodsLSAE::DownInverseGauss(const Matrix& A, const Matrix& b)
{
	Matrix x(b.GetLines(), b.GetColumns());
	int n = b.GetLines();
	for (int k = 0; k < n; ++k)
	{
		double d = 0;
		for (int j = k - 1; j >= 0; --j)
		{
			double s = A[k][j] * x[j][0];
			d = d + s;
		}
		x[k][0] = (b[k][0] - d) / A[k][k];
	}
	return x;
}

Matrix MethodsLSAE::NouseholderMethod(const Matrix& A, const Matrix& b)
{
	int n = A.GetLines();
	Matrix Q = Matrix::Ones(n);
	Matrix R(A);
	for (int j = 0; j < n - 1; ++j)
	{
		Matrix z = Matrix::Ort(n - j);
		Matrix y(n - j, 1);
		for (int i = j, k = 0; i < n; ++k, ++i)
		{
			y[k][0] = R[i][j];
		}
		double yNorm = MathOper::EuclideanNorm(y);
		double rhoNorm = MathOper::EuclideanNorm(y-z*yNorm);
		Matrix w = (y - z * yNorm) / rhoNorm;
		Matrix E = Matrix::Ones(w.GetLines());
		Matrix Qi = (E - ((w * w.Transposition()) * 2)).MatrixEdition(R.GetLines());
		Q = Q * Qi;
		R = Qi * R;
	}
	Matrix y = Q.Transposition() * b;
	Matrix x = UpInverseGauss(R, y);
	return x;
}

Matrix MethodsLSAE::SimpleIterationMethod(const Matrix& An, const Matrix& bn, int& counter, double eps)
{
	Matrix A(An);
	double mu = 0;
	int n = A.GetLines();
	Matrix B(n, n);
	double BNorm = 1;

	Matrix b(bn);
	bool transp = true;
	for (int i = 1; i < 3; ++i)
	{
		mu = 1 / Norm(i, A);
		B = Matrix::Ones(n) - (A * mu);
		BNorm = Norm(i, B);
		if (BNorm < 1) break;
		if (i == 2 && BNorm >= 1 && transp)
		{
			if (transp)
			{
				transp = false;
				i = 1;
			}
			b = An.Transposition() * b;
			A = An.Transposition() * A;
		}
	}

	double exprBNorm = BNorm/(1-BNorm);
	Matrix c = b * mu;
	Matrix xk(n,1);
	Matrix xk1(c);
	do
	{
		counter++;
		xk = xk1;
		xk1 = B * xk + c;
	} while (MathOper::EuclideanNorm(A * xk1 - b) > eps);
	return xk1;
}

Matrix MethodsLSAE::SeidelMethod(const Matrix& An, const Matrix& bn, int& counter, double eps)
{
	Matrix A(An);
	Matrix b(bn);
	if (!A.DiagonalPredominance())
	{
		A = An.Transposition() * A;
		b = An.Transposition() * b;
	}
	int n = A.GetLines();
	Matrix c = Matrix(n, n);
	Matrix d = Matrix(n, 1);
	for (int i = 0; i < n; ++i)
	{
		d[i][0] = b[i][0] / A[i][i];
		for (int j = 0; j < n; ++j)
		{
			if (i != j) c[i][j] = -A[i][j] / A[i][i];
		}
	}
	Matrix xk = d;
	Matrix xk1 = d;
	do
	{
		counter++;
		for (int i = 0; i < n; ++i)
		{
			xk1[i][0] = 0;
			for (int j = 0; j < n; ++j)
			{
				if (j < i)
				{
					xk1[i][0] += c[i][j] * xk1[j][0];
				}
				else if (j > i)
				{
					xk1[i][0] += c[i][j] * xk[j][0];
				}
			}
			xk1[i][0] += d[i][0];
		}
		xk = xk1;
	} while (MathOper::EuclideanNorm(A*xk1-b)>eps);
	return xk1;
}