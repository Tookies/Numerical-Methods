#include "MethodsLSAE.h"

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
	M = M - L;
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
		for (int i = j, int k = 0; i < n; ++k, ++i)
		{
			y[k][0] = R[i][j];
		}
		double yNorm = MathOper::EuclideanNorm(y);
		double rhoNorm = MathOper::EuclideanNorm(y-z*yNorm);
		Matrix w = (y - z * yNorm) / rhoNorm;
		Matrix E = Matrix::Ones(w.GetLines());
		Matrix Qi = (E - ((w * w.Transposition()) * 2)).MatrixEdition(R.GetLines());
		if (j == 0) Q = Qi;
		R = Q * R;
	}
	Matrix x = UpInverseGauss(R, Q.Transposition() * b);
	return x;
}

Matrix MethodsLSAE::SimpleIterationMethod(const Matrix& A, const Matrix& b)
{
	int n = A.GetLines();
	double mu = 1 / MathOper::EuclideanNorm(A);
	Matrix c = b * mu;
	Matrix B = Matrix::Ones(n) - (A * mu);
	double BNorm = MathOper::EuclideanNorm(B);
	double exprBNorm = BNorm/(1-BNorm);
	Matrix xk(n,1);
	Matrix xk1(c);
	do
	{
		xk = xk1;
		xk1 = B * xk + c;
	} while (exprBNorm * (MathOper::EuclideanNorm(xk1 - xk)) > eps);
	return xk1;
}

Matrix SeidelMethod(const Matrix& A, const Matrix& b)