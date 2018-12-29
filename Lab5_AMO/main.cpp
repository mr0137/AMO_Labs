#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#define A 0.0
#define B 7.0

double f(double x) 
{
	return ((x*x - sin(x))*(cos(2.5 * x)));
}

double* Single_method(double **Data, int N) 
{
	int i, j, k, m;
	double temp1, temp2, temp3;
	double* res = new double[N];

	for (i = 0; i < N; i++) 
	{
		temp1 = Data[i][i];
		for (j = i; j < N + 1; j++) 
		{
			Data[i][j] /= temp1;
		}
		for (k = i + 1; k < N; k++)
		{
			temp2 = Data[k][i];
			for (m = i; m < N + 1; m++) 
			{
				Data[k][m] -= Data[i][m] * temp2;
			}
		}
	}
	for (i = N - 1; i >= 0; i--) 
	{
		temp3 = 0;
		for (j = N - 1; j > i; j--) 
		{
			temp3 += Data[i][j] * res[j];
		}
		res[i] = Data[i][N] - temp3;
	}
	return res;
}

double Legandre(int n, double x) 
{
	double Ln1, Ln = x, Ln_1 = 1;
	if (n == 0) return Ln_1;
	if (n == 1) return Ln;

	int i = 1;
	while (i < n) 
	{
		Ln1 = (1.0 * (2 * i + 1) / (i + 1)) * x * Ln - (1.0 * i / (i + 1)) * Ln_1;
		Ln_1 = Ln;
		Ln = Ln1;
		++i;
	}
	return Ln;
}

double func(int N, double x, int k1, int k2)
{
	double t = (2 * x - A - B) / (B - A);
	if (k2 == N)
	{
		return  f(x)*Legandre(k1, t);
	}
	return  Legandre(k1, t)*Legandre(k2, t);
}

double Simpson_method(int k1, int k2, int n, int N) 
{
	const double h = (B - A) / n;
	double s1 = 0, s2 = 0;
	int i;
	for (i = 1; i < n; i += 2) 
	{
		s1 += func(N, A + i*h, k1, k2);
		s2 += func(N, A + (i + 1)*h, k1, k2);
	}
	return h / 3 * (f(A) + 4 * s1 + 2 * s2);
}


double Integral(int k1, int k2, int N) 
{
	double In, I2n, eps = 1e-8;
	int n = (int)((B - A) / sqrt(sqrt(eps)));

	In = Simpson_method(k1, k2, n, N);
	I2n = Simpson_method(k1, k2, 2 * n, N);

	while (fabs(In - I2n) / 15 > eps) 
	{
		In = I2n;
		n *= 2;
		I2n = Simpson_method(k1, k2, n, N);
	}
	return I2n;
}

double *Data_(int N) 
{
	double **Data = new double*[N];
	int i, j;
	for (i = 0; i < N + 1; i++)
	{
		Data[i] = new double[N + 1];
	}

	for (i = 0; i < N; i++) 
	{
		for (j = 0; j < N + 1; j++) 
		{
			Data[i][j] = Data[j][i] = Integral(i, j, N);
		}
		Data[i][N] = Integral(i, N, N);
	}
	
	return Single_method(Data, N);;
}

double Polinom(double x, int n, double *Data)
{
	double  t = (2 * x - B - A) / (B - A), function = 0;
	int i;
	for (i = 0; i < n; i++)
	{
		function += Data[i] * Legandre(i, t);
	}
	return function;
}

int Polinom(double eps) 
{
	double y, x, Error, step = 0.2, *Data;
	int n = 2, k = (int)(B - A) / step, i;

	do {
		Error = 0;
		n *= 2;
		x = A;
		for (i = 0; i <= k; ++i) 
		{
			Data = Data_(n);
			y = f(x) - Polinom(x, n, Data);
			Error += y * y;
			x += step;
		}
		Error = sqrt(Error / (k + 1));
	} while (Error > eps);

	return n;
}

int main()
{
	setlocale(LC_ALL, "Rus");
	ofstream  stream;
	double eps = 0.1, x, res, *Data;
	int N = 0;

	cout << "ָה¸ע נאסק¸ע... " << endl;
	N = Polinom(eps);

	system("cls");
	cout << "\t\tN = " << N << endl;
	stream.open("Data.csv");
	cout << "|==========|=====================|" << endl;
	for (x = A; x <= B; x += 0.2) 
	{
		Data = Data_(N);
		res = Polinom(x, N, Data);
		stream << fixed << x << ";" << res << endl;
		cout << "|" << setw(10) << setprecision(3) << x << "|" << setw(21) << setprecision(10) << res << "|" << endl;
	}
	cout << "|==========|=====================|" << endl;
	stream.close();
	system("pause");
	return 0;
}