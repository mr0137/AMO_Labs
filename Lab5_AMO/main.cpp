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

double* Single_method(double** Data, int N)
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

double Simpson_method_p(int k1, int k2, int N, double H)
{
	double s1 = 0, s2 = 0;
	int i;

	for (i = 1; i < N; i++)
	{
		if (i % 2 == 0)
		{
			s2 += Legandre(k1, A + i * H) * Legandre(k2, A + i * H);
		}
		else
		{
			s1 += Legandre(k1, A + i * H) * Legandre(k2, A + i * H);
		}
	}

	return (H / 3) * (s1 * 4 + s2 * 2 + Legandre(k1, A) * Legandre(k2, A) + Legandre(k1, A + N * H) * Legandre(k2, A + N * H));
}


double Simpson_method_f(int k1, int N, double H)
{
	double s1 = 0, s2 = 0;
	int i;

	for (i = 1; i < N; i++)
	{
		if (i % 2 == 0)
		{
			s2 += Legandre(k1, A + i * H) * f(A + i * H);
		}	
		else
		{
			s1 += Legandre(k1, A + i * H) * f(A + i * H);
		}
	}
		
	return (H / 3) * (s1 * 4 + s2 * 2 + Legandre(k1, A) * f(A) + Legandre(k1, A + N * H) * f(A + N * H));
}

double Integral_p(int k1, int k2)
{
	double In, I2n = 0, eps = 1e-8;
	int n = (int)((B - A) / sqrt(sqrt(eps)));
	double h = double(B - A) / n;

	I2n = Simpson_method_p(k1, k2, n, h);

	do {
		In = I2n;
		n *= 2;
		h /= 2;
		I2n = Simpson_method_p(k1, k2, n, h);
	} while (fabs((In - I2n)/(I2n)) > 15* eps);

	return I2n;
}

double Integral_f(int k1)
{
	double In, I2n = 0, eps = 1e-8;
	int n = (int)((B - A) / sqrt(sqrt(eps)));
	double h = double(B - A) / n;

	I2n = Simpson_method_f(k1, n, h);

	do {
		In = I2n;
		n *= 2;
		h /= 2;
		I2n = Simpson_method_f(k1, n, h);
	} while (fabs((In - I2n)/(I2n * 15)) > 15 *  eps);

	return I2n;
}

double Polinom(double x, int n, double *Data)
{

	double function = 0;
	int i;
	for (i = 0; i < n; i++)
	{
		function += Data[i] * Legandre(i, x);
	}
	return function;
}

double* Data_(int N)
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
			Data[i][j] = Data[j][i] = Integral_p(i, j);
		}
		Data[i][N] = Integral_f(i);
	}

	return Single_method(Data, N);
}
int Polinom(double eps)
{
	
	double x, Error, step = 0.2, *Data;
	int n = 0, k = (int)(B - A) / step, i;
	double h = double(B - A) / k;

	do {
		n++;
		Data = Data_(n);
		Error = 0;
		for (i = 0; i <= k; i++)
		{
			x = Polinom(i * h + A, n, Data) - f(i * h + A);
			Error += x * x;
		}
		Error = sqrt(Error / (k + 1));
		//cout << "Error = " << Error << "N = " << n << endl;
		
	} while (Error > eps);

	return n;
}

int main()
{
	setlocale(LC_ALL, "Rus");
	ofstream  stream;
	
	double eps = 0.1, x, res, *Data;
	int N = 0, i;
	double h;

	cout << "ָה¸ע נאסק¸ע... " << endl;
	N = Polinom(eps);
	system("cls");
	h = (B - A) / N;
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