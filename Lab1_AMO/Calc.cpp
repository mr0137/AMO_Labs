#include "Calc.h"

void getCos_1(double x, double eps, double &R, int &n, double &U, double &result) 
{
	do
	{
		n++;
		result += U;
		U = -(pow(x, 2) / (2 * n * (2 * n - 1))) * U;
	} while (abs(U) > eps);
	R = U;
}

void First(double x, int &n_forTask2)
{	
	cout << "|=======|=====|=====================|=====================|" << endl;
	cout << "|" << setw(7) << setprecision(5) << "eps" << "|" << setw(5) << setprecision(4) << "n" << "|" << setw(21) << "fluff" << "|" << setw(21) << "R" << "|" << endl;
	cout << "|=======|=====|=====================|=====================|" << endl;
	double R, result, eps, U, fluff;
	int n;
	x = abs(x);
	while (x >= 2 * PI)
	{
		x = x - 2 * PI;
	}

	for (eps = pow(10, -2); eps >= pow(10, -14); eps *= pow(10, -3)) 
	{
		result = n = 0.0;
		U = 1.0;
		getCos_1(x, eps, R, n, U, result);
		fluff = abs(cos(x) - result);
		cout <<"|"<<setw(7) << setprecision(5) << eps <<"|"<< setw(5) << setprecision(4) << n <<"|"<< setw(21) << setprecision(14) << fluff <<"|"<< setw(21) << R <<"|"<< endl;
		if (eps == pow(10, -8))
		{
			n_forTask2 = n;
		}
	}
	cout << endl;
}

void getCos_2(double x, double &R, int n, double &result, double &U) 
{
	int i;
	U = 1;

	for (i = 1; i <= n; i++) 
	{
		result += U;
		U = -(pow(x, 2) / (2 * i * (2 * i - 1))) * U;
	}
	R = U;
}

void Second(int n, double a, double h) 
{
	cout << "|=======|===========================|=====================|" << endl;
	cout <<"|"<< setw(7) << setprecision(4) <<"X" <<"|"<< setw(27) << setprecision(14) << "fluff" <<"|"<< setw(21) << "R" <<"|"<< endl;
	cout << "|=======|===========================|=====================|" << endl;
	int i;
	double x, result, U, R;
	double fluff;
	for (i = 0; i <= 10; i++) 
	{
		x = a + h * i;
		x = abs(x);

		while (x >= 2 * PI)
		{
			x = x - 2 * PI;
		}

		result = 0.0;
		U = 1.0;
		getCos_2(x, R, n, result, U);

		fluff = abs(cos(x) - result);
		R = U;
		cout << "|" << setw(7) << setprecision(4) << a + h * i << "|" << setw(27) << setprecision(14) << fluff << "|" << setw(21) << R << "|" << endl;
	}
	cout << "|=======|===========================|=====================|" << endl;
}