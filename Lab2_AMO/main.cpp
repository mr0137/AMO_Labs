#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>

using namespace std;

double function(double x)
{
	return log(2 + sin(x)) - (x / 10) - 0.5;
}

double Dfunction(double x)
{
	return (cos(x) / (2 + sin(x)) - 0.1);
}

void iteration(double a, double b)
{
	double m, M;
	double x, eps, x1;
	double lamda;
	double q;

	cout << "|=======|===============|=====================|" << endl;
	cout << "|" << setw(7) << setprecision(5) << "eps" << "|" << setw(15) << setprecision(10) << "x" << "|" << setw(21) << "fluff" << "|" << endl;
	cout << "|=======|===============|=====================|" << endl;

	if (fabs(Dfunction(a)) > fabs(Dfunction(b)))
	{
		M = Dfunction(a);
		m = Dfunction(b);
	}
	else
	{
		m = Dfunction(a);
		M = Dfunction(b);
	}
	q = 1 - fabs(m / M);
	lamda = 1 / M;
	for (eps = 0.01; eps >= 1e-14; eps *= 1e-3)
	{
		x = (a + b) / 2;
		do
		{
			x1 = x;
			x = x1 - lamda * function(x1);
		} while (fabs(x1 - x) > ((1 - q) / q * eps));
		cout << "|" << setw(7) << setprecision(5) << eps << "|" << setw(15) << setprecision(10) << x << "|" << setw(21) << setprecision(10) << fabs((x1 - x) * q / (1 - q)) << "|" << endl;

	}
	cout << "|=======|===============|=====================|" << endl;
}
void bisection(double a, double b)
{
	double x, eps;
	cout << "|=======|===============|=====================|" << endl;
	cout << "|" << setw(7) << setprecision(5) << "eps" << "|" << setw(15) << setprecision(10) << "x" << "|" << setw(21) << "fluff" << "|" << endl;
	cout << "|=======|===============|=====================|" << endl;
	for (eps = 0.01; eps >= 1e-14; eps *= 1e-3)
	{
		while (b - a >= 2 * eps)
		{
			x = (a + b) / 2;
			if (function(x) * function(a) < 0)
			{
				b = x;
			}
			else if (function(x) * function(b) < 0)
			{
				a = x;
			}
		}
		x = (a + b) / 2;
		cout << "|" << setw(7) << setprecision(5) << eps << "|" << setw(15) << setprecision(10) << x << "|" << setw(21) << setprecision(10) << fabs((b - a)/2) << "|" << endl;
	}
	cout << "|=======|===============|=====================|" << endl;
}
void comparison(double a, double b, double counterIT, double counterB)
{
	double x, eps, x1;
	double m, M;
	double lamda;
	double q;
	cout << "|=======|===============|==========|==========|" << endl;
	cout << "|" << setw(7) << setprecision(5) << "eps" << "|" << setw(15) << setprecision(10) << "x" << "|" << setw(10) << "IT" << "|" << setw(10) << "B" << "|" << endl;
	cout << "|=======|===============|==========|==========|" << endl;
	if (fabs(Dfunction(a)) > fabs(Dfunction(b)))
	{
		M = Dfunction(a);
		m = Dfunction(b);
	}
	else
	{
		m = Dfunction(a);
		M = Dfunction(b);
	}
	q = 1 - fabs(m / M);
	lamda = 1 / M;
	for (eps = 0.01; eps >= 1e-14; eps *= 1e-3)
	{
		counterIT = 0;
		counterB = 0;
		/*=============================*/
		x1 = 0;
		x = (a + b) / 2;
		do
		{
			x1 = x;
			x = x1 - lamda * function(x1);
			counterIT++;
		} while (fabs(x1 - x) > ((1 - q) / q * eps));

		/*=============================*/
	
		while (b - a >= 2 * eps)
		{
			x = (a + b) / 2;
			if (function(x) * function(a) < 0)
			{
				b = x;
			}
			else if (function(x) * function(b) < 0)
			{
				a = x;
			}
			counterB++;
		/*=============================*/
		}
		cout << "|" << setw(7) << setprecision(5) << eps << "|" << setw(15) << setprecision(10) << "x" << "|" << setw(10) << counterIT << "|" << setw(10) << counterB << "|" << endl;
		cout << "|=======|===============|==========|==========|" << endl;
	}
}

int main()
{
	/*
	Корни, высчитанные "ручками":
	x1 = -2.36775971459153
	x2 = -0.435066997696302
	x3 = 2.93000063277712
	*/
	cout << "|\t\tIteration method" << setw(15) << "|" << endl;

	int counterIT = 0;
	int counterB = 0;
	iteration(-3.0, -2.0);
	iteration(-1.0, 0.0);
	iteration(2.5, 3);

	cout << "|\t\tBisection method" << setw(15) << "|" << endl;
	bisection(-3.0, -2.0);
	bisection(-1.0, 0.0);
	bisection(2.5, 3.0);

	cout << "|\t\tSpeed" << setw(26) << "|" << endl;
	comparison(-3.0, -2.0, counterIT, counterB);

	system("pause");
	return 0;
}