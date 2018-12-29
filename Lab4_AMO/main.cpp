#define _USE_MATH_DEFINES 
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <iostream>

#define A -1
#define B 10

using namespace std;

/*
	Вариант 8
*/

double y2(double x)
{
	return (x + 2) * pow(M_E, x);
}

double getH(double eps)
{
	double x = A;
	double m2 = 0;
	double h;

	while (x <= B)
	{
		if (y2(x) > m2)
		{
			m2 = y2(x);
		}
			x += eps;
	}

	h = sqrt(12 * eps / (B - A) / m2);
	return h;
}

double func(double x) 
{
	return x * pow(M_E,x);
}

double pfunc(double x) 
{
	return (x - 1) * pow(M_E, x);
}

double integrate(int n, double h) 
{
	double In = 0.0;

	for (int i = 1; i < n; i++) {
		In += func(A + i * h);
	}

	In += (func(A) + func(A + n * h)) / 2;

	return In*h;
}

double runge(double eps, double h, int n) 
{
	double In, I2n;
	
	In = integrate(n, h);
	I2n = integrate(2 * n, h);
	
	while (fabs(In - I2n) > 3 * eps)
	{
		In = I2n;
		n *= 2;
		h /= 2;
		I2n = integrate(n, h);
	}

	return I2n;
}

void first(double &h, double eps, int &n, double &fluff)
{
	double value;
	h = getH(eps);
	n = (int)((B - A) / h) + 1;
	value = pfunc(B) - pfunc(A);
	fluff = fabs(integrate(n, h) - value);
	cout << "\t\t\t Task 1" << endl;
	cout << "|===========|===============|=====================|=====================|" << endl;
	cout << "|" << setw(11) <<  "eps" << "|" << setw(15) << "h" << "|" << setw(21) << "Value" << "|" << setw(21) << "fluff" << "|" << endl;
	cout << "|===========|===============|=====================|=====================|" << endl;
	cout << "|" << setw(11) << setprecision(5) << eps << "|" << setw(15) << setprecision(10) << h << "|" << setw(21) << setprecision(10) << value << "|" << setw(21) << setprecision(10) << fluff << "|" << endl;
	cout << "|===========|===============|=====================|=====================|" << endl;
}

void second(double h, double eps, int n, double fluff)
{
	n =(int)(B - A) / sqrt(eps) + 1;
	h = double((B - A)) / n;
	double integral = runge(fluff, h, n);
	cout << "\t\t\t Task 2" << endl;
	cout << "\t\t Integral = " << integral << endl;


	cout << "|===========|===============|=====================|" << endl;
	cout << "|" << setw(11) << "eps" << "|" << setw(15) << "h" << "|" << setw(21) << "fluff" << "|" << endl;
	cout << "|===========|===============|=====================|" << endl;
	cout << "|" << setw(11) << setprecision(5) << fluff << "|" << setw(15) << setprecision(10) << h << "|" << setw(21) << setprecision(10) << fabs(integral - (pfunc(B) - pfunc(A))) << "|" << endl;
	cout << "|===========|===============|=====================|" << endl;
}

int main() 
{
	setlocale(LC_ALL, "Rus");
	double h=0, fluff=0, eps = 1e-5;
	int n;
	first(h, eps, n, fluff);
	second(h, eps, n, fluff);
 	system("pause");
	return 0;
}