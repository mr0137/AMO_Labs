#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <iomanip>

#pragma warning (disable : 4244)
using namespace std;

#define A  0.0
#define B  1.0
#define p_y0  sqrt(2)
#define p_z0  sqrt(2)/2


double fy(double x, double y, double z) 
{
	return z;
}

double fz(double x, double y, double z) 
{
	return (y - 1/(y*y*y));
}

double yo(double x) 
{
	return sqrt(1 + pow(M_E,2*x));
}

void R_Kutta(double &x, double &y, double &z, double h) 
{
	double k1, k2, k3, k4, g1, g2, g3, g4;

	k1 = h * fz(x, y, z);
	g1 = h * fy(x, y, z);

	k2 = h * fz(x + h / 2.0, y + g1 / 2.0, z + k1 / 2.0);
	g2 = h * fy(x + h / 2.0, y + g1 / 2.0, z + k1 / 2.0);

	k3 = h * fz(x + h / 2.0, y + g2 / 2.0, z + k2 / 2.0);
	g3 = h * fy(x + h / 2.0, y + g2 / 2.0, z + k2 / 2.0);

	k4 = h * fz(x + h, y + g3, z + k3);
	g4 = h * fy(x + h, y + g3, z + k3);

	x += h;
	z += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
	y += (g1 + 2.0 * g2 + 2.0 * g3 + g4) / 6.0;
}

void R_Kutta_full(double &x, double &y, double &z, double h) 
{
	int i, n = (B - A) / h;

	for (i = 0; i < n; i++)
	{
		R_Kutta(x, y, z, h);
	}
}

void Calculate(double x, double &y, double z, double fluff, double &h2)
{
	double y2, y1 = y;;

	do {
		y2 = y1;
		x = A;
		y = p_y0;
		z = p_z0;
		h2 /= 2;

		R_Kutta_full(x, y, z, h2);

		y1 = y;
	} while (fabs(y2 - y1) > 15 * fluff);
}

int main()
{
	setlocale(LC_ALL, "Rus");
	double x = A, y = p_y0, z = p_z0, fluff, h2 = 2.0, h = (B - A) / 10.0;
	int i, n = (B - A) / h;

	R_Kutta_full(x, y, z, h);
	cout << "|===|==========|==========|===========|==========|" << endl;
	cout << "|" << setw(3)  << "X" << "|" << setw(10) << "Y func" << "|" << setw(10) << "Y R_Kutta" << "|" << setw(11) << "Fluff" << "|" << setw(10) << "H" << "|" << endl;
	cout << "|===|==========|==========|===========|==========|" << endl;
	cout << "|" << setw(3) << setprecision(3) << x << "|" << setw(10) << setprecision(9) << yo(x) << "|" << setw(10) << setprecision(9) << y << "|" << setw(11) << setprecision(6) << (fluff = fabs(yo(x) - y)) << "|" << setw(10) << setprecision(10) << h << "|" << endl;
	cout << "|===|==========|==========|===========|==========|" << endl;

	R_Kutta_full(x, y, z, h2);
	
	Calculate(x, y, z, fluff, h2);

	cout << "\t\tRuange" << endl;
	cout << "|===|==========|==========|===========|==========|" << endl;
	cout << "|" << setw(3) << "X" << "|" << setw(10) << "Y func" << "|" << setw(10) << "Y R_Kutta" << "|" << setw(11) << "Fluff" << "|" << setw(10) << "H" << "|" << endl;
	cout << "|===|==========|==========|===========|==========|" << endl;
	cout << "|" << setw(3) << setprecision(3) << x << "|" << setw(10) << setprecision(9) << yo(x) << "|" << setw(10) << setprecision(9) << y << "|" << setw(11) << setprecision(6) << (fabs(yo(x) - y)) << "|" << setw(10) << setprecision(10) << h2 << "|" << endl;
	cout << "|===|==========|==========|===========|==========|" << endl;
	
	x = A;
	y = p_y0;
	z = p_z0;

	cout << "\t\tValue Table" << endl;
	cout << "|===|==========|==========|" << endl;
	cout << "|" << setw(3) << "X" << "|" << setw(10) << "Y func" << "|" << setw(10) << "YR_Kutta" << "|" << endl;
	cout << "|===|==========|==========|" << endl;

	for (i = 0; i < n; i++) 
	{
		R_Kutta(x, y, z, h);
		cout << "|" << setw(3) << setprecision(3) << x << "|" << setw(10) << setprecision(9) << yo(x) << "|" << setw(10) << setprecision(9) << y << "|" << endl;
	}

	cout << "|===|==========|==========|" << endl;
	system("pause");
	return 0;
}