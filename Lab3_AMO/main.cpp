/*
| 12  6  2 16 | 148 |
| 20 56 18 17 | 218 |
| 18  0 34 15 | 230 |
|  2  5 17 17 | 144 |
1) Схема с выбором главного елемента
2) Метод простой итерации
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#define n 4
#define m 5
#define eps 0.0001

using namespace std;

void copy(const double buf_matr[n][m], double matr[n][m])
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			matr[i][j] = buf_matr[i][j];
		}
	}
}

void print(const double x[n], const double matr[n][m]) 
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		sum = 0;
		
		for (int j = 0; j < m-1; j++)
		{
			sum += matr[i][j] * x[j];
			cout << "(" << matr[i][j] << " * " << x[j] << ")";
			if (j < m - 2) cout << " +";
			cout << " ";
		}
		cout << " = " <<sum  << endl;
		
	}
	for (int i = 0; i < n; i++)
	{
		cout << "X" << i + 1 << " = " << setw(7) << setprecision(5) << x[i] << endl;
	}
		
}

void print_matr(const double matr[n][m]) 
{
	for (int i = 0; i < n; i++) 
	{
		cout << "|=======|=======|=======|=======|=======|" << endl;
		for (int j = 0; j < m; j++)
		{
			cout << "|" << setw(7) << setprecision(5) <<matr[i][j];
		}
		cout << "|" << endl;
	}
	cout << "|=======|=======|=======|=======|=======|" << endl;
}

void Main_elem(double matr[n][m], double *res) 
{
	double x[n][m];
	double kf[n - 1];
	int i, j, z, k;
	for (k = 0; k < n; k++) 
	{
		double max = fabs(matr[0][0]);

		int im = 0, jm = 0;

		for (i = 0; i < n; i++)
		{
			for (j = 0; j < m - 1; j++)
			{
				if (fabs(matr[i][j]) > max)
				{
					im = i;
					jm = j;
					max = matr[i][j];
				}
			}
				
		}

		if (max == 0)
		{
			goto last;
		}

		for (i = 0, j = 0; i < n; i++) 
		{
			if (i == im) continue;
			kf[j] = -matr[i][jm] / matr[im][jm];
			j++;
		}

		for ( i = 0, z = 0; i < n; i++) 
		{
			if (i == im) continue;
			for (j = 0; j < m; j++)
			{
				matr[i][j] += matr[im][j] * kf[z];
			}	
			z++;
		}
		last:
		for (j = 0; j < m; j++) 
		{
			x[k][j] = matr[im][j];
			matr[im][j] = 0;
		}

	}

	print_matr(x);

	for (i = n - 1; i >= 0; i--) 
	{
		for (j = 0; j < m - 1; j++) 
		{
			if (x[i][j] != 0) 
			{
				res[j] = x[i][m - 1] / x[i][j];
				for (k = 0; k < n; k++) 
				{
					if (k == i) continue;
					x[k][m - 1] -= x[k][j] * res[j];
					x[k][j] = 0;
				}
				continue;
			}
		}
	}
}

void Simple_iter(double matr[n][m], double x[n])
{
	int i, j, it = 0;
	double temp2, q = 0.0, max;
	double temp[n], beta[n], alfa[n][m];

	for (i = 0; i < n; i++)
	{
		temp[i] = 0.0;
	}
	
	for (i = 0; i < n; i++)
	{
		beta[i] = matr[i][m - 1] / matr[i][i];
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i != j)
				alfa[i][j] = -matr[i][j] / matr[i][i];
			else
				alfa[i][j] = 0;
		}
		x[i] = beta[i];
	}
	
	for (i = 0; i < n; i++)
	{
		max = 0.0;
		for (j = 0; j < n; j++)
		{
			max = fabs(alfa[i][j]) + max;
		}
			
		if (max > q) q = max;
	}
	do
	{
		it++;
		for (i = 0; i < n; i++)
		{
			temp[i] = x[i];
		}

		for (i = 0; i < n; i++)
		{
			temp2 = 0.0;
			for (j = 0; j < n; j++)
			{
				if (i != j)
				{
					temp2 = alfa[i][j] * x[j] + temp2;
				}
					
			}	
			x[i] = beta[i] + temp2;
		}

		max = 0.0;
		temp2 = 0.0;
		for (i = 0; i < n; i++)
		{
			temp2 = fabs(temp[i] - x[i]);
			if (temp2 > max)
			{
				max = temp2;
			}
		}
	} while (max > fabs((1 - q) * eps / q));
	cout << "For eps = " << eps << " num of iterations is = " << it << endl;
}



int main() 
{
	setlocale(LC_ALL, "Rus");
	const double buf_matr[n][m] =
	{
		{12, 6, 2, 16, 148},
		{20, 56, 18, 17, 218},
		{16, 0, 34, 15, 230},
		{2, 5, 17, 17, 144}
	};

	/*
		Уравнение (1) и (4) не подходят к условию, исправляем :
		| 12  6  2 16 | 148 |  ---------->   | 12  6   2 16 | 148 |					  | 12  6   2   16 | 148 | (1) = (1) - (4) | 26 -10   0  -3 |  90 |
		| 20 56 18 17 | 218 |(2) = (2) - (3) |  2 56 -16  2 | -12 | --------------->  |  2 56 -16    2 | -12 | --------------->|  2  56 -16   2 | -12 |
		| 18  0 34 15 | 230 |(3) = (3) - (1) |  6 -6  32 -1 |  82 |					  |  6 -6  32   -1 |  82 |				   |  6  -6  32  -1 |  82 |
		|  2  5 17 17 | 144 |                |  2  5  17 17 | 144 |(4) = 2 *(4) - (3) | -2 16   2   35 | 206 |				   | -2  16   2  35 | 206 |
		*Я достаточно долго мучал эту СЛАР и это лучшее, к чему смог прийти.
		**Корни в двух методах между собой не сошлись...
	*/
	const double new_matr[n][m] =
	{
		{26,  -10,    0,    -3,    90},
		{ 2,   56,  -16,     2,   -12},
		{ 6,   -6,   32,    -1,    82},
		{-2,   16,    2,    35,   206}
	};
	double matr[n][m], x[n];
	copy(buf_matr, matr);
	cout << "Схема с выбором главного елемента:" << endl;
	cout << "Система:" << endl;
	print_matr(matr);
	cout << "Переделанная система:" << endl;
	Main_elem(matr, x);
	print(x, buf_matr);

	copy(new_matr, matr);
	cout << "Метод простой итерации:" << endl;
	cout << "Система:" << endl;
	print_matr(matr);
	cout << "Переделанная система:" << endl;
	Simple_iter(matr,x);
	print(x, new_matr);
	system("pause");
	return 0;
}