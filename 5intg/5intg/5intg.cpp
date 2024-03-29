#include "pch.h"
#include <iostream>
#include <iomanip>

using namespace std;

double f(double x) {
	return x*x;
}

double f2(double x, double y) {
	return x * x / (1 + y * y);
}

double trapn(double f(double), double a, double b, unsigned int n) {
	double h = (b - a) / n;
	double sum = 0;
	double x = a + h;
	for (int i = 1; i < n; i++) {
		sum += 2 * f(x);
		x += h;
	}
	sum += f(b) + f(a);
	sum *= h / 2;
	return sum;
}

double trap(double f(double), double a, double b, double eps) {
	if (a > b) swap(a, b);
	if (b - a <= DBL_EPSILON) return 0;
	unsigned int n = 2;
	double sum2 = trapn(f, a, b, n);
	double sum;
	double eps3 = eps * 3;
	do {
		sum = sum2;
		n *= 2;
		sum2 = trapn(f, a, b, n);
	} while (fabs(sum - sum2) >= eps3);
	return sum2;
}

double simpsm(double f(double), double a, double b, int m) {

	double h = (b - a) / m / 2;
	double h2 = h * 2;

	double sum = 0;
	double tmp = 0;

	double x = a + h;
	for (int i = 0; i < m; i++) {
		tmp +=  f(x);
		x += h2;
	}

	sum += 4 * tmp;
	tmp = 0;
	x = a + h2;
	for (int i = 1; i < m; i++) {
		tmp += f(x);
		x += h2;
	}

	sum += 2 * tmp;
	sum += f(a) + f(b);
	sum *= h / 3;
	return sum;
}

double simps(double f(double), double a, double b, double eps) {
	bool swapped = false;
	if (a > b) {
		swap(a, b);
		swapped = true;
	}
	if (b - a <= DBL_EPSILON) return 0;
	unsigned int m = 2;
	double sum2 = simpsm(f, a, b, m);
	double sum;
	double eps15 = eps * 15;
	do {
		sum = sum2;
		m *= 2;
		sum2 = simpsm(f, a, b, m);
	} while (fabs(sum - sum2) >= eps15);
	if (swapped) sum2 = -sum2;
	return sum2;
}

double kubsimpsMN(double f(double, double), double a, double b, double c, double d, int M, int N) {
	int N2 = 2 * N;
	int M2 = 2 * M;
	double hx = (b - a) / N2;
	double hy = (d - c) / M2;
	double sum = 0;

	for (int i2 = 0; i2 < N2; i2+=2) {
		for (int j2 = 0; j2 < M2; j2+=2) {
			sum += f(a + hx * i2, c + hy * j2);
			sum += 4 * f(a + hx * (i2 + 1), c + hy * j2);
			sum += f(a + hx * (i2 + 2), c + hy * j2);
			sum += 4 * f(a + hx * i2, c + hy * (j2 + 1));
			sum += 16 * f(a + hx * (i2 + 1), c + hy * (j2 + 1));
			sum += 4 * f(a + hx * (i2 + 2), c + hy * (j2 + 1));
			sum += f(a + hx * i2, c + hy * (j2 + 2));
			sum += 4 * f(a + hx * (i2 + 1), c + hy * (j2 + 2));
			sum += f(a + hx * (i2 + 2), c + hy * (j2 + 2));

		}
	}

	sum *= hx * hy / 9;
	return sum;
}

double kubsimps(double f(double, double), double a, double b, double c, double d, double eps) {
	bool swapped1 = false;
	if (a > b) {
		swap(a, b);
		swapped1 = true;
	}
	bool swapped2 = false;
	if (c > d) {
		swap(c, d);
		swapped2 = true;
	}
	if (swapped1 && swapped2) swapped1 = swapped2 = false;
	if (fabs(b - a) <= DBL_EPSILON) return 0;
	if (fabs(d - c) <= DBL_EPSILON) return 0;
	unsigned int M = 2;
	double sum2 = kubsimpsMN(f, a, b, c, d, M, M);
	double sum;
	double eps15 = eps * 15;
	do {
		sum = sum2;
		M *= 2;
		sum2 = kubsimpsMN(f, a, b, c, d, M, M);
	} while (fabs(sum - sum2) >= eps15);
	if (swapped1 || swapped2) sum2 = -sum2;
	return sum2;
}

int main() {
	setlocale(LC_ALL, "rus");
	cout << "Формула трапеций            " << setprecision(15) << trap(f, 0, 1, 1e-6) << endl;
	cout << "Формула Симпсона            " << setprecision(15) << simps(f, 0, 1, 1e-6) << endl;
	cout << "Кубатурная формула Симпсона " << setprecision(15) << kubsimps(f2, 0, 4, 1, 2, 1e-10) << endl;
}