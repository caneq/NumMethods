#include "pch.h"
#include "EulerMethods.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

typedef double(*func)(const vector<double> &, double);

void plot(const vector<vector<double>>& explicitPoints, const vector<vector<double>>& implicitPoints) {
	ofstream imp("implicitPoints.txt");
	ofstream exp("explicitPoints.txt");

	for (int i = 0; i < explicitPoints[0].size(); i++) {
		for (int j = 0; j < explicitPoints.size(); j++) {
			exp << explicitPoints[j][i] << " ";
		}
		exp << "\n";
	}

	for (int i = 0; i < implicitPoints[0].size(); i++) {
		for (int j = 0; j < implicitPoints.size(); j++) {
			imp << implicitPoints[j][i] << " ";
		}
		imp << "\n";
	}

	imp.close();
	exp.close();
}

double f1(const vector<double> &yk, double tk) {
	//return -yk[0] * yk[1] + sin(tk) / tk;
	//return -yk[0] * yk[1] + 1 - tk*tk/6 + tk*tk*tk*tk/25;
	//return yk[1] - (2 * yk[0] + tk * yk[1])*yk[0];
	return  2 * yk[1] * yk[2];
}

double f2(const vector<double> &yk, double tk) {
	//return -yk[1] * yk[1] + tk/(1 + tk*tk);
	//return exp(yk[0]) - (yk[0] + 2 * yk[1])*yk[0];
	return 1.5* yk[0] * yk[2];
}

double f3(const vector<double> &yk, double tk) {
	//return -yk[1] * yk[1] + tk / (1 + tk * tk);
	//return exp(yk[0]) - (yk[0] + 2 * yk[1])*yk[0];
	return -1*yk[0] * yk[1];
}

double A[3][3];
double b[3];

double hf1(const vector<double> &yk, double tk) {
	return  A[0][0] * yk[0] + A[0][1] * yk[1] + A[0][2] * yk[2] - b[0];
}
double hf2(const vector<double> &yk, double tk) {
	return  A[1][0] * yk[0] + A[1][1] * yk[1] + A[1][2] * yk[2] - b[1];
}
double hf3(const vector<double> &yk, double tk) {
	return  A[2][0] * yk[0] + A[2][1] * yk[1] + A[2][2] * yk[2] - b[2];
}

void hard() {
	double lambda[3] = { -100,-20,-30 };

	A[0][0] = (2 * lambda[0] + 4 * lambda[1]) / 6;
	A[0][1] = A[0][2] = A[1][0] = A[2][0] = 2 * (lambda[0] - lambda[1]) / 6;
	A[1][1] = A[2][2] = (2 * lambda[0] + lambda[1] + 3 * lambda[2]) / 6;
	A[1][2] = A[2][1] = (2 * lambda[0] + lambda[1] - 3 * lambda[2]) / 6;
	b[0] = -(4 * lambda[0] + 2 * lambda[1]) / 6;
	b[1] = -(4 * lambda[0] - lambda[1] - 9 * lambda[2]) / 6;
	b[2] = -(4 * lambda[0] - lambda[1] + 9 * lambda[2]) / 6;

	vector<func> f; f.push_back(hf1); f.push_back(hf2); f.push_back(hf3);
	vector<double> u0; 	u0.push_back(10); u0.push_back(22), u0.push_back(9);
	vector<vector<double>> explicitPoints, implicitPoints;
	vector<double> point0; point0.push_back(0);
	for (auto i : u0) {
		point0.push_back(i);
	}
	explicitPoints.push_back(point0);
	implicitPoints.push_back(point0);
	double T = 1;
	explicitEulerMethod(f, T, 1e-2, 1, u0, explicitPoints);
	implicitEulerMethod(f, T, 1e-2, 0.001, 1, u0, implicitPoints);
	plot(explicitPoints, implicitPoints);
}

void usual() {
	vector<func> f; f.push_back(f1); f.push_back(f2); f.push_back(f3);
	vector<double> u0; 	u0.push_back(1.0); u0.push_back(1), u0.push_back(1);
	vector<vector<double>> explicitPoints, implicitPoints;
	vector<double> point0; point0.push_back(0);
	for (auto i : u0) {
		point0.push_back(i);
	}
	explicitPoints.push_back(point0);
	implicitPoints.push_back(point0);
	explicitEulerMethod(f, 1, 1e-3, 0.1, u0, explicitPoints);
	implicitEulerMethod(f, 1, 1e-3, 0.001, 0.1, u0, implicitPoints);
	plot(explicitPoints, implicitPoints);
}


int main() {
	hard();
	//usual();

}