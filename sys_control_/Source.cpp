#include "Eigen/Dense"
#include "Eigen/Core"
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES

#include <math.h>
using namespace Eigen;
using namespace std;

void Euler(MatrixXd dX, MatrixXd &X, MatrixXd u, double tt, double dt, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D) {
	MatrixXd f1 = A*X + B*u;
	X = X + f1*dt;
}

void RungeKutta(MatrixXd dX, MatrixXd &X, MatrixXd u, double tt, double dt, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D) {
	MatrixXd k1 = A*X + B*u;
	MatrixXd k2 = A*(X + 0.5*k1*dt) + B*u;
	MatrixXd k3 = A*(X + 0.5*k2*dt) + B*u;
	MatrixXd k4 = A*(X + k3*dt) + B*u;
	MatrixXd k = (k1 + 2.0*k2 + 2.0*k3 + k4)*dt / 6.0;
	X = X + k;
}

int main() {
	double	k = 1.0;
	double	m = 0.1;
	double	c = 0.1;

	MatrixXd A(2, 2);
	A(0, 0) = 0;
	A(0, 1) = 1;
	A(1, 0) = -k / m;
	A(1, 1) = -c / m;
	MatrixXd B(2, 1);
	B(0, 0) = 1;
	B(1, 0) = 1 / m;

	MatrixXd C(1, 2);
	C(0, 0) = 1;
	C(0, 1) = 0;
	MatrixXd D(1, 1);
	D(0, 0) = 0;

	cout << "A=" << endl;
	cout << A << endl;
	cout << "B=" << endl;
	cout << B << endl;
	cout << "C=" << endl;
	cout << C << endl;
	cout << "D=" << endl;
	cout << D << endl;

	double dt = 0.01;
	double tt = 0.0;
	MatrixXd X(2, 1);
	X(0, 0) = 10;
	X(1, 0) = 0;
	MatrixXd dX(2, 1);
	dX(0, 0) = 0;
	dX(1, 0) = 0;
	MatrixXd u(1, 1);
	u(0, 0) = 0;
	MatrixXd Y(1, 1);
	Y(0, 0) = 0;


	double freq = 1.0;
	ofstream ofs("test.csv");
	ofs << "time," << "y" << endl;

	for (int i = 0; i < 1000; i++) {

		//u(0, 0) = sin(2 * M_PI*freq*tt);
		//Euler(dX, X, u, tt, dt, A, B, C, D);
		RungeKutta(dX, X, u, tt, dt, A, B, C, D);
		Y = C*X;

		//cout << tt << X(0, 0) << "," << X(1, 0) << endl;
		//ofs << tt << "," << u(0, 0) << "," << X(0, 0) << "," << X(1, 0) << endl;

		ofs << tt << "," << Y(0, 0) << endl;

		tt += dt;
	}




	return 0;
}