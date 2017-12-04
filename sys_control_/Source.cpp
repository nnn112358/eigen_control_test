#include "Eigen/Dense"
#include "Eigen/LU"
#include "Eigen/Core"
#include <iostream>
using namespace Eigen;
using namespace std;



void rk4(MatrixXd &dX, MatrixXd &X, MatrixXd &u, double &dt,double &tt,MatrixXd &A, MatrixXd &B, MatrixXd &C, MatrixXd &D) {


	MatrixXd dX1 = A*X + B*u;
	MatrixXd dX2 = A*X +B*u;



	X = dX*dt + X;
	cout << dX << endl;
}


int main(){
	double	k = 3.0;
	double	m = 0.1;
	double	c = 0.01;


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
	cout <<D << endl;
	


	double dt = 0.01;
	double tt = 0.0;
	MatrixXd X(2, 1);
	X(0, 0) = 1;
	X(1, 0) = 0;
	MatrixXd dX(2, 1);
	dX(0, 0) = 0;
	dX(1, 0) = 0;
	MatrixXd u(1, 1);
	u(0, 0) = 0;


	for (int i=0;i<1000;i++) {
		rk4(dX, X, u, tt,dt,A, B, C, D);
	
	
		tt += dt;
	}






	return 0;
}