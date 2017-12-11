#include "Eigen/Dense"
#include "Eigen/Core"
#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES

#include <math.h>
using namespace Eigen;
using namespace std;

void RungeKutta(MatrixXd dX, MatrixXd &X, MatrixXd u, double tt, double dt, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D);
MatrixXd RiccatiSolver(MatrixXd A, MatrixXd B, MatrixXd Q, MatrixXd R);

int main() {
	// Parameters defining the system
	double m = 250.0;         // system mass
	double 	k = 40.0;          // spring constant
	double  b = 60.0;          // damping constant

	MatrixXd A(3, 3);
	A << 1, -1, 1, 1, -k / m, -b / m, 1, 1, 1;
	MatrixXd B(3, 1);
	B << 0, -1 / m, 1;

	MatrixXd C(1, 3);
	C << 1, 0, 1;

	MatrixXd D(1, 3);
	D << 1, 0, 0;

	MatrixXd K(1, 3);		//FeedBask Gain
	K << 15.80864382, -3.63298165, 7.85453193;

	double dt = 0.01;
	double tt = 0.0;
	MatrixXd X(3, 1);
	X << 10, 0, 0;
	MatrixXd dX(3, 1);
	dX << 0, 0, 0;

	MatrixXd u(1, 1);
	u << 0;
	MatrixXd Y(1, 1);
	Y << 0;
	MatrixXd R(3, 3), Q(3, 3);
	Q << 300, 0,0,0, 0, 60,0,0,0;
	R << 1;

	MatrixXd P = RiccatiSolver(A, B, Q, R);

	cout << "P : " << endl << P << endl;

	ofstream ofs("test.csv");
	ofs << "time," << ",u,y,x0,x1,x2" << endl;
	for (int i = 0; i < 1000; i++) {
		RungeKutta(dX, X, u, tt, dt, A, B, C, D);
		Y = C*X;
		u = -K*X;
		ofs << tt << "," << u(0, 0) << "," << Y(0, 0) << "," << X(0, 0) << "," << X(1, 0) << "," << X(2, 0) << endl;
		tt += dt;
	}
	return 0;

}


void RungeKutta(MatrixXd dX, MatrixXd &X, MatrixXd u, double tt, double dt, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D) {
	MatrixXd k1 = A*X + B*u;
	MatrixXd k2 = A*(X + 0.5*k1*dt) + B*u;
	MatrixXd k3 = A*(X + 0.5*k2*dt) + B*u;
	MatrixXd k4 = A*(X + k3*dt) + B*u;
	MatrixXd k = (k1 + 2.0*k2 + 2.0*k3 + k4)*dt / 6.0;
	X = X + k;
}

//���b�J�`������������
MatrixXd RiccatiSolver(MatrixXd A, MatrixXd B, MatrixXd Q, MatrixXd R) {

	int n = (int)A.rows();

	//�n�~���g���s��𐶐�����
	MatrixXd H(2 * n, 2 * n);
	H << A, -B*R.inverse()*B.transpose(), -Q, -A.transpose();

	//�ŗL�l�ƌŗL�x�N�g�������߂�
	EigenSolver<Matrix<double, 4, 4>> es(H);
	if (es.info() != Success) abort();

	//�x�N�g���Z�b�gu,v�������o��
	vector<VectorXd> v_set, u_set;
	for (int i = 0; i<2 * n; i++) {
		if (es.eigenvalues().real()(i) < 0) {
			v_set.push_back(es.eigenvectors().real().block(0, i, n, 1));
			u_set.push_back(es.eigenvectors().real().block(n, i, n, 1));
		}
	}

	int num = (int)v_set.size();
	MatrixXd v(n, num), u(n, num);
	for (int i = 0; i<num; i++) {
		v.block(0, i, n, 1) = v_set[i];
		u.block(0, i, n, 1) = u_set[i];
	}

	//��P�����߂�
	MatrixXd P = u * v.inverse();
	return P;
}