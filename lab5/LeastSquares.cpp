
#include "LeastSquares.h"

Matrix delta(Matrix &A, const Matrix &P, const Matrix &w) {
	Matrix Qx = cofactorMatrix(A, P);

	return -1 * Qx * A.trans() * P * w;
}

Matrix delta(Matrix &A, const Matrix &w) {
	Matrix Qx = cofactorMatrix(A);

	return -1 * Qx * A.trans() * w;
}

bool belowTolerances(const Matrix &delta, const Matrix &tolerances) {
	for (unsigned int i = 0; i < delta.getrows(); i++) {
		if (delta.at(i, 0) > tolerances.at(i, 0))
			return false;
	}

	return true;
}

Matrix misclosure(const Matrix &obs, const Matrix &est) {
	return est - obs;
}

Matrix residuals(Matrix &A, const Matrix &delta, const Matrix &w) {
	return A * delta + w;
}

double aposteriori(Matrix &v, const Matrix &P, double dof) {
	return (v.trans() * P * v).at(0, 0) / dof;
}

double aposteriori(Matrix &v, double dof) {
	return (v.trans() * v).at(0, 0) / dof;
}

Matrix cofactorMatrix(Matrix &A, const Matrix &P) {
	Matrix N = A.trans() * P * A;
	return N.inv();
}

Matrix cofactorMatrix(Matrix &A) {
	Matrix N = A.trans() * A;
	return N.inv();
}

Matrix unknownCovariance(Matrix &A, const Matrix &P, double aposteriori) {
	Matrix Qx = cofactorMatrix(A, P);
	return aposteriori * Qx;
}

Matrix unknownCovariance(Matrix &A, double aposteriori) {
	Matrix Qx = cofactorMatrix(A);
	return aposteriori * Qx;
}

Matrix observeCovariance(Matrix &A, const Matrix &P) {
	Matrix Cx = cofactorMatrix(A, P);
	return A.trans() * Cx * A;
}

Matrix observeCovariance(Matrix &A) {
	Matrix Cx = cofactorMatrix(A);
	return A * Cx * A.trans();
}

Matrix residualCovariance(const Matrix &C, Matrix &A, const Matrix &P) {
	Matrix Cl = observeCovariance(A, P);
	return C - Cl;
}

Matrix residualCovariance(Matrix &A) {
	Matrix Cl = observeCovariance(A);
	return eye(Cl.getrows()) - Cl;
}

Matrix correlation(Matrix &A) {
	Matrix Qx = cofactorMatrix(A);
	int n = Qx.getrows();

	Matrix corr;
	corr.resize(n, n);

	for (unsigned int j = 0; j < Qx.getrows(); j++) {
		for (int i = 0; i < n; i++) {
			corr[j][i] = Qx[j][i] / sqrt(Qx[j][j] * Qx[i][i]);
		}
	}

	return corr;
}

Matrix eye(int size) {
	Matrix eye;
	eye.resize(size, size);

	for (int i = 0; i < size; i++) {
		eye[i][i] = 1;
	}

	return eye;
}