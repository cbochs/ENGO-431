#include "AbsoluteOrientation.h"

AbsoluteOrientation::AbsoluteOrientation(const vector<Point3D>& coords_object, const vector<Point3D>& coords_model) {
	this->num_points = coords_object.size();
	this->coords_object = coords_object;
	this->coords_model = coords_model;
}

void AbsoluteOrientation::computeOrientation(const Point3D &_T, const Angles &_ang, double _lambda) {
	T = _T;
	Angles ang = _ang;
	lambda = _lambda;
	M.rotate(ang);

	Matrix obs = convertToVector(coords_object);
	Matrix est;

	double threshold = 1e-6;
	do {
		A = absoluteA();
		est = absoluteCond();
		w = misclosure(obs, est);
		del = delta(A, w);
		
		ang.omega += del.at(0, 0);
		ang.phi += del.at(1, 0);
		ang.kappa += del.at(2, 0);
		lambda += del.at(3, 0);
		T.x += del.at(4, 0);
		T.y += del.at(5, 0);
		T.z += del.at(6, 0);

		M.rotate(ang);
	} while (del.maxAbsElem() > threshold);

	v = residuals(A, del, w);
}

Matrix AbsoluteOrientation::absoluteA() {

	Matrix A(num_points * 3, 7);

	Point3D pm;

	double omega = M.getOmega();
	double phi = M.getPhi();
	double kappa = M.getKappa();

	for (unsigned int i = 0; i < num_points; i++) {
		pm = coords_model[i];

		// define sine and cosine values for omega, phi, and kappa
		// so to reduce computation time
		vector<double> sin_vals = { sin(omega), sin(phi), sin(kappa) };
		vector<double> cos_vals = { cos(omega), cos(phi), cos(kappa) };
		
		computeRowX(A, 3 * i, pm, sin_vals, cos_vals);
		computeRowY(A, 3 * i + 1, pm, sin_vals, cos_vals);
		computeRowZ(A, 3 * i + 2, pm, sin_vals, cos_vals);
	}

	return A;
}

Matrix AbsoluteOrientation::absoluteCond() {

	Matrix cond(num_points * 3, 1);

	Point3D pm;

	for (unsigned int i = 0; i < num_points; i++) {
		pm = coords_model[i];
		pm.transform(lambda, M, T);

		cond.at(3 * i, 0) = pm.x;
		cond.at(3 * i + 1, 0) = pm.y;
		cond.at(3 * i + 2, 0) = pm.z;
	}

	return cond;
}

void AbsoluteOrientation::computeRowX(Matrix &A, int row, const Point3D pm, const vector<double> &sinv, const vector<double> &cosv) {
	A.at(row, 0) = lambda * (pm.y * (-sinv.at(0) * sinv.at(2) + cosv.at(0) * sinv.at(1) * cosv.at(2)) +
							 pm.z * ( cosv.at(0) * sinv.at(2) + sinv.at(0) * sinv.at(1) * cosv.at(2)));

	A.at(row, 1) = lambda * (-pm.x * sinv.at(1)	* cosv.at(2) +
							  pm.y * sinv.at(0) * cosv.at(1) * cosv.at(2) -
							  pm.z * cosv.at(0) * cosv.at(1) * cosv.at(2));

	A.at(row, 2) = lambda * (-pm.x *  cosv.at(1) * sinv.at(2) +
							  pm.y * (cosv.at(0) * cosv.at(2) - sinv.at(0) * sinv.at(1) * sinv.at(2)) +
							  pm.z * (sinv.at(0) * cosv.at(2) + cosv.at(0) * sinv.at(1) * sinv.at(2)));

	A.at(row, 3) = pm.x * M.at(0, 0) + pm.y * M.at(0, 1) + pm.z * M.at(0, 2);

	A.at(row, 4) = 1;
	A.at(row, 5) = 0;
	A.at(row, 6) = 0;
}

void AbsoluteOrientation::computeRowY(Matrix &A, int row, const Point3D pm, const vector<double> &sinv, const vector<double> &cosv) {
	A.at(row, 0) = lambda * (pm.y * (-sinv.at(0) * cosv.at(2) - cosv.at(0) * sinv.at(1) * sinv.at(2)) +
							 pm.z * ( cosv.at(0) * cosv.at(2) - sinv.at(0) * sinv.at(1) * sinv.at(2)));

	A.at(row, 1) = lambda * (pm.x * sinv.at(1) * sinv.at(2) -
							 pm.y * sinv.at(0) * cosv.at(1) * sinv.at(2) +
							 pm.z * cosv.at(0) * cosv.at(1) * sinv.at(2));

	A.at(row, 2) = lambda * (-pm.x *   cosv.at(1) * cosv.at(2) +
							  pm.y * (-cosv.at(0) * sinv.at(2) - sinv.at(0) * sinv.at(1) * cosv.at(2)) +
							  pm.z * (-sinv.at(0) * sinv.at(2) + cosv.at(0) * sinv.at(1) * cosv.at(2)));

	A.at(row, 3) = pm.x * M.at(1, 0) + pm.y * M.at(1, 1) + pm.z * M.at(1, 2);

	A.at(row, 4) = 0;
	A.at(row, 5) = 1;
	A.at(row, 6) = 0;
}

void AbsoluteOrientation::computeRowZ(Matrix &A, int row, const Point3D pm, const vector<double> &sinv, const vector<double> &cosv) {
	A.at(row, 0) = lambda * (-pm.y * cosv.at(0) * cosv.at(1) -
							  pm.z * sinv.at(0) * cosv.at(1));

	A.at(row, 1) = lambda * (pm.x * cosv.at(1) +
							 pm.y * sinv.at(0) * sinv.at(1) -
							 pm.z * cosv.at(0) * sinv.at(1));

	A.at(row, 2) = 0;

	A.at(row, 3) = pm.x * M.at(2, 0) + pm.y * M.at(2, 1) + pm.z * M.at(2, 2);

	A.at(row, 4) = 0;
	A.at(row, 5) = 0;
	A.at(row, 6) = 1;
}

Matrix AbsoluteOrientation::getA() {
	return A;
}

Matrix AbsoluteOrientation::getMisclosure() {
	return w;
}

Matrix AbsoluteOrientation::getDelta() {
	return del;
}

Matrix AbsoluteOrientation::getResiduals() {
	return v;
}

vector<Point3D> AbsoluteOrientation::getObjectCoords() {
	return coords_object;
}

vector<Point3D> AbsoluteOrientation::getModelCoords() {
	return coords_model;
}

RotationMatrix AbsoluteOrientation::getM() {
	return M;
}

Point3D AbsoluteOrientation::getT() {
	return T;
}

double AbsoluteOrientation::getScale() {
	return lambda;
}