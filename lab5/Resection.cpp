
#include "Resection.h"

Resection::Resection(const vector<Point3D> &coords_object, const vector<Point2D> &coords_image, double c) {
	this->num_points = coords_object.size();
	this->coords_object = coords_object;
	this->coords_image = coords_image;
	this->c = c;
}

void Resection::computeResection(const Point3D &_T, const Angles &_ang, const Matrix &tolerances) {
	T = _T;
	Angles ang = _ang;
	M.rotate(ang);

	Matrix obs = convertToVector(coords_image);
	Matrix est;

	double threshold = 1e-6;
	do {
		A = resectionA();
		est = resectionCond();
		w = misclosure(obs, est);
		del = delta(A, w);

		T.x += del.at(0, 0);
		T.y += del.at(1, 0);
		T.z += del.at(2, 0);
		ang.omega += del.at(3, 0);
		ang.phi += del.at(4, 0);
		ang.kappa += del.at(5, 0);

		M.rotate(ang);
	} while (!belowTolerances(del, tolerances));

	v = residuals(A, del, w);
}

Matrix Resection::resectionA() {

	Matrix A(num_points * 2, 6);

	Point3D op;

	double omega = M.getOmega();
	double phi = M.getPhi();
	double kappa = M.getKappa();

	for (unsigned int i = 0; i < num_points; i++) {
		op = coords_object[i];

		// define sine and cosine values for omega, phi, and kappa
		// so to reduce computation time
		vector<double> sin_vals = { sin(omega), sin(phi), sin(kappa) };
		vector<double> cos_vals = { cos(omega), cos(phi), cos(kappa) };

		computeRowX(A, 2 * i, op, sin_vals, cos_vals);
		computeRowY(A, 2 * i + 1, op, sin_vals, cos_vals);
	}

	return A;
}

Matrix Resection::resectionCond() {

	Matrix cond(num_points * 2, 1);

	Point3D op;

	for (unsigned int i = 0; i < num_points; i++) {
		op = coords_object[i];
		
		cond.at(2 * i, 0) = -c * U(op) / W(op);
		cond.at(2 * i + 1, 0) = -c * V(op) / W(op);
	}

	return cond;
}

void Resection::computeRowX(Matrix &A, int row, const Point3D &op, const vector<double> &sinv, const vector<double> &cosv) {
	double u = U(op);
	double v = V(op);
	double w = W(op);
	double coeff = -c / pow(w, 2);
	Point3D diff = op.difference(T);

	A.at(row, 0) = coeff * (M.at(2, 0) * u - M.at(0, 0) * w);
	A.at(row, 1) = coeff * (M.at(2, 1) * u - M.at(0, 1) * w);
	A.at(row, 2) = coeff * (M.at(2, 2) * u - M.at(0, 2) * w);

	A.at(row, 3) = coeff * (diff.y * (u * M.at(2, 2) - w * M.at(0, 2)) -
							diff.z * (u * M.at(2, 1) - w * M.at(0, 1)));

	A.at(row, 4) = coeff * (diff.x * (-w * sinv.at(1) * cosv.at(2)				- u * cosv.at(1)) +
							diff.y * ( w * sinv.at(0) * cosv.at(1) * cosv.at(2) - u * sinv.at(0) * sinv.at(1)) +
							diff.z * (-w * cosv.at(0) * cosv.at(1) * cosv.at(2) + u * cosv.at(0) * sinv.at(1)));

	A.at(row, 5) = -c * v / w;
}

void Resection::computeRowY(Matrix &A, int row, const Point3D &op, const vector<double> &sinv, const vector<double> &cosv) {
	double u = U(op);
	double v = V(op);
	double w = W(op);
	double coeff = -c / pow(w, 2);
	Point3D diff = op.difference(T);

	A.at(row, 0) = coeff * (M.at(2, 0) * v - M.at(1, 0) * w);
	A.at(row, 1) = coeff * (M.at(2, 1) * v - M.at(1, 1) * w);
	A.at(row, 2) = coeff * (M.at(2, 2) * v - M.at(1, 2) * w);

	A.at(row, 3) = coeff * (diff.y * (v * M.at(2, 2) - w * M.at(1, 2)) -
							diff.z * (v * M.at(2, 1) - w * M.at(1, 1)));

	A.at(row, 4) = coeff * (diff.x * ( w * sinv.at(1) * sinv.at(2)				- v * cosv.at(1)) +
							diff.y * (-w * sinv.at(0) * cosv.at(1) * sinv.at(2) - v * sinv.at(0) * sinv.at(1)) +
							diff.z * ( w * cosv.at(0) * cosv.at(1) * sinv.at(2) + v * cosv.at(0) * sinv.at(1)));

	A.at(row, 5) = c * u / w;
}

double Resection::U(const Point3D &point) {
	Point3D diff = point.difference(T);
	return M.at(0, 0) * diff.x + M.at(0, 1) * diff.y + M.at(0, 2) * diff.z;
}

double Resection::V(const Point3D &point) {
	Point3D diff = point.difference(T);
	return M.at(1, 0) * diff.x + M.at(1, 1) * diff.y + M.at(1, 2) * diff.z;
}

double Resection::W(const Point3D &point) {
	Point3D diff = point.difference(T);
	return M.at(2, 0) * diff.x + M.at(2, 1) * diff.y + M.at(2, 2) * diff.z;
}

Matrix Resection::getA() {
	return A;
}

Matrix Resection::getMisclosure() {
	return w;
}

Matrix Resection::getDelta() {
	return del;
}

Matrix Resection::getResiduals() {
	return v;
}

vector<Point3D> Resection::getObjectCoords() {
	return coords_object;
}

vector<Point2D> Resection::getImageCoords() {
	return coords_image;
}

RotationMatrix Resection::getM() {
	return M;
}

Point3D Resection::getT() {
	return T;
}

double Resection::getFocalLength() {
	return c;
}