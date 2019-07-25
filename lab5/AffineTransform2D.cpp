
#include "AffineTransform2D.h"

AffineTransform2D::AffineTransform2D(vector<Point2D> &coords_from, vector<Point2D> &coords_to) {
	setCoordinates(coords_from, coords_to);
}

void AffineTransform2D::setCoordinates(vector<Point2D> &coords_from, vector<Point2D> &coords_to) {
	this->coords_from = coords_from;
	this->coords_to = coords_to;

	num_points = coords_to.size();

	computeTransformation();
}

void AffineTransform2D::computeTransformation() {
	affineA();

	Matrix observations(2 * num_points, 1);
	for (unsigned int i = 0; i < num_points; i++) {
		observations[2 * i][0] = coords_to[i].x;
		observations[2 * i + 1][0] = coords_to[i].y;
	}

	del = delta(A, -1 * observations);
	v = residuals(A, del, -1 * observations);

	params.a = del.at(0, 0);
	params.b = del.at(1, 0);
	params.dx = del.at(2, 0);
	params.c = del.at(3, 0);
	params.d = del.at(4, 0);
	params.dy = del.at(5, 0);

	params.sx = sqrt(pow(params.a, 2) + pow(params.c, 2));
	params.sy = sqrt(pow(params.b, 2) + pow(params.d, 2));
	params.theta = atan2(params.c, params.a);
	params.delta = atan2(params.a * params.b + params.c * params.d, params.a * params.d - params.b * params.c);
}

void AffineTransform2D::affineA() {
	A.resize(2 * num_points, 6);
	A.clear();

	for (unsigned int i = 0; i < num_points; i++) {
		A[2 * i][0] = coords_from[i].x;
		A[2 * i][1] = coords_from[i].y;
		A[2 * i][2] = 1.0;
		A[2 * i + 1][3] = coords_from[i].x;
		A[2 * i + 1][4] = coords_from[i].y;
		A[2 * i + 1][5] = 1.0;
	}
}

Matrix AffineTransform2D::getA() {
	return A;
}

Matrix AffineTransform2D::getDelta() {
	return del;
}

Matrix AffineTransform2D::getResiduals() {
	return v;
}

vector<Point2D> AffineTransform2D::getFromCoords() {
	return coords_from;
}

vector<Point2D> AffineTransform2D::getToCoords() {
	return coords_to;
}

AffineParams AffineTransform2D::getParams() {
	return params;
}