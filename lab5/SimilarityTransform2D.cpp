
#include "SimilarityTransform2D.h"

SimilarityTransform2D::SimilarityTransform2D(vector<Point2D> &coords_from, vector<Point2D> &coords_to) {
	setCoordinates(coords_from, coords_to);
}

void SimilarityTransform2D::setCoordinates(vector<Point2D> &coords_from, vector<Point2D> &coords_to) {
	this->coords_from = coords_from;
	this->coords_to = coords_to;

	num_points = coords_to.size();

	computeTransformation();
}

void SimilarityTransform2D::computeTransformation() {
	similarityA();

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
	params.dy = del.at(3, 0);

	params.theta = atan2(params.b, params.a);
	params.lambda = sqrt(pow(params.a, 2) + pow(params.b, 2));
}

void SimilarityTransform2D::similarityA() {
	A.resize(2 * num_points, 4);
	A.clear();

	for (unsigned int i = 0; i < num_points; i++) {
		A[2 * i][0] = coords_from[i].x;
		A[2 * i][1] = -coords_from[i].y;
		A[2 * i][2] = 1.0;

		A[2 * i + 1][0] = coords_from[i].y;
		A[2 * i + 1][1] = coords_from[i].x;
		A[2 * i + 1][3] = 1.0;
	}
}

Matrix SimilarityTransform2D::getA() {
	return A;
}

Matrix SimilarityTransform2D::getDelta() {
	return del;
}

Matrix SimilarityTransform2D::getResiduals() {
	return v;
}

vector<Point2D> SimilarityTransform2D::getFromCoords() {
	return coords_from;
}

vector<Point2D> SimilarityTransform2D::getToCoords() {
	return coords_from;
}

SimilarityParams SimilarityTransform2D::getParams() {
	return params;
}