#include "RelativeOrientation.h"

RelativeOrientation::RelativeOrientation(const vector<Point2D> &coords_left, const vector<Point2D> &coords_right, double c) {
	this->num_points = coords_left.size();
	this->coords_left = increaseDimension(coords_left, -c);
	this->coords_right = increaseDimension(coords_right, -c);
	this->c = c;
}

void RelativeOrientation::computeOrientation(const Point3D &_B, const Angles &_ang) {
	B = _B;
	Angles ang = _ang;
	M.rotate(ang.omega, ang.phi, ang.kappa);

	Matrix w, del;

	double threshold = 1e-6;
	do {
		coplanarityA();
		w = coplanarityCond();
		del = delta(A, w);

		B.y += del[0][0];
		B.z += del[1][0];
		ang.omega += del[2][0];
		ang.phi += del[3][0];
		ang.kappa += del[4][0];

		M.rotate(ang);
	} while (del.maxAbsElem() > threshold);

	computeModelSpace();
}

void RelativeOrientation::coplanarityA() {

	A.resize(num_points, 5);

	Point3D pl, pr;
	Matrix partials(3, 3);

	double omega = M.getOmega();
	double phi = M.getPhi();
	double kappa = M.getKappa();

	for (unsigned int i = 0; i < num_points; i++) {
		pl = coords_left[i];
		pr = coords_right[i];

		pr.rotateBy(M.trans());

		// differential wrt By
		A[i][0] = pr.x * pl.z - pl.x * pr.z;

		// differential wrt Bz
		A[i][1] = pl.x * pr.y - pr.x * pl.y;

		// rotation differentials (setup the matrix)
		partials[0][0] = B.x;  partials[0][1] = B.y;  partials[0][2] = B.z;
		partials[1][0] = pl.x; partials[1][1] = pl.y; partials[1][2] = pl.z;

		// differential wrt omega
		partials[2][0] = 0;
		partials[2][1] = -pr.z;
		partials[2][2] = pr.y;

		A[i][2] = determinant(partials);

		// differential wrt phi
		partials[2][0] = -pr.y * sin(omega) + pr.z * cos(omega);
		partials[2][1] = pr.x * sin(omega);
		partials[2][2] = -pr.x * cos(omega);

		A[i][3] = determinant(partials);

		// differential wrt kappa
		partials[2][0] = -pr.y * cos(omega) * cos(phi) - pr.z * sin(omega) * cos(phi);
		partials[2][1] = pr.x * cos(omega) * cos(phi) - pr.z * sin(phi);
		partials[2][2] = pr.x * sin(omega) * cos(phi) + pr.y * sin(phi);

		A[i][4] = determinant(partials);
	}
}

Matrix RelativeOrientation::coplanarityCond() {

	Matrix cond;
	cond.resize(num_points, 1);

	Point3D pl, pr;
	Matrix condition(3, 3);

	for (unsigned int i = 0; i < num_points; i++) {
		pl = coords_left[i];
		pr = coords_right[i];

		pr.rotateBy(M.trans());

		condition[0][0] = B.x;  condition[0][1] = B.y;  condition[0][2] = B.z;
		condition[1][0] = pl.x; condition[1][1] = pl.y; condition[1][2] = pl.z;
		condition[2][0] = pr.x; condition[2][1] = pr.y; condition[2][2] = pr.z;

		cond[i][0] = determinant(condition);
	}

	return cond;
}

void RelativeOrientation::computeModelSpace() {

	coords_model.clear();
	parallax.clear();

	Point3D pl, pr;

	for (unsigned int i = 0; i < num_points; i++) {
		pl = coords_left[i];
		pr = coords_right[i];

		pr.rotateBy(M.trans());

		// explicitly solve for coefficients because Matrix inv() does not support
		// non positive definite matrices
		double lambda = (B.x * pr.z - B.z * pr.x) / (pl.x * pr.z - pl.z * pr.x);
		double mu = (B.x * pl.z - B.z * pl.x) / (pl.x * pr.z - pl.z * pr.x);
		
		pl.scaleBy(lambda);

		pr.scaleBy(mu);
		pr.translateBy(B);

		Point3D mp; // model point
		mp.x = pl.x;
		mp.y = (pl.y + pr.y) / 2;
		mp.z = pl.z;
		coords_model.push_back(mp);

		Point3D pp = pr.difference(pl); // parallax point
		parallax.push_back(pp);
	}
}

Matrix RelativeOrientation::getA() {
	return A;
}

vector<Point3D> RelativeOrientation::getLeftCoords() {
	return coords_left;
}

vector<Point3D> RelativeOrientation::getRightCoords() {
	return coords_right;
}

vector<Point3D> RelativeOrientation::getModelCoords() {
	return coords_model;
}

vector<Point3D> RelativeOrientation::getParallax() {
	return parallax;
}

RotationMatrix RelativeOrientation::getM() {
	return M;
}

Point3D RelativeOrientation::getB() {
	return B;
}

double RelativeOrientation::getFocalLength() {
	return c;
}

double RelativeOrientation::determinant(Matrix &mat) {
	if (mat.getcols() != 3 || mat.getrows() != 3) {
		std::cout << "Incorrect sized matrix. Exiting..." << std::endl;
		exit(1);
	}

	return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
			mat[1][0] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
			mat[2][0] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);
}