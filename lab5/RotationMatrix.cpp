#include "RotationMatrix.h"

RotationMatrix::RotationMatrix() {
	resize(3, 3);
	clear();
}

RotationMatrix::RotationMatrix(double omega, double phi, double kappa) : RotationMatrix() {
	rotate(omega, phi, kappa);
}

void RotationMatrix::rotate(double rx, double ry, double rz) {
	Angles angle(rx, ry, rz);
	rotate(angle);
}

void RotationMatrix::rotate(Angles angle) {
	clear();
	rotateAboutX(angle.omega);
	rotateAboutY(angle.phi);
	rotateAboutZ(angle.kappa);
}

void RotationMatrix::rotateAboutX(double rx) {
	Matrix rot(3, 3);
	rot[0][0] = 1;
	rot[1][1] = cos(rx);
	rot[1][2] = sin(rx);
	rot[2][1] = -rot[1][2];
	rot[2][2] = rot[1][1];

	(*this) = rot * (*this);
}

void RotationMatrix::rotateAboutY(double ry) {
	Matrix rot(3, 3);
	rot[0][0] = cos(ry);
	rot[0][2] = -sin(ry);
	rot[1][1] = 1;
	rot[2][0] = -rot[0][2];
	rot[2][2] = rot[0][0];

	(*this) = rot * (*this);
}

void RotationMatrix::rotateAboutZ(double rz) {
	Matrix rot(3, 3);
	rot[0][0] = cos(rz);
	rot[0][1] = sin(rz);
	rot[1][0] = -rot[0][1];
	rot[1][1] = rot[0][0];
	rot[2][2] = 1;

	(*this) = rot * (*this);
}

RotationMatrix RotationMatrix::trans() {
	RotationMatrix temp;
	for (unsigned int i = 0; i < 3; i++)
		for (unsigned int j = 0; j < 3; j++)
			temp.at(j, i) = this->at(i, j);
	return temp;
}

void RotationMatrix::clear() {
	Matrix::clear();
	for (unsigned int i = 0; i < 3; i++) {
		this->at(i, i) = 1;
	}
}

double RotationMatrix::getOmega() {
	return atan2(-this->at(2, 1), this->at(2, 2));
}

double RotationMatrix::getPhi() {
	return asin(this->at(2, 0));
}

double RotationMatrix::getKappa() {
	return atan2(-this->at(1, 0), this->at(0, 0));
}

Angles RotationMatrix::getAngles() {
	return Angles(getOmega(), getPhi(), getKappa());
}

const RotationMatrix& RotationMatrix::operator= (const Matrix& mat1) {
	if (this == &mat1)
		return *this;

	resize(3, 3);

	for (unsigned int i = 0; i < 3; i++)
		for (unsigned int j = 0; j < 3; j++)
			this->at(i, j) = mat1.at(i, j);

	return *this;
}

const RotationMatrix operator* (const RotationMatrix& mat1, const RotationMatrix& mat2) {

	RotationMatrix temp;
	temp.at(0, 0) = 0;
	temp.at(1, 1) = 0;
	temp.at(2, 2) = 0;
	for (unsigned int i = 0; i < mat1.getrows(); i++)
		for (unsigned int j = 0; j < mat2.getcols(); j++)
			for (unsigned int k = 0; k < mat1.getcols(); k++)
				temp[i][j] += mat1[i][k] * mat2[k][j];

	return temp;
}