#pragma once

#include "Matrix.h"

struct Angles {
	double omega; // rotation about x-axis [rad]
	double phi;   // rotation about y-axis [rad]
	double kappa; // rotation about z-axis [rad]

	Angles() : omega(0), phi(0), kappa(0) {}
	Angles(double _omega, double _phi, double _kappa) :
		omega(_omega), phi(_phi), kappa(_kappa) {}
};

class RotationMatrix : public Matrix {
public:
	/** RotationMatrix
	 * the constructor of this class
	 * 
	 * @param omega - rotation about the x-axis [rad]
	 * @param phi	- rotation about the y-axis [rad]
	 * @param kappa - rotation about the z-axis [rad]
	 */
	RotationMatrix();
	RotationMatrix(double omega, double phi, double kappa);

	/** rotate
	 * rotate the matrix about all three (3) axis
	 * 
	 * R = R3(kappa) * R2(phi) * R1(omega)
	 * 
	 * @param rx	- rotation about the x-axis [rad]
	 * @param ry	- rotation about the y-axis [rad]
	 * @param rz	- rotation about the z-axis [rad]
	 * @param angle - containing all three rotation angles [rad]
	 */
	void rotate(double rx, double ry, double rz);
	void rotate(Angles angle);

	void rotateAboutX(double rx); // rotate current matrix about its x-axis
	void rotateAboutY(double ry); // rotate current matrix about its y-axis
	void rotateAboutZ(double rz); // rotate current matrix about its z-axis

	RotationMatrix trans(); // returns the transpose (inverse) rotation matrix
	void clear(); // clear the rotation matrix

	double getOmega();  // returns x-axis rotation [rad]
	double getPhi();    // returns y-axis rotation [rad]
	double getKappa();  // returns z-axis rotation [rad]
	Angles getAngles(); // returns all angles [rad]

	const RotationMatrix& operator= (const Matrix& mat1);
	friend const RotationMatrix operator* (const RotationMatrix& mat1, const RotationMatrix& mat2);
private:
};