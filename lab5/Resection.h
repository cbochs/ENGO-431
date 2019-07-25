#pragma once

#include "LeastSquares.h"
#include "Point.h"
#include "RotationMatrix.h"
#include "Matrix.h"

class Resection {
public:
	Resection(const vector<Point3D> &coords_object, const vector<Point2D> &coords_image, double c);

	/** computeOrientation
	 * Computes all parameters in the resection
	 *
	 * @param _T	  - the point of expansion for the translation vector
	 * @param _ang	  - the point of expansion for the rotation angles
	 */
	void computeResection(const Point3D &_T, const Angles &_ang, const Matrix &tolerances);

	Matrix getA();
	Matrix getMisclosure();
	Matrix getDelta();
	Matrix getResiduals();

	vector<Point3D> getObjectCoords();
	vector<Point2D> getImageCoords();

	RotationMatrix getM();
	Point3D getT();
	double getFocalLength();
private:
	Matrix A, w, del, v;

	unsigned int num_points;
	vector<Point3D> coords_object;
	vector<Point2D> coords_image;

	RotationMatrix M; // rotates from model to object space
	Point3D T;		  // translation vector from object to model space
	double c;		  // focal length

	/** resectionA
	 * Determines a 2n-by-7 design matrix for resection
	 *
	 * observations: [Xi, Yi]
	 * unknowns: [Tx, Ty, Tz, omega, phi, kappa]
	 *
	 * collinearity condition:
	 * [Xi]			[m11 m12 m13] [X - Tx]	 [U]
	 * [Yi] = (1/s) [m21 m22 m33] [Y - Ty] = [V]
	 * [Zi]			[m31 m32 m33] [Z - Tz]   [W]
 	 *
	 * @return - the 2n-by-7 design matrix
	*/
	Matrix resectionA();

	/** resectionCond
	 * Determines a 2n-by-1 vector estimating the parametric function for resection
	 *
	 * observations: [Xi, Yi]
	 * unknowns: [Tx, Ty, Tz, omega, phi, kappa]
	 *
	 * collinearity condition:
	 * [Xi]			[m11 m12 m13] [X - Xc]	 [U]
	 * [Yi] = (1/s) [m21 m22 m33] [Y - Yc] = [V]
	 * [Zi]			[m31 m32 m33] [Z - Zc]	 [W]
	 *
	 * @return - the 2n-by-1 vector of estimates
	 */
	Matrix resectionCond();

	/** computeRowX
	 * updates the i-th row of the design matrix, A, for X with all its partial derivative values
	 *
	 * A(Xi) = d{Tx, Ty, Tz, omega, phi, kappa}
	 *
	 * @param A	   - the current design matrix being calculated (manipulated)
	 * @param row  - the X row of which to manipulate in the A matrix
	 * @param op   - the object point coordinate corresponding to the current i-th section of the A matrix
	 * @param sinv - a vector containing all the sine values in order of sin{omega, phi, kappa}
	 * @param cosv - a vector containing all the cosine values in order of cos{omega, phi, kappa}
	 */
	void computeRowX(Matrix &A, int row, const Point3D &op, const vector<double> &sinv, const vector<double> &cosv);

	/** computeRowY
	 * updates the i-th row of the design matrix, A, for Y with all its partial derivative values
	 *
	 * A(Yi) = d{Tx, Ty, Tz, omega, phi, kappa}
	 *
	 * @param A	   - the current design matrix being calculated (manipulated)
	 * @param row  - the Y row of which to manipulate in the A matrix
	 * @param op   - the object point coordinate corresponding to the current i-th section of the A matrix
	 * @param sinv - a vector containing all the sine values in order of sin{omega, phi, kappa}
	 * @param cosv - a vector containing all the cosine values in order of cos{omega, phi, kappa}
	 */
	void computeRowY(Matrix &A, int row, const Point3D &op, const vector<double> &sinv, const vector<double> &cosv);

	/** U
	 * computes the U value of the coplanarity condition for a given point
	 * 
	 * @param point - some object point
	 * 
	 * @return		- the output U value
	 */
	double U(const Point3D &point);

	/** V
	* computes the V value of the coplanarity condition for a given point
	*
	* @param point - some object point
	*
	* @return		- the output V value
	*/
	double V(const Point3D &point);

	/** W
	* computes the W value of the coplanarity condition for a given point
	*
	* @param point - some object point
	*
	* @return		- the output W value
	*/
	double W(const Point3D &point);
};