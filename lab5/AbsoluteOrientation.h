#pragma once

#include "RotationMatrix.h"
#include "LeastSquares.h"
#include "Point.h"

class AbsoluteOrientation {
public:
	/** AbsoluteOrientation
	 * the constructor of this class; sets up the object and model coordinates, respectively
	 * 
	 * @param coords_object - the object coordinates of the selected points
	 * @param coords_model	- the model coordinates of the selected points
	 */
	AbsoluteOrientation(const vector<Point3D> &coords_object, const vector<Point3D> &coords_model);

	/** computeOrientation
	 * Computes all parameters in the absolute orientation adjustment
	 * 
	 * @param _T	  - the point of expansion for the translation vector
	 * @param _ang	  - the point of expansion for the rotation angles
	 * @param _lambda - the point of expansion for the scale
	 */
	void computeOrientation(const Point3D &_T, const Angles &_ang, double _lambda);

	Matrix getA();
	Matrix getMisclosure();
	Matrix getDelta();
	Matrix getResiduals();

	vector<Point3D> getObjectCoords();
	vector<Point3D> getModelCoords();

	RotationMatrix getM();
	Point3D getT();
	double getScale();
private:
	Matrix A, w, del, v;

	unsigned int num_points;
	vector<Point3D> coords_object;
	vector<Point3D> coords_model;

	RotationMatrix M; // rotates from model to object space
	Point3D T;		  // translation vector from object to model space
	double lambda;	  // scale from model to object space

	/** absoluteA
	 * Determines a 3n-by-7 design matrix for an absolute orientation LS
	 * 
	 * observations: [Xo, Yo, Zo]
	 * unknowns: [omega, phi, kappa, s, Tx, Ty, Tz]
	 * 
	 * parametric function:
	 * [Xo]		[m11 m12 m13] [Xm]	 [Tx]
	 * [Yo] = s [m21 m22 m33] [Ym] + [Ty]
	 * [Zo]		[m31 m32 m33] [Zm]	 [Tz]
	 * 
	 * @return - the 3n-by-7 design matrix
	 */
	Matrix absoluteA();

	/** absoluteCond
	 * Determines a 3n-by-1 vector estimating the parametric function for absolute orientation
	 *
	 * observations: [Xo, Yo, Zo]
	 * unknowns: [omega, phi, kappa, s, Tx, Ty, Tz]
	 *
	 * parametric function:
	 * [Xo]		[m11 m12 m13] [Xm]	 [Tx]
	 * [Yo] = s [m21 m22 m33] [Ym] + [Ty]
	 * [Zo]		[m31 m32 m33] [Zm]	 [Tz]
	 * 
	 * @return - the 3n-by-1 vector of estimates
	 */
	Matrix absoluteCond();

	/** computeRowX
	 * updates the i-th row of the design matrix, A, for X with all its partial derivative values
	 * 
	 * A(Xi) = d{omega, phi, kappa, lambda, Tx, Ty, Tz}
	 * 
	 * @param A	   - the current design matrix being calculated (manipulated)
	 * @param row  - the X row of which to manipulate in the A matrix
	 * @param pm   - the model point coordinate corresponding to the current i-th section of the A matrix
	 * @param sinv - a vector containing all the sine values in order of sin{omega, phi, kappa}
	 * @param cosv - a vector containing all the cosine values in order of cos{omega, phi, kappa}
	 */
	void computeRowX(Matrix &A, int row, const Point3D pm, const vector<double> &sinv, const vector<double> &cosv);

	/** computeRowY
	 * updates the i-th row of the design matrix, A, for Y with all its partial derivative values
	 *
	 * A(Yi) = d{omega, phi, kappa, lambda, Tx, Ty, Tz}
	 *
	 * @param A	   - the current design matrix being calculated (manipulated)
	 * @param row  - the Y row of which to manipulate in the A matrix
	 * @param pm   - the model point coordinate corresponding to the current i-th section of the A matrix
	 * @param sinv - a vector containing all the sine values in order of sin{omega, phi, kappa}
	 * @param cosv - a vector containing all the cosine values in order of cos{omega, phi, kappa}
	 */
	void computeRowY(Matrix &A, int row, const Point3D pm, const vector<double> &sinv, const vector<double> &cosv);

	/** computeRowY
	 * updates the i-th row of the design matrix, A, for Z with all its partial derivative values
	 *
	 * A(Zi) = d{omega, phi, kappa, lambda, Tx, Ty, Tz}
	 *
	 * @param A	   - the current design matrix being calculated (manipulated)
	 * @param row  - the Z row of which to manipulate in the A matrix
	 * @param pm   - the model point coordinate corresponding to the current i-th section of the A matrix
	 * @param sinv - a vector containing all the sine values in order of sin{omega, phi, kappa}
	 * @param cosv - a vector containing all the cosine values in order of cos{omega, phi, kappa}
	 */
	void computeRowZ(Matrix &A, int row, const Point3D pm, const vector<double> &sinv, const vector<double> &cosv);
};