
#pragma once

#include "RotationMatrix.h"
#include "LeastSquares.h"
#include "Point.h"

class RelativeOrientation {
public:
	/** RelativeOrientation
	 * the constructor of this class; converts the set of 2D coordinates into 3D points
	 * by applying the focal length; the LEFT IMAGE defines the datum
	 * 
	 * (x, y) -> (x, y, -c)
	 * 
	 * @param coords_left  - the 2D set of coordinates for the left image
	 * @param coords_right - the 2D set of coordinates for the right image
	 * @param c			   - the focal length of the images
	 */
	RelativeOrientation(const vector<Point2D> &coords_left, const vector<Point2D> &coords_right, double c);

	/** computeOrientation
	 * computes all parameters in the relative orientation adjustment
	 * 
	 * @param _B   - the point of expansion for the base vector
	 * @param _ang - the point of expansion for the rotation angles
	 */
	void computeOrientation(const Point3D &_B, const Angles &_ang);

	Matrix getA();

	vector<Point3D> getLeftCoords();
	vector<Point3D> getRightCoords();

	vector<Point3D> getModelCoords();
	vector<Point3D> getParallax();

	RotationMatrix getM();
	Point3D getB();
	double getFocalLength();

private:
	Matrix A;

	unsigned int num_points;
	vector<Point3D> coords_left; // constructs reference frame
	vector<Point3D> coords_right;

	vector<Point3D> coords_model;
	vector<Point3D> parallax;

	RotationMatrix M; // rotates from model to image space
	Point3D B;
	double c; // focal length

	/** coplanarityA
	 * Computes a n-by-5 design matrix for a relative orientaion transformation
	 *
	 * coplanarity condition: B . (Vl x (M' * Vr)) = 0
	 */
	void coplanarityA();

	/** coplanarityCond
	 * Computes a n-by-1 vector of the coplanarity condition equations
	 *
	 * coplanarity condition: B . (Vl x (M' * Vr)) = 0
	 *
	 * @return - the desired output equation vector
	 */
	Matrix coplanarityCond();

	/** computeModelSpace
	 * Determines the final model space coordinates as well as using the RO parallax using
	 * the space intersected coordinates of each image
	 */
	void computeModelSpace();

	/** determinant
	 * Computes the determinant of a strictly 3-by-3 matrix
	 *
	 * @param mat - the matrix in which we want the determinant of
	 *
	 * @return	  - the determinant of the matrix
	 */
	double determinant(Matrix &mat);
};


