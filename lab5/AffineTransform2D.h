#pragma once

#include "LeastSquares.h"
#include "Point.h"

struct AffineParams {
	double a, b, c, d;
	double dx, dy;
	double sx, sy;
	double theta, delta;
};

class AffineTransform2D {
public:
	/** AffineTransform2D
	 * the constructor of this class; will compute the affine transform
	 * with the given coordinates
	 *
	 * @param coords_from - the starting coordinates for the transformation
	 * @param coords_to	  - the expected final coordinates after the transformation
	 */
	AffineTransform2D(vector<Point2D> &coords_from, vector<Point2D> &coords_to);

	/** setCoordinates
	 * sets the coordinates to their respective values and computes the
	 * affine transform
	 *
	 * @param coords_from - the starting coordinates for the transformation
	 * @param coords_to	  - the expected final coordinates after the transformation
	 */
	void setCoordinates(vector<Point2D> &coords_from, vector<Point2D> &coords_to);

	Matrix getA();
	Matrix getDelta();
	Matrix getResiduals();

	vector<Point2D> getFromCoords();
	vector<Point2D> getToCoords();

	AffineParams getParams();
private:
	Matrix A, del, v;

	unsigned int num_points;
	vector<Point2D> coords_from; // considered as KNOWNS in the LS adjustment
	vector<Point2D> coords_to;	 // considered as OBSERVATIONS in the LS adjustment

	AffineParams params;

	/** computeTransformation
	 * Computes all parameters in a linear affine transformation
	 */
	void computeTransformation();

	/** affine_A
	 * Determines an N-by-6 design matrix for a linear affine transformation
	 */
	void affineA();
};