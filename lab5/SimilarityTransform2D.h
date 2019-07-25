#pragma once

#include "LeastSquares.h"
#include "Point.h"

struct SimilarityParams {
	double a, b;
	double dx, dy;
	double lambda, theta;
};

class SimilarityTransform2D {
public:
	/** SimilarityTransform2D
	 * the constructor of this class; will compute the similarity transform
	 * with the given coordinates
	 * 
	 * @param coords_from - the starting coordinates for the transformation
	 * @param coords_to	  - the expected final coordinates after the transformation
	 */
	SimilarityTransform2D(vector<Point2D> &coords_from, vector<Point2D> &coords_to);

	/** setCoordinates
	 * sets the coordinates to their respective values and computes the 
	 * similarity transform
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
	
	SimilarityParams getParams();
private:
	Matrix A, del, v;

	unsigned int num_points;
	vector<Point2D> coords_from; // considered as KNOWNS in the LS adjustment
	vector<Point2D> coords_to;	 // considered as OBSERVATIONS in the LS adjustment

	SimilarityParams params;

	/** computeTransformation
	 * Computes all parameters in a linear similarity transformation
	 */
	void computeTransformation();

	/** similarityA
	 * Determines an N-by-4 design matrix for a linear similarity transformation
	 */
	void similarityA();
};

