 
#include<iostream>
#include<string>
#include<fstream>

#include "LeastSquares.h"
#include "Point.h"
#include "RotationMatrix.h"
#include "Matrix.h"
#include "SimilarityTransform2D.h"
#include "Resection.h"
#include "Lab5.h"

int main() {

	// Check if using test data from lecture slides
	bool testing = isTesting();

	vector<Point2D> coords_image = read2DPoints("Enter the image points filename : ");
	vector<Point3D> coords_object = read3DPoints("Enter the control points filename: ");

	string tolerance_filename;
	cout << "Enter the tolerances filename: ";
	getline(cin, tolerance_filename);

	Matrix tolerances;
	tolerances.read(tolerance_filename.c_str());

	const double H = 1845.45;
	double c = 153.358;

	if (testing)
		c = 152.15;

	SimilarityTransform2D similarity(coords_image, decreaseDimension(coords_object));
	SimilarityParams params = similarity.getParams();

	Resection resection(coords_object, coords_image, c);
	resection.computeResection(Point3D(params.dx, params.dy, H), Angles(0, 0, params.theta), tolerances);

	printStatistics(resection);
	printParameters(resection);
	
	std::cout << "Done!" << std::endl;
	std::cin.get();
	return 0;
}