
#include "Lab5.h"

void printStatistics(Resection &resection) {
	Matrix v = resection.getResiduals();
	Matrix A = resection.getA();
	Matrix corr = correlation(A);
	Matrix Cx = unknownCovariance(A, (0.0119*0.0119));
	Matrix Cv = residualCovariance(A);
	Matrix redund(Cv.getrows(), 1);
	for (unsigned int i = 0; i < redund.getrows(); i++) {
		redund.at(i, 0) = Cv.at(i, i);
	}

	v.print(3, 12, "Residuals");
	corr.print(3, 12, "Correlation");
	Cx.print(3, 12, "Variance-Covariance");
	redund.print(3, 12, "Redundancy");

	v.print("residuals.txt", 3);
	corr.print("correlation.txt", 3);
	Cx.print("covariance.txt", 9);
	redund.print("redundancy.txt", 3);
}

void printParameters(Resection &resection) {
	RotationMatrix M = resection.getM();
	Point3D T = resection.getT();

	vector<Point3D> params;
	params.push_back(T);
	params.push_back(Point3D(rad2deg(M.getOmega()), rad2deg(M.getPhi()), rad2deg(M.getKappa())));

	printPoints(params, "params.txt", 6);
}

bool isTesting() {
	string ans;
	std::cout << "Using test data (Y/N)? ";
	std::getline(std::cin, ans);
	if (ans[0] == 'Y')
		return true;
	else
		return false;
}

double rad2deg(double rad) {
	return rad * (180 / pi);
}

double deg2rad(double deg) {
	return deg * (pi / 180);
}