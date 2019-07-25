#include "Point.h"

vector<Point3D> increaseDimension(const vector<Point2D> &points2D) {
	return increaseDimension(points2D, 0);
}

vector<Point3D> increaseDimension(const vector<Point2D> &points2D, double value) {
	vector<Point3D> points3D;

	for (Point2D p : points2D) {
		points3D.push_back(Point3D(p.x, p.y, value));
	}
	
	return points3D;
}

vector<Point2D> decreaseDimension(const vector<Point3D> &points3D) {
	vector<Point2D> points2D;

	for (Point3D p : points3D) {
		points2D.push_back(Point2D(p.x, p.y));
	}

	return points2D;
}

vector<Point3D> transformPoints(const vector<Point3D> &points, double scale, const RotationMatrix &R, const Point3D &T) {
	vector<Point3D> transformed;

	for (Point3D p : points) {
		p.transform(scale, R, T);
		transformed.push_back(p);
	}

	return transformed;
}

Matrix convertToVector(const vector<Point2D> &points) {
	Matrix vect(points.size() * 2, 1);

	for (unsigned int i = 0; i < points.size(); i++) {
		vect.at(2 * i, 0) = points[i].x;
		vect.at(2 * i + 1, 0) = points[i].y;
	}

	return vect;
}

Matrix convertToVector(const vector<Point3D> &points) {
	Matrix vect(points.size() * 3, 1);

	for (unsigned int i = 0; i < points.size(); i++) {
		vect.at(3 * i, 0) = points[i].x;
		vect.at(3 * i + 1, 0) = points[i].y;
		vect.at(3 * i + 2, 0) = points[i].z;
	}

	return vect;
}

GeodeticPoint findPoint(const vector<GeodeticPoint> &points, string id) {
	for (GeodeticPoint gp : points) {
		if (gp.id == id)
			return gp;
	}
	return GeodeticPoint();
}

Point2D findPoint(const vector<Point2D> &points, string id) {
	for (Point2D p : points) {
		if (p.id == id)
			return p;
	}
	return Point2D();
}

Point3D findPoint(const vector<Point3D> &points, string id) {
	for (Point3D p : points) {
		if (p.id == id)
			return p;
	}
	return Point3D();
}

vector<GeodeticPoint> findPoints(const vector<GeodeticPoint> &points, string id) {
	vector<GeodeticPoint> found;
	
	for (GeodeticPoint gp : points) {
		if (gp.id == id)
			found.push_back(gp);
	}
	
	return found;
}

vector<Point2D> findPoints(const vector<Point2D> &points, string id) {
	vector<Point2D> found;

	for (Point2D p : points) {
		if (p.id == id)
			found.push_back(p);
	}

	return found;
}

vector<Point3D> findPoints(const vector<Point3D> &points, string id) {
	vector<Point3D> found;

	for (Point3D p : points) {
		if (p.id == id)
			found.push_back(p);
	}

	return found;
}

vector<GeodeticPoint> readGeodeticPoints(string prompt) {
	string filename;
	std::cout << prompt;
	getline(std::cin, filename);

	ifstream infile;
	infile.open(filename, ifstream::in);

	if (!infile.is_open())
	{
		std::cout << "Error: read - Error opening file" << std::endl;
		exit(1);
	}

	vector<GeodeticPoint> points;

	GeodeticPoint p;
	while (infile >> p.id >> p.lat >> p.lon >> p.h) {
		points.push_back(p);
	}

	infile.close();

	return points;
}

vector<Point2D> read2DPoints(string prompt) {
	string filename;
	std::cout << prompt;
	std::getline(std::cin, filename);

	ifstream infile;
	infile.open(filename, ifstream::in);

	if (!infile.is_open())
	{
		std::cout << "Error: read - Error opening file" << std::endl;
		exit(1);
	}

	vector<Point2D> points;

	Point2D p;
	while (infile >> p.id >> p.x >> p.y) {
		points.push_back(p);
	}

	infile.close();

	return points;
}

vector<Point3D> read3DPoints(string prompt) {
	string filename;
	std::cout << prompt;
	std::getline(std::cin, filename);

	ifstream infile;
	infile.open(filename, ifstream::in);

	if (!infile.is_open())
	{
		std::cout << "Error: read - Error opening file" << std::endl;
		exit(1);
	}

	vector<Point3D> points;

	Point3D p;
	while (infile >> p.id >> p.x >> p.y >> p.z) {
		points.push_back(p);
	}

	infile.close();

	return points;
}

void printPoints(const vector<Point2D> &points, const char *filename, int prec, int width) {
	ofstream outfile;
	outfile.open(filename);

	if (!outfile.is_open())
	{
		cout << "Error: print - Error opening file" << endl;
		exit(1);
	}

	outfile.precision(prec);
	outfile.setf(std::ios::showpoint);
	outfile.setf(std::ios::fixed, std::ios::floatfield);

	for (Point2D p : points) {
		outfile << p.id << "\t";
		outfile << setw(width) << p.x << "\t" << p.y << endl;
	}
}

void printPoints(const vector<Point3D> &points, const char *filename, int prec, int width) {
	ofstream outfile;
	outfile.open(filename);

	if (!outfile.is_open())
	{
		cout << "Error: print - Error opening file" << endl;
		exit(1);
	}

	outfile.precision(prec);
	outfile.setf(std::ios::showpoint);
	outfile.setf(std::ios::fixed, std::ios::floatfield);

	for (Point3D p : points) {
		outfile << p.id << "\t";
		outfile << setw(width) << p.x << "\t" << p.y << "\t" << p.z << endl;
	}
}