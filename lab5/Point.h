#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Matrix.h"
#include "RotationMatrix.h"

struct GeodeticPoint {
	string id;
	double lat, lon, h;

	GeodeticPoint() : id(), lat(0), lon(0), h(0) {}
	GeodeticPoint(double _lat, double _lon, double _h) : id(), lat(_lat), lon(_lon), h(_h) {}
	GeodeticPoint(string _id, double _lat, double _lon, double _h) : id(_id), lat(_lat), lon(_lon), h(_h) {}
};

struct Point2D {
	string id;
	double x, y;

	Point2D() : id(), x(0), y(0) {}
	Point2D(double _x, double _y) : id(), x(_x), y(_y) {}
	Point2D(string _id, double _x, double _y) : id(_id), x(_x), y(_y) {}
};

struct Point3D {
	string id;
	double x, y, z;

	Point3D() : id(), x(0), y(0), z(0) {}
	Point3D(double _x, double _y, double _z) : id(), x(_x), y(_y), z(_z) {}
	Point3D(string _id, double _x, double _y, double _z) : id(_id), x(_x), y(_y), z(_z) {}

	/** scaleBy
	 * scales the 3D point by a given float value
	 * 
	 * @param scale - the scale factor in which to multiply each component {x, y, z} by
	 */
	void scaleBy(double scale) {
		x *= scale;
		y *= scale;
		z *= scale;
	}

	/** translateBy
	 * translates the 3D point by a given translation 3D point
	 * 
	 * @param T - the translation point in which to add to each component {x, y, z}
	 */
	void translateBy(const Point3D &T) {
		x += T.x;
		y += T.y;
		z += T.z;
	}

	/** rotateBy
	 * rotates the 3D point by a given rotation matrix
	 * 
	 * @param R - the rotation matrix in which to rotate the entire 3D point by
	 */
	void rotateBy(const RotationMatrix &R) {
		Matrix vec(3, 1);
		vec[0][0] = x;
		vec[1][0] = y;
		vec[2][0] = z;

		vec = R * vec;

		x = vec[0][0];
		y = vec[1][0];
		z = vec[2][0];
	}

	/** transform
	 * combined scale, rotation, and translation transormation to the 3D point
	 * 
	 * @param scale - the scale factor in which to multiply each component {x, y, z} by
	 * @param R	    - the rotation matrix in which to rotate the entire 3D point by
	 * @param T		- the translation point in which to add to each component {x, y, z}
	 */
	void transform(double scale, const RotationMatrix &R, const Point3D &T) {
		rotateBy(R);
		scaleBy(scale);
		translateBy(T);
	}

	// difference = {x, y, z} - p.{x, y, z}
	Point3D difference(const Point3D p) const {
		return Point3D(x - p.x, y - p.y, z - p.z);
	}
};

/** increaseDimension
 * increases a set of 2D points into the third dimension. If no z-value is provided
 * then all points will lie in the plane z=0
 * 
 * @param points2D - the set of 2D points to increase the dimension of
 * @param value	   - the z-value specifying the height of all points
 * 
 * @return		   - a set of 3D points residing in the plane z=const.
 */
vector<Point3D> increaseDimension(const vector<Point2D> &points2D);
vector<Point3D> increaseDimension(const vector<Point2D> &points2D, double value);

/** decreaseDimension
 * projects a set of 3D points down into the plane z=0 by removing their z-component
 * 
 * @param points3D - the set of 3D points to decrease the dimension of
 *
 * @return - a set of 2D points residing in the 2D plane (z=0)
 */
vector<Point2D> decreaseDimension(const vector<Point3D> &points3D);

/** transformPoints
 * does a complete 3D transformation a set of 3D points
 * 
 * @param points - the desired 3D points to be transformed
 * @param lambda - the scale of the transformation
 * @param R		 - the 3D rotation matrix of the transformation
 * @param T		 - the 3D translation vector of the transformaiton
 *
 * @return		 - the fully transformed set of 3D coordinates
 */
vector<Point3D> transformPoints(const vector<Point3D> &points, double lambda, const RotationMatrix &R, const Point3D &T);

/** convertToVector
 * converts any points to a vector of alternating [x,y,{z}] pairs
 * 
 * @param  - the set of points in which to convert to a vector
 * 
 * @return - the vector containing alternating [x,y,{z}] pairs
 */
Matrix convertToVector(const vector<Point2D> &points);
Matrix convertToVector(const vector<Point3D> &points);

/** findPoint
 * finds the FIRST OCCURANCE of a {geodetic, 2D, 3D} point with a given id
 * 
 * @param points - the set of points of which to search. Can be either geodetic, 2D, or 3D points
 * @param id	 - the id of the desired point
 * 
 * @return		 - the first occurring point or a newly generated blank point if not found
 */
GeodeticPoint findPoint(const vector<GeodeticPoint> &points, string id);
Point2D findPoint(const vector<Point2D> &points, string id);
Point3D findPoint(const vector<Point3D> &points, string id);

/** findPoints
 * finds ALL OCCURANCES of a {geodetic, 2D, 3D} points with a given id
 * 
 * @param points - the set of points of which to search. Can be either geodetic, 2D, or 3D points
 * @param id	 - the id of the desired points
 * 
 * @return		 - all found occurances of the points with a given id
 */
vector<GeodeticPoint> findPoints(const vector<GeodeticPoint> &points, string id);
vector<Point2D> findPoints(const vector<Point2D> &points, string id);
vector<Point3D> findPoints(const vector<Point3D> &points, string id);

/** readGeodeticPoints
 * reads in a set of geodetic points from a text file
 * 
 * DOEST NOT REQUIRE A HEADER CONTAINING THE NUMBER OF POINTS
 * 
 * @param prompt - a simple prompt to ask the user for the required filename
 */
vector<GeodeticPoint> readGeodeticPoints(string prompt);

/** read2DPoints
 * reads in a set of 2D points from a text file
 * 
 * DOEST NOT REQUIRE A HEADER CONTAINING THE NUMBER OF POINTS
 * 
 * @param prompt - a simple prompt to ask the user for the required filename
 */
vector<Point2D> read2DPoints(string prompt);

/** read3DPoints
* reads in a set of 3D points from a text file
*
* DOEST NOT REQUIRE A HEADER CONTAINING THE NUMBER OF POINTS
*
* @param prompt - a simple prompt to ask the user for the required filename
*/
vector<Point3D> read3DPoints(string prompt);

/** printPoints
 * prints either a given vector of points (2D or 3D) to a specified file
 * 
 * @param points   - a vector of 2D or 3D point objects which you wish to print to a file
 * @param filename - the filename of of which to print the points to
 * @param prec	   - the decimal precision for the output points (optional)
 * @param width	   - the allocated width of each printed number (optional)
 */
void printPoints(const vector<Point2D> &points, const char *filename, int prec = 3, int width = 12);
void printPoints(const vector<Point3D> &points, const char *filename, int prec = 3, int width = 12);
