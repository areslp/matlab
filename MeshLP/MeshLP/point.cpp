#include "point.h"

#include <cmath>
#include <iostream>

using namespace std;

double Point::epsilon = 0.00001;

Point::Point(double x, double y, double z) :
x(x), y(y), z(z)
{
	
}

Point::Point(const Point& o) {
	x = o.x;
	y = o.y;
	z = o.z;
}



Point& Point::operator=(const Point& rhs) {
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;

	return *this;
}


void Point::normalize_point() {
	double norm = sqrt(x*x + y*y + z*z);

	if (norm != 0.0) {
		x /= norm;
		y /= norm;
		z /= norm;
	} else {
		x = y = z = 0.0;
	}
}


Point cross_product(Point a, Point b, Point c) {
	double u[3], v[3];

	u[0] = b.x - a.x;
	u[1] = b.y - a.y;
	u[2] = b.z - a.z;

	v[0] = c.x - a.x;
	v[1] = c.y - a.y;
	v[2] = c.z - a.z;

	double x,y,z;

	x = u[1]*v[2] - u[2]*v[1];
	y = u[2]*v[0] - u[0]*v[2];
	z = u[0]*v[1] - u[1]*v[0];

	return Point(x, y, z);
}

double dot_product(Point a, Point b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}


Point centroid(vector<Point> verts) {
	int size = verts.size();
	int i;

	double x = 0, y = 0, z = 0;
	for(i=0; i<size; i++) {
		x += verts[i].x / double(size);
		y += verts[i].y / double(size);
		z += verts[i].z / double(size);
	}

	return Point(x, y, z);
}


double triangleArea(const Point& a, const Point& b, const Point& c) {
	//actually returning Double the area
	return (b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y);
}


Point barycentricCoords(const Point& a, const Point& b, const Point& c,
								const Point& p) {

	//one needs to be careful about orientation here
	
	double areaABC = triangleArea(a,b,c);

	double areaPBC = triangleArea(p,b,c);
	
	double areaPCA = triangleArea(p,c,a);

	double alpha = areaPBC/areaABC;
	double beta = areaPCA/areaABC;
	double gamma = 1 - (alpha+beta);
	
	return Point(alpha, beta, gamma);
}



