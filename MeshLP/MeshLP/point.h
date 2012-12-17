#ifndef POINT_H
#define POINT_H

//points

#include <vector>

using namespace std;


class Point {
	friend Point cross_product(Point a, Point b, Point c);
	friend double dot_product(Point a, Point b);
	friend Point centroid(vector<Point> verts);
	
	public:
		double x,y,z;

		Point(double x, double y, double z);

		Point(const Point& o);
		Point& operator=(const Point& rhs);

		void normalize_point();
	
		bool operator==(Point &p) { 
			return ((x <= p.x + epsilon) &&  (x >= p.x - epsilon) &&
					  (y <= p.y + epsilon) &&  (y >= p.y - epsilon) &&
					  (z <= p.z + epsilon) &&  (z >= p.z - epsilon)
					 );
		}

		bool operator<(const Point& p) const {
			return ((x < p.x) || (x == p.x && y < p.y));
		}

	private:
		static double epsilon;
};

double triangleArea(const Point& a, const Point& b, const Point& c);
Point barycentricCoords(const Point& a, const Point& b, const Point& c,
								const Point& p);

#endif










