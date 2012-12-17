#ifndef OFFOBJ_H
#define OFFOBJ_H


#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "point.h"

using namespace std;








//class for colors
class RGBA {
	public:
		RGBA() {	
		}

		RGBA(const RGBA& o) {
			r = o.r; g = o.g; b = o.b; a = o.a;
		}

		RGBA& operator=(const RGBA& o) {
			r = o.r; g = o.g; b = o.b; a = o.a;
			return *this;
		}

		bool operator==(RGBA& o) {
			return (r == o.r) && (g == o.g) && (b == o.b) && (a == o.a);
		}

		double r,g,b,a;
};



class OffObj {
	public:
		OffObj();
		OffObj(const OffObj& o);

		OffObj& operator=(const OffObj& o);
		bool operator==(OffObj& o);

		void writeOff(ofstream & theFile);
		void updateBBox(double x, double y, double z);
		void orient_facets();
		void calculate_facet_normals();
		void calculate_vert_normals();
		
		vector<Point> vertices;
		vector<vector<int> > facets;
		vector<RGBA> vert_colors;
		vector<RGBA> facet_colors;		
		vector<vector<int> > topo_verts;
		vector<vector<int> > topo_facets;
		vector<Point> facet_normals;
		vector<Point> vert_normals;

		bool is_coff;

		Point bbmin, bbmax;
};




class OffFileReader {
	public:
		OffFileReader(string& filename);
		~OffFileReader();
		
		string trim(string& in);
		void split(string& in, const string& delims, vector<string>& tokens);
		bool get_next_line(vector<string>& tokens);
		bool get_next_off_object(OffObj& o);

		bool bad;

	private:
		ifstream input;
};



#endif




