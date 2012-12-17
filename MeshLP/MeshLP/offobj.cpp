#include "offobj.h"

#include <deque>
#include <algorithm>
using namespace std;


OffObj::OffObj() : 
	bbmin(9999999,9999999,9999999), 	bbmax(-9999999,-9999999,-9999999), is_coff(false)
{

}

OffObj::OffObj(const OffObj& o) :
	bbmin(9999999,9999999,9999999), 	bbmax(-9999999,-9999999,-9999999) 
{
	vertices = o.vertices;
	topo_verts = o.topo_verts;
	facets = o.facets;
	vert_colors = o.vert_colors;
	facet_colors = o.facet_colors;	
	facet_normals = o.facet_normals;
	vert_normals = o.vert_normals;
	is_coff = o.is_coff;

	bbmin = o.bbmin;
	bbmax = o.bbmax;
}

OffObj& OffObj::operator=(const OffObj& o) {
	vertices = o.vertices;
	topo_verts = o.topo_verts;
	facets = o.facets;
	vert_colors = o.vert_colors;
	facet_colors = o.facet_colors;
	facet_normals = o.facet_normals;
	vert_normals = o.vert_normals;
	is_coff = o.is_coff;

	bbmin = o.bbmin;
	bbmax = o.bbmax;

	return *this;
}

bool OffObj::operator==(OffObj& o) {
	return (vertices.size() == o.vertices.size()) &&
		(facets.size() == o.facets.size());
}





void OffObj::writeOff(ofstream & theFile) {
	if (this->is_coff) {
		theFile << "COFF ";
	} else {
		theFile << "OFF ";
	}
	theFile << vertices.size() << " ";
	theFile << facets.size() << " ";
	theFile << 0 << endl;

	int i,j;
	for(i=0; i<vertices.size(); i++) {
		Point p = vertices[i];
		theFile << p.x << " " << p.y << " " << p.z << endl;
		if (this->is_coff) {
			RGBA color = vert_colors[i];
			theFile << " " << color.r;
			theFile << " " << color.g;
			theFile << " " << color.b;
			theFile << " " << color.a;
		}
	}

	for(i=0; i<facets.size(); i++) {
		vector<int> facet = facets[i];

		theFile << facet.size();
		for(j=0; j<facet.size(); j++) {
			theFile << " " << facet[j];
		}

		RGBA color = facet_colors[i];
		if (color.r < 0.0) {
			//no colors, don't output them
		} else {
			theFile << " " << color.r;
			theFile << " " << color.g;
			theFile << " " << color.b;
			theFile << " " << color.a;
		}

		theFile << endl;
	}
}

void OffObj::updateBBox(double x, double y, double z) {
	if (x < bbmin.x) {
		bbmin.x = x;
	}
	if (y < bbmin.y) {
		bbmin.y = y;
	}
	if (z < bbmin.z) {
		bbmin.z = z;
	}
	if (x > bbmax.x) {
		bbmax.x = x;
	}
	if (y > bbmax.y) {
		bbmax.y = y;
	}
	if (z > bbmax.z) {
		bbmax.z = z;
	}
}




void OffObj::orient_facets() {
	vector<bool> marked_facets;
	vector<vector<int> > topo_facets;
	int i,j;

	for(i=0; i<facets.size(); i++) {
		marked_facets.push_back(false);

		vector<int> facet = facets[i];
		vector<int> topo_facet;

		for(j=0; j<facet.size(); j++) {
			int v1, v2;
			v1 = facet[j];
			v2 = facet[(j+1)%facet.size()];

			//find two intersections
			int l,m;
			int adjacent_facet = 0;
			bool found = false;
			vector<int> topo_vert1 = topo_verts[v1];
			vector<int> topo_vert2 = topo_verts[v2];
			for(l=0; l<topo_vert1.size() && !found; l++) {
				for(m=0; m<topo_vert2.size() && !found; m++) {
					if (topo_vert1[l] == topo_vert2[m] && topo_vert1[l] != i) {
						adjacent_facet = topo_vert1[l];
						found = true;
					}
				}
			}

			if (found) {
				topo_facet.push_back(adjacent_facet);
			}
		}

		topo_facets.push_back(topo_facet);
	}

	cerr << "finishing building facet neighbor list" << endl;



	//first is the facet,
	//second is the facet to orient with
	deque<pair<int, int> > facet_queue;
	for(i=0; i<topo_facets[0].size(); i++) {
		facet_queue.push_back(make_pair<int, int>(topo_facets[0][i], 0));
		marked_facets[topo_facets[0][i]] = true;
	}
	marked_facets[0] = true;

	int flip_count = 0;
	int total_count = 1;
	while(facet_queue.size() != 0) {
		pair<int, int> facet_info = facet_queue[0];
		facet_queue.pop_front();

		//push adjacent facets onto the queue;
		for(i=0; i<topo_facets[facet_info.first].size(); i++) {
			int neighbor = topo_facets[facet_info.first][i];
			if (!marked_facets[neighbor]) {
				facet_queue.push_back(make_pair<int, int>(neighbor, facet_info.first));
				marked_facets[neighbor] = true;
			}
		}		
		vector<int> facet1 = facets[facet_info.first];
		vector<int> facet2 = facets[facet_info.second];

		//compare orientation of facet1 with facet2
		bool flip_facet = false;
		if (facet1[0] == facet2[0]) {
			if (facet1[1] == facet2[1]) {
				//flip
				flip_facet = true;
			} else if (facet1[2] == facet2[2]) {
				//flip
				flip_facet = true;
			} else {
				//keep same
			}
		} else if (facet1[0] == facet2[1]) {
			if (facet1[1] == facet2[2]) {
				//flip
				flip_facet = true;
			} else if (facet1[2] == facet2[0]) {
				//flip
				flip_facet = true;
			} else {
				//keep same
			}
		} else if (facet1[0] == facet2[2]) {
			if (facet1[1] == facet2[0]) {
				//flip
				flip_facet = true;
			} else if (facet1[2] == facet2[1]) {
				//flip
				flip_facet = true;
			} else {
				//keep same
			}
		} else {
			//facet1[1] and facet1[2] are the shared verts
			if (facet1[1] == facet2[0]) {
				if (facet1[2] == facet2[1]) {
					//flip
					flip_facet = true;
				} else {
					//keep same
				}
			} else if (facet1[1] == facet2[1]) {
				if (facet1[2] == facet2[2]) {
					//flip
					flip_facet = true;
				} else {
					//keep same
				}
			} else if (facet1[1] == facet2[2]) {
				if (facet1[2] == facet2[0]) {
					//flip
					flip_facet = true;
				} else {
					//keep same
				}
			}
		}

		if (flip_facet) {
			int tmp = facet1[0];
			facet1[0] = facet1[1];
			facet1[1] = tmp;

			facets[facet_info.first] = facet1;

			//flip normal
			Point normal = facet_normals[facet_info.first];

			normal.x = -normal.x;
			normal.y = -normal.y;
			normal.z = -normal.z;

			facet_normals[facet_info.first] = normal;

			flip_count++;
		}
		total_count++;
	}

	cerr << "flipped: " << flip_count << " of " << total_count << " facets." << endl;

}


void OffObj::calculate_facet_normals() {
	int i;
	for(i=0; i<facets.size(); i++) {
		vector<int> facet = facets[i];

		//calculate normal for facet
		//use first 3 points
		Point normal(0, 0, 0);

		if (facet.size() > 2) {
			Point	pt0 = vertices[facet[0]];
			Point pt1 = vertices[facet[1]];
			Point pt2 = vertices[facet[2]];

			normal = cross_product(pt0, pt1, pt2);

			normal.normalize_point();
		} else {
			//an edge or a point
		}

		facet_normals[i] = normal;
	}
}


//assumes facet normals exist!
void OffObj::calculate_vert_normals() {
	//calculate vert normals
	int i,j;
	for(i=0; i<vertices.size(); i++) {
		vector<Point> normal_list;
		for(j=0; j<topo_verts[i].size(); j++) {
			normal_list.push_back(facet_normals[topo_verts[i][j]]);
		}

		Point normal = centroid(normal_list);
		normal.normalize_point();

		vert_normals.push_back(normal);
	}
}























OffFileReader::OffFileReader(string& filename) {
	input.open(filename.c_str());

	if (!input) {
		cerr << "File '" << filename << " does not exist" << endl;
		bad = true;
	} else {
		bad = false;
	}
}

OffFileReader::~OffFileReader() {
	input.close();
}

string OffFileReader::trim(string& in) {
	//remove leading whitespace

	string ret = in;

	string::iterator p = ret.begin();

	while (p != ret.end() && *p == ' ') {
		p = ret.erase(p);
	}	

	return ret;	
}


void OffFileReader::split(string& in, const string& delims, vector<string>& tokens) {
	tokens.clear();

	string::iterator p1 = in.begin();
	string::iterator p2 = in.begin();

	while(p2 != in.end()) {
		if (delims.find(*p2) != string::npos) {
			//there's a delimiter here
			if (p1 == p2) {
				//ignore the character
				p1++;
			} else {
				//p1 to p2 is a substring
				string s(p1, p2);
				tokens.push_back(s);
				//advance p1
				p1 = p2;
				p1++;
			}				
		}

		p2++;
	}

	if (p1 != p2) {
		//p1 to p2 is a substring
		string s(p1, p2);
		tokens.push_back(s);
	}
}

bool OffFileReader::get_next_line(vector<string>& tokens) {
	string s1, s2;
	bool done = false;
	bool end = false;

	while (!done && !end) {
		getline(input, s1, '\n');
		s2 = trim(s1);

		if (input.eof()) {
			end = true;
		} else if (s2.length() == 0) {
			//an empty line, skip
		} else if (s2[0] == '#') {
			//a comment line, skip
		} else {
			split(s2, "{ }\t", tokens);
			//if (tokens[0] != "appearance") {
				if (tokens.size() != 0) {
					done = true;
				}
			//} else {
				//skip appearance
		//	}
		}			
	}


	if (end) {
		return false;
	} else {
		return true;
	}		
}

bool OffFileReader::get_next_off_object(OffObj& o) {
	vector<string> tokens;
	bool found_off = false;


	o.facets.clear();
	o.facet_colors.clear();
	o.vertices.clear();
	o.topo_verts.clear();
	o.facet_normals.clear();
	o.vert_normals.clear();


	cerr.flags(ios::fixed);
	cerr.precision(2);


	int num_verts=0, num_facets=0, num_outliers=0;
	int i=0,j=0;


	while(!found_off && get_next_line(tokens)) {
		if (tokens[0] == "OFF") {
			o.is_coff = false;
			//found an off
			found_off = true;
			if (tokens.size() == 4) {
				num_verts = atoi(tokens[1].c_str());
				num_facets = atoi(tokens[2].c_str());
				num_outliers = atoi(tokens[3].c_str());
			} else if (tokens.size() == 1) {
				//get the next line
				get_next_line(tokens);
				if (tokens.size() == 3) {
					num_verts = atoi(tokens[0].c_str());
					num_facets = atoi(tokens[1].c_str());
					num_outliers = atoi(tokens[2].c_str());
				}
			} else {
				//format error
				cerr << "OFF file format, missing number of verts/facets/outliers" << endl;
				return false;
			}


			//get verts

			for(i=0; i<num_verts; i++) {
				//get geom vertices
				double x = 0, y = 0, z = 0;
				get_next_line(tokens);
				x = atof(tokens[0].c_str());
				y = atof(tokens[1].c_str());
				z = atof(tokens[2].c_str());
				o.vertices.push_back(Point(x, y, z));
				vector<int> topo_vert;
				o.topo_verts.push_back(topo_vert);

				o.updateBBox(x, y, z);

				if (i%500 == 0) {
					cerr << "Vertices " << 100*double(i)/double(num_verts) << "%\r";
				}
			}
			cerr << "Vertices " << 100*double(i)/double(num_verts) << "%" << endl;
			cerr << num_verts << " vertices loaded." << endl;



			//get facets

			for(i=0; i<num_facets; i++) {
				//get geom vertices
				int verts_in_facet = 0;
				get_next_line(tokens);
				verts_in_facet = atoi(tokens[0].c_str());

				vector<int> facet;
				for(j=0; j<verts_in_facet; j++) {
					int which_vert = 0;
					which_vert = atoi(tokens[j+1].c_str());
					facet.push_back(which_vert);			
				}

				o.facets.push_back(facet);

				for(j=0; j<facet.size(); j++) {
					o.topo_verts[facet[j]].push_back(o.facets.size()-1);
				}

				//get rgba for facet
				double r = 0, g = 0, b = 0, a = 0;

				Point	p = o.vertices[facet[0]];
				if (tokens.size() == 1+verts_in_facet) {
					//some default color
					/*	
						float delta[3] = {0};
						delta[0] = bbmax.x - bbmin.x;
						delta[1] = bbmax.y - bbmin.y;
						delta[2] = bbmax.z - bbmin.z;
						r = p.x/delta[0]/2.0 + 1.0;
						g = p.y/delta[1]/2.0 + 1.0;
						b = p.z/delta[2]/2.0 + 1.0;
						a = 1.0f;
					 */	
					//r = 0.45f; g = 0.5f; b = 0.9f; a = 1.0f;
					r = 0.35f; g = 0.8f; b = 0.9f; a = 1.0f;
				} else {
					//read them
					if (tokens.size() == 5+verts_in_facet) {
						r = atof(tokens[verts_in_facet+1].c_str());
						g = atof(tokens[verts_in_facet+2].c_str());
						b = atof(tokens[verts_in_facet+3].c_str());
						a = atof(tokens[verts_in_facet+4].c_str());
					} else {
						//error
						cerr << "invalid facet specification" << endl;
						return false;
					}
				}

				RGBA facet_color;
				facet_color.r = r;
				facet_color.g = g;
				facet_color.b = b;
				facet_color.a = a;
				o.facet_colors.push_back(facet_color);

				//just give a default
				Point normal(0, 0, 0);
				o.facet_normals.push_back(normal);

				if (i%1000 == 0) {
					cerr << "Facets " << 100*double(i)/double(num_facets) << "%\r";
				}

			}
			cerr << "Facets " << 100*double(i)/double(num_facets) << "%" << endl;
			cerr << num_facets << " facets loaded." << endl << endl;
		} else if (tokens[0] == "COFF") {
			//found an off
			o.is_coff = true;
			found_off = true;
			if (tokens.size() == 4) {
				num_verts = atoi(tokens[1].c_str());
				num_facets = atoi(tokens[2].c_str());
				num_outliers = atoi(tokens[3].c_str());
			} else if (tokens.size() == 1) {
				//get the next line
				get_next_line(tokens);
				if (tokens.size() == 3) {
					num_verts = atoi(tokens[0].c_str());
					num_facets = atoi(tokens[1].c_str());
					num_outliers = atoi(tokens[2].c_str());
				}
			} else {
				//format error
				cerr << "OFF file format, missing number of verts/facets/outliers" << endl;
				return false;
			}


			//get verts

			for(i=0; i<num_verts; i++) {
				//get geom vertices
				double x = 0, y = 0, z = 0;
				get_next_line(tokens);
				x = atof(tokens[0].c_str());
				y = atof(tokens[1].c_str());
				z = atof(tokens[2].c_str());
				o.vertices.push_back(Point(x, y, z));
				vector<int> topo_vert;


				//get rgba for facet
				double r = 0, g = 0, b = 0, a = 0;

			
				if (tokens.size() == 3) {
					r = 0.45f; g = 0.5f; b = 0.9f; a = 1.0f;
				} else if (tokens.size() == 7) {
					r = atof(tokens[3].c_str());
					g = atof(tokens[4].c_str());
					b = atof(tokens[5].c_str());
					a = atof(tokens[6].c_str());
				} else {
					cerr << "invalid vertex specified" << endl;
					return false;
				}

				RGBA vert_color;
				vert_color.r = r;
				vert_color.g = g;
				vert_color.b = b;
				vert_color.a = a;
				o.vert_colors.push_back(vert_color);


				o.topo_verts.push_back(topo_vert);

				o.updateBBox(x, y, z);

				if (i%500 == 0) {
					cerr << "Vertices " << 100*double(i)/double(num_verts) << "%\r";
				}
			}
			cerr << "Vertices " << 100*double(i)/double(num_verts) << "%" << endl;
			cerr << num_verts << " vertices loaded." << endl;



			//get facets

			for(i=0; i<num_facets; i++) {
				//get geom vertices
				int verts_in_facet = 0;
				get_next_line(tokens);
				verts_in_facet = atoi(tokens[0].c_str());

				vector<int> facet;
				for(j=0; j<verts_in_facet; j++) {
					int which_vert = 0;
					which_vert = atoi(tokens[j+1].c_str());
					facet.push_back(which_vert);			
				}

				o.facets.push_back(facet);

				for(j=0; j<facet.size(); j++) {
					o.topo_verts[facet[j]].push_back(o.facets.size()-1);
				}

				//get rgba for facet
				double r = 0, g = 0, b = 0, a = 0;

				Point	p = o.vertices[facet[0]];
				if (tokens.size() == 1+verts_in_facet) {
					//some default color
					/*	
						float delta[3] = {0};
						delta[0] = bbmax.x - bbmin.x;
						delta[1] = bbmax.y - bbmin.y;
						delta[2] = bbmax.z - bbmin.z;
						r = p.x/delta[0]/2.0 + 1.0;
						g = p.y/delta[1]/2.0 + 1.0;
						b = p.z/delta[2]/2.0 + 1.0;
						a = 1.0f;
					 */	
					r = 0.45f; g = 0.5f; b = 0.9f; a = 1.0f;
				} else {
					//read them
					if (tokens.size() == 5+verts_in_facet) {
						r = atof(tokens[verts_in_facet+1].c_str());
						g = atof(tokens[verts_in_facet+2].c_str());
						b = atof(tokens[verts_in_facet+3].c_str());
						a = atof(tokens[verts_in_facet+4].c_str());
					} else {
						//error
						cerr << "invalid facet specification" << endl;
						return false;
					}
				}

				RGBA facet_color;
				facet_color.r = r;
				facet_color.g = g;
				facet_color.b = b;
				facet_color.a = a;
				o.facet_colors.push_back(facet_color);

				//just give a default
				Point normal(0, 0, 0);
				o.facet_normals.push_back(normal);

				if (i%1000 == 0) {
					cerr << "Facets " << 100*double(i)/double(num_facets) << "%\r";
				}

			}
			cerr << "Facets " << 100*double(i)/double(num_facets) << "%" << endl;
			cerr << num_facets << " facets loaded." << endl << endl;
		}
	}


	return found_off;
}





